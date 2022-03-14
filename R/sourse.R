
pauc<-function(auc,n = 100,n.plus=0.5){
  
  np<- n.plus
  if(n.plus<1){
    if(n.plus<0){
      np<- round(0.5*n)
    }else{
      np<-round(n.plus*n)
    }
  }
  nm<- n-np
  
  return(min(c(2*pwilcox(nm*np*(1-auc), nm, np),1)))
}


##############################################################
### computing p-value of AUC #################################
##############################################################

compute.p.auc<-function(auc, n, n.plus=0.5,alternative="two.sided"){
  p<-NA
  if(is.null(auc)||is.na(auc)){return(p)}
  if(alternative=="greater"){
    p<-pauc(auc,n,n.plus)
  }
  
  if(alternative=="less"){
    p<-1-pauc(1-auc,n,n.plus)
  }
  if(alternative=="two.sided"){
    p1<-pauc(auc,n,n.plus)
    p2<-1-pauc(1-auc,n,n.plus)
    p<-p1+p2
  }
  p
}


compute.auc<-function(x,labs, pos=NULL,alternative="two.sided",significance=T){
  #######x: measurements, labs: class labels, pos: positive class, alternative: greater, less, two.sided, significance: whether to compute the p-value, to.expect: whether to compute the expected AUC-value
  labs<- as.factor(labs)
  if(is.null(pos)){
    lv<-levels(as.factor(labs))
    pos<-lv[1]
  }
  is.pos<-TRUE
  nona.xi<-which(!is.na(x))
  
  #####Computing AUC-value
  if(length(nona.xi)==0){
    res<-data.frame(AUC=NA,Positive.Correlation=is.pos)
  }else{	
    nona.level<-unique(labs[nona.xi])
    if(length(nona.level)==1){
      is.pos<-pos==nona.level
      n.pos<-length(which(labs==nona.level))
      res<-data.frame(AUC=NA,Positive.Correlation=is.pos)
    }else{  
      
      pred <- ROCR::prediction(as.numeric(x), labs)
      aucv <- ROCR::performance(pred,"tpr", "fpr", measure = "auc")
      aucval <- attr(aucv,"y.values")[[1]]
      
      ##### Judging the positive correlation          
      if (aucval<0.5){
        aucval <- 1-aucval
        is.pos <- FALSE  
####if auc-value >=0.5, positive correlation is true; otherwise, positive correlation is false
      }
      res<-data.frame(AUC=aucval, Positive.Correlation=is.pos)
    }
  }
  
  if(!significance){
    return(res)
  }
  
  ##### Counting the sample sizes: the whole sample size and the positive sample size
  x1<-x[nona.xi]
  labs1<-labs[nona.xi]
  n.x<-length(x1)
  pos.objs<-which(labs1==pos)
  
  if(is.pos){
    pos.x1<-x1[pos.objs]
  }else{
    pos.x1<-x1[-pos.objs]
  }  
  n.plus<-length(pos.x1)
  
  ###### Computing p-value of AUC   
  pval<-compute.p.auc(aucval,n.x,n.plus,alternative=alternative)
  res<-data.frame(res, P.value=pval)

  return(res)
}

computeAUCs <- function(dattable,significance=T, alternative="two.sided"){####dattable ist character und der letzte Spalt ist diagnose .###pos: die selbst difinierte Plus-Klasse 
  labs<- as.factor(dattable[,ncol(dattable)]) 
 
  number_of_cores <- detectCores()-1
  cl <- makePSOCKcluster(number_of_cores) #not to overload your computer
  registerDoParallel(cl)
  getDoParWorkers()
  
  xx<- data.frame(dattable[,-ncol(dattable)])
  
  pfms <- foreach(i = 1:ncol(xx), .export=c("compute.auc","compute.p.auc" ,"pauc"), 
                  .packages=c('ROCR'))%dopar% {
                    
                    result <- compute.auc(xx[,i],labs)
                    
                    return(result)
                  }
  
  stopCluster(cl)
  pfm<- do.call(rbind.data.frame, pfms)
  
  return(pfm)
  
  
}

generate.random.bm <- function(n,n.plus=0.5,score.range=NULL,m)
{
  np <- n.plus
  
  labs <- factor(c(rep("+",np),rep("-",n-np)))
  
  scores <- matrix(rep(1,nrow = n, ncol=m))
  if (is.null(score.range)){
    scores <- matrix(runif(n*m), nrow=n, ncol=m)
  }else{
    scores <- matrix(sample(score.range,n*m,replace=T,nrow=n , ncol=m))
  }
  return(data.frame(scores,labs))
}
######################################################
### selecting top n features: ########################
### weight means the selection probability ###########
######################################################
#######################################################################################
### selecting features according to AUC-value or the AUC and interclass correlation ###
### output: the column numbers of features and weight: the selection probability    ### 
#######################################################################################
###############################################################
######### select  a data include only the top features ########
###########  ##################################  ##############
select.n.features<-function(pfm,n.feats){
  n.feats<-min(n.feats,length(pfm))
  n.feats <- as.numeric(n.feats)
  pfm.level<-unique(pfm)  
  m<-length(pfm.level)
  w<-rep(0,length(pfm))
  
  if(length(pfm)==m){
    selfeat<-order(pfm, decreasing = T)[1:n.feats]
    w[selfeat]<-1
    return(list(features=selfeat,weight=w))
  }
  
  selfeat<-rep(0,n.feats)
  k<-0
  pfm.level<-sort(pfm.level)
  while(k<n.feats){
    xi<-which(pfm==pfm.level[m])
    n.xi<-length(xi)
    w[xi]<-1
    res.selfeat<-n.feats-k
    if(n.xi<=res.selfeat){
      selfeat[c((k+1):(k+n.xi))]<-xi
      k<-k+n.xi
    }else{
      w[xi]<-res.selfeat/n.xi
      selfeat[c((k+1):n.feats)]<-xi[sample(n.xi,res.selfeat)]
      k<-n.feats
    }
    m<-m-1
  }
  return(list(features=selfeat,weight=w))
}

select.features<- function(dattable,method="auc",max.no.features=20,is.fast=T, k.MI=3){
  #new
  if(is.null(nrow(dattable))){
    return(NULL)
  }
  
  inddat <- 1:nrow(dattable)
  labs<-dattable[,ncol(dattable)]
  
  switch(method, 
         auc={aucs<-computeAUCs(dattable)
         vals<-aucs[,1]   
         }, 
         acc={ vals<-compute.accuracies(dattable, is.fast)     
         }, 
         entropy={vals<-compute.entropies(dattable,is.fast)*(-1)
         },
         kNN.MI={vals<-apply(dattable[,-ncol(dattable)], 2, kNN.MI, dattable[,ncol(dattable)],k=k.MI)
         }
  )
  
  sf<-select.n.features(vals, max.no.features)
  
  return(sf)
}

selectData<- function(data, No.features, method="auc",is.fast=T, k.MI=3){
  lv<-levels(as.factor(data[,ncol(data)]))
  pos<-lv[1]
  selectedFeature<- select.features(data,method=method,max.no.features=No.features,is.fast=T, k.MI=3)
  features<- selectedFeature$features
  
  
  selectDataF<- matrix(0, nrow(data), length(features))

  for(i in 1:length(features)){
    a<-data[,features[i]]
    
    selectDataF[,i]<-a
  }
  selectDataF<- as.data.frame(selectDataF)
  selectDataF$labels<- data[,ncol(data)]
  return(selectDataF)
}

##########################################################################################
###################### HAUCA Curve   ##########################################
#################################################################################


plotLines<-function(x, y, xlab=NULL, ylab=NULL, main=NULL, titles=NULL,  shape="vh", opacity=NULL) {
  if(is.vector(y)) y<-matrix(y, nrow=length(x))
  
  if(is.null(titles)) titles<-paste(colnames(y))
  if(is.null(opacity)) opacity<-rep(1, ncol(y))
  
  
  if(length(x)!=nrow(y) | ncol(y)!=length(titles) | ncol(y)!=length(opacity))
    stop("The length of 'x', the row number of 'y' must be identical, as well as
         the column number of 'y', the length of 'titles' and 'opacity' !")
  

  xa <- list(
    title = xlab,
    titlefont = FALSE
  )
  ya <- list(
    title = ylab,
    titlefont = FALSE
  )

  
  p<-plotly::plot_ly(x=x) %>% plotly::layout(title=main,xaxis = xa, yaxis = ya)
  
  
  for(i in 1:ncol(y)) {
    p<-plotly::add_lines(p, y=y[, i], line=list(shape=shape), opacity=opacity[i], name=titles[i])
  }
  return(p)
}

HAUCA.table<- function(data){
  
  a<- computeAUCs(data)

  correctd_P.value<- a[,3]* (ncol(data)-1)

  d<- data.frame(a,correctd_P.value)

  return(d)
}

HAUCA.curve<- function(data, method = "auc",   from.auc= 0.5, to.auc= 1, increment= 0.01, 
                       confidence.interval= 0.95, xlab="AUC values", ylab="Number of Features", main="HAUCA Curve", 
                       titles=NULL, shape="vh",opacity=NULL){
  
 
  if(from.auc<0 | from.auc>1 | to.auc<0 |to.auc>1){
    stop("AUC not between 0 &1 ")
  }
  
  levs<- levels(as.factor(data[,ncol(data)]))

  pos.labs<- subset(data, data[,ncol(data)] == levs[1])
  
  rand.Data<- generate.random.bm(nrow(data),nrow(pos.labs),score.range = NULL, ncol(data[,-ncol(data)]))
  
  seqs<-seq(from=from.auc, to=to.auc, by=increment)
  if (method == "auc"){
    ## real Data
    real_aucs<- computeAUCs(data)
    real.Data<-sapply(seqs, function(x) sum(real_aucs$AUC>x))
    
    
    ## random Data
    rand_aucs<- computeAUCs(rand.Data)
    random.Data<-sapply(seqs, function(x) sum(rand_aucs$AUC>x))
    
  }else{
    ## real Data
    real_entp<- compute.entropies(data)
    real.Data<-sapply(seqs, function(x) sum(real_entp$AUC<x))
    
    ## random Data
    rand_entp<- compute.entropies(rand.Data)
    random.Data<-sapply(seqs, function(x) sum(rand-entp$AUC<x))
  }
  ### Confidence Interval for AUC
  confedince.Interval<-c(1:length(seqs))
  for(i in 1: length(seqs)){
    
    p<- compute.p.auc(seqs[i], nrow(data),nrow(pos.labs))
    q<-qbinom(confidence.interval, ncol(data)-1, p)
    confedince.Interval[i]<- q
  }
  auc_frame<- data.frame(real.Data, random.Data, confedince.Interval)
  return(plotLines(seqs, auc_frame , xlab=xlab, ylab= ylab, main=main, titles=titles, shape= shape, opacity=opacity))
}
############################################################################
############### FUNCTIONS FOR COMBINING BIOMARKERS #######################
#####################################################################

##### compute covariance for selected features ###########
########  ###############################   ##############
difference.sd<- function(x,labs,method="sd"){
  labs<- as.factor(labs)
  dat<- data.frame(x,labs)
  
  levs <- levels(dat[,ncol(dat)])
  
  d1 <- subset(dat,dat[,ncol(dat)]==levs[1])
  d2 <- subset(dat,dat[,ncol(dat)]==levs[2])
  
  sd1<- sd(d1[,1])
  sd2<- sd(d2[,1])
  return(c(sd1,sd2))
}

difference.mean<- function(x,labs){
  labs<- as.factor(labs)
  dat<- data.frame(x,labs)
  dat[,ncol(dat)]<- as.factor(dat[,ncol(dat)])
  levs <- levels(dat[,ncol(dat)])
  
  d1 <- subset(dat,as.factor(dat[,ncol(dat)])==levs[1])
  d2 <- subset(dat,as.factor(dat[,ncol(dat)])==levs[2])
  
  mean1<- mean(d1[,1])
  mean2<- mean(d2[,1])
  
  return(abs(mean1-mean2))
}

difference.means <- function(dat,method="mean"){
  dat[,ncol(dat)]<- as.factor(dat[,ncol(dat)])
  levs <- levels(dat[,ncol(dat)])
  d1 <- subset(dat,dat[,ncol(dat)]==levs[1])
  d2 <- subset(dat,dat[,ncol(dat)]==levs[2])
  
  ### Here you could for instance use the median as an alternative to the mean.
  
  means1 <- apply(d1[,-ncol(d1)],2,method) 
  means2 <- apply(d2[,-ncol(d2)],2,method)
  
  return(abs(means1-means2))
}

#########################################################
##### Intra-class variance covariance matrices #########
#########################################################
cor.cov<- function(x,y, method="kendall"){
  if(method=="kendall"){
    k.cor<- cor(x,y,method="kendall")
    cov<- k.cor*sd(x)*sd(y)
    
  }
  if(method=="pearson"){
    p.cor<- cor(x,y,method="pearson")
    cov<- p.cor*sd(x)*sd(y)
    
  }
  if(method=="corrected.pearson"){
    k.cor<- cor(x,y,method="kendall")
    corct.p.cor<- sin(pi*k.cor/2)
    cov<- corct.p.cor*sd(x)*sd(y)
    
  }
  
  return(cov)
}

cov.mat<- function(data,method.cov="standard",method.cor= "corrected.pearson"){
  levs <- levels(data[,ncol(data)])
  d1 <- subset(data,data[,ncol(data)]==levs[1])
  d2 <- subset(data,data[,ncol(data)]==levs[2])
  
  if (method.cov== "covRob"){
    cov1 <- covRob(d1[,-ncol(d1)],corr = correlation)
    cov2 <- covRob(d2[,-ncol(d2)],corr = correlation)
    
  }
  if (method.cov== "standard"){
    mat.1<- matrix(0,ncol(data[,-ncol(data)]),ncol(data[,-ncol(data)]))
    colnames(mat.1)<- colnames(data[,-ncol(data)])
    rownames(mat.1)<- colnames(data[,-ncol(data)])
    
    mat.2<- matrix(0,ncol(data[,-ncol(data)]),ncol(data[,-ncol(data)]))
    colnames(mat.2)<- colnames(data[,-ncol(data)])
    rownames(mat.2)<- colnames(data[,-ncol(data)])
    
    if (method.cor!="No Correlation"){
      for(i in 1:ncol(data[,-ncol(data)])){
        for(j in 1: ncol(data[,-ncol(data)])){
          mat.1[i,j]<- cor.cov(d1[,i], d1[,j], method=method.cor)
          mat.2[i,j]<- cor.cov(d2[,i], d2[,j], method=method.cor)
          
        }
      }
      
    }
    if(method.cor=="No Correlation"){
      for(i in 1:ncol(data[,-ncol(data)])){
        mat.1[i,i]<- var(d1[,i])
        mat.2[i,i]<- var(d2[,i])
        
      }
    }
  }
  else{
    stop("the method is either 'standard' or 'covRob' from robust package ")
  }
  return(list(sigma0=mat.1,sigma1=mat.2))
}

auc.combi <- function(dat,method.mean="mean",method.cov="standard",method.cor= "pearson"){
  delta.mu <- difference.means(dat,method=method.mean)
  covam <- cov.mat(dat,method.cov=method.cov, method.cor= method.cor)
  
  return(pnorm(sqrt(t(delta.mu)%*%make.positive.definite(solve(covam$sigma0+covam$sigma1))%*%delta.mu)))
}
#######################################################
##### compute combined auc for selected features ###########
########  #################################   ##############
#

combinedAUCs<- function(data,no.features,method.combine="Theoretical",method.mean="mean",method.cov="standard", neg.class=NULL, method.cor= "corrected.pearson"){ 
  ## data: data frame with the outcome in the last column
  ## no.features: number of feature which are wanted to be combined 
  ## method.mean: either mean or median 
  ## method.cov: either Standard devision or covRob
  ## neg.class: The class which considerd negative in outcome, used with the combination based on real dat
  ## correlation: either taking the correlation into consideration or not
  ## method.cor: method of correlation, either pearson, kendall, or corrected.pearson according to kendall correlation 
  data[,ncol(data)]<- as.factor(data[,ncol(data)])
  combauc<-c(1:no.features)
  labels<- as.factor(data[,ncol(data)])
  
  dat<- selectData(data,no.features)
  
  ## correct the directions of Biomarkers
  
  aucs<- computeAUCs(dat)
  
  for(i in 1: nrow(aucs)){
    
    ifelse(aucs[i,2]== T, dat[,i]<- dat[,i], dat[,i]<- -1 * dat[,i])
  }
  
  if (method.combine== "Theoretical"){
    
    ## first Bio
    md<- difference.mean(dat[,1],dat[,ncol(dat)])
    sd<- difference.sd(dat[,1],dat[,ncol(dat)])
    auc1<- pnorm(md/sqrt(sd[1]^2+ + sd[2]^2))
    
    combauc[1]<- auc1
    
    ## combining 
    for(i in 2:no.features){
      comb.Bio<- dat[,1:i]
      frame<- data.frame(comb.Bio, labels)
      
      comb_auc<- auc.combi(frame, method.mean, method.cov,  method.cor)
      combauc[i]<-comb_auc
    }
  }
  
  if (method.combine== "RealData"){
    
    ## normalized the features by mean and sd
    daN<- subset(dat, dat[,ncol(dat)] == neg.class)

    if(nrow(daN)== 0){
      stop("please check the negative class if it is written correctly and please without quotation marks")
    }
    me<- apply(daN[,-ncol(daN)],2,mean)
    sigma<- apply(daN[,-ncol(daN)],2,sd)
    
    norm.data<- matrix(0, nrow(dat),ncol(dat)-1)
    colnames(norm.data)<- colnames(dat[,-ncol(dat)])
    rownames(norm.data)<- rownames(dat[,-ncol(dat)])
    
    for (i in 1:ncol(dat)-1){
      d<- as.numeric(dat[,i])
      da<- (d - me[i]) / sigma[i]
      
      norm.data[,i]<- da
    }
    norm.data<- data.frame(norm.data, dat$labels)
    
    
    #### 
    ## compute the AUCs when combining the features
    comb<- rep(0,no.features)
    
    auc1<- compute.auc(norm.data[,1], norm.data[,ncol(norm.data)])[1]
    comb[1]<- auc1
    
    for (i in 2:no.features){
      newbio<- apply(norm.data[,1:i], 1, sum)
      aucComb<- compute.auc(newbio, norm.data[,ncol(norm.data)])[1]
      comb[i]<- aucComb
      
    }
    
    combauc<- as.numeric(comb)
  }
  return(combauc)
}


##########################################################################################
theor_combinedAUCs_sameOrder<- function(data, no.features,method.mean="mean",
                            method.cov="standard",method.cor= "corrected.pearson"){
  combauc<-c(1:no.features)
  data[,ncol(data)]<- as.factor(data[,ncol(data)])
  
  ## correct the directions of Biomarkers
  
  aucs<- computeAUCs(data)
  
  for(i in 1: nrow(aucs)){
    
    ifelse(aucs[i,2]== T, data[,i]<- data[,i], data[,i]<- -1 * data[,i])
  }
  ## first Bio
  md<- difference.mean(data[,1],data[,ncol(data)])
  sd<- difference.sd(data[,1],data[,ncol(data)])
  auc1<- pnorm(md/sqrt(sd[1]^2 + sd[2]^2))
  
  combauc[1]<- auc1
  
  ## start combining
  for(i in 2:no.features){
    comb.Bio<- data[,1:i]
    frame<- data.frame(comb.Bio, data[,ncol(data)])
    comb.auc<- auc.combi(frame, method.mean, method.cov,  method.cor)
    
    combauc[i]<-comb.auc
  }
  return(combauc)
}

boot.comb.normalized<- function(dat,no.features,neg.class){
  dat[,ncol(dat)]<- as.factor(dat[,ncol(dat)])
  aucs<- computeAUCs(dat)
  
  for(i in 1: nrow(aucs)){
    
    ifelse(aucs[i,2]== T, dat[,i]<- dat[,i], dat[,i]<- -1 * dat[,i])
  }
  
  
  #### 
  ## normalized the features by mean and sd
  daN<- subset(dat, dat[,ncol(dat)] == neg.class)
  me<- apply(daN[,-ncol(daN)],2,mean)
  sigma<- apply(daN[,-ncol(daN)],2,sd)
  
  
  norm.data<- matrix(0, nrow(dat),ncol(dat)-1)
  colnames(norm.data)<- colnames(dat[,-ncol(dat)])
  rownames(norm.data)<- rownames(dat[,-ncol(dat)])
  
  for (i in 1:ncol(dat)-1){
    d<- as.numeric(dat[,i])
    da<- (d - me[i]) / sigma[i]
    
    norm.data[,i]<- da
  }
  norm.data<- data.frame(norm.data, dat[,ncol(dat)])
  

  ## compute the AUCs when combining the features
  comb<- rep(0,no.features)
  auc1<- compute.auc(norm.data[,1], norm.data[,ncol(norm.data)])[1]
  comb[1]<- auc1
  
  for (i in 2:no.features){
    newbio<- apply(norm.data[,1:i], 1, sum)
    aucComb<- compute.auc(newbio, norm.data[,ncol(norm.data)])[1]
    comb[i]<- aucComb
    
  }
  
  return(as.numeric(comb))
}


Bootstrapping<- function(data, no.features, no.simulations,method.combine="Theoretical",method.mean="mean",
                         method.cov="standard",method.cor= "corrected.pearson",neg.class=NULL,lowCI= c(0.025,0.05,0.1)){
 
  combb<- matrix(0,no.simulations,no.features)
  
  lower.ci<- matrix(0, no.features, length(lowCI)+1)
  colnames(lower.ci)<- c("real.Data", paste("confidence.interval",lowCI))
  
  dat<- selectData(data,no.features)
  
  if(method.combine=="Theoretical"){
    combaucs<- combinedAUCs(dat,no.features=no.features,method.combine=method.combine,
                method.mean=method.mean,method.cov=method.cov,method.cor=method.cor)
    
    lower.ci[,1]<- combaucs
    for(i in 1:no.simulations){
      sampdata<-sample(nrow(dat),nrow(dat),replace=T)
      samp<- dat[sampdata,]
   
      comb<- theor_combinedAUCs_sameOrder(samp, no.features,method.mean,method.cov,method.cor) 
      
      combb[i,]<-comb
    }
  }
  
  if(method.combine=="RealData"){
    combaucs<- combinedAUCs(dat,no.features=no.features,method.combine=method.combine,
                            neg.class=neg.class)
    
    lower.ci[,1]<- combaucs
    
    for(i in 1:no.simulations){
      sampdata<-sample(nrow(dat),nrow(dat),replace=T)
      samp<- dat[sampdata,]
      
      comb<- boot.comb.normalized(samp, no.features, neg.class)
      
      combb[i,]<-comb
    }
  }
  
  for(j in 1: length(lowCI)){
    for(i in 1:no.features){
      lower.q<- quantile(combb[,i], lowCI[j])
      
      lower.ci[i,j+1]<- lower.q
    }
  }

  return(data.frame(lower.ci))
}

plot_Bootsrapping_th<- function(data, no.features, no.simulations,
method.mean="mean",method.cov="standard",method.cor="corrected.pearson",
        lowCI= c(0.025,0.05,0.1),xlab="no.of combining features", ylab="AUC value",
        main="Biomarkers Combination",titles=NULL, shape="l",opacity=NULL){

  
    y<- Bootstrapping(data,no.features=no.features, no.simulations=no.simulations,
                      method.combine="Theoretical",method.mean=method.mean,
                      method.cov=method.cov,method.cor=method.cor, lowCI= lowCI)
 
  x<- c(1:no.features)
  
  return(plotLines(x, y, xlab=xlab, ylab= ylab, main=,main, titles=titles,  shape=shape, opacity=opacity))
}

plot_Bootsrapping_rd<- function(data, no.features, no.simulations,negclass="n",
lowCI= c(0.025,0.05,0.1),xlab="no.of combining features", ylab="AUC value",
  main="Biomarkers Combination",titles=NULL, shape="l",opacity=NULL){
 
    y<- Bootstrapping(data,no.features=no.features, no.simulations=no.simulations,
                      method.combine="RealData", neg.class=negclass,lowCI= lowCI)
  
  
  x<- c(1:no.features)
  
  return(plotLines(x, y, xlab=xlab, ylab= ylab, main=,main, titles=titles,  shape=shape, opacity=opacity))
}



##############################################################################################################
###############   compute sample size   ##########################################################
###############################################################################################################
sampleSize<- function(wAUC, n, nplus, no.DataFeatures){
  
  if(n == 0){
    stop("n can't be 0")
  }
  ratio<- nplus/n
  while(T){
    pValue<-compute.p.auc(wAUC, n, nplus)
    correctedPValue<- pValue* no.DataFeatures
    if(is.null(correctedPValue)){
      stop(correctedPValue)
    }
    if(is.null(wAUC) || is.null(n) || is.null(nplus)){
      stop('Input is missing')
    }
    
    if(correctedPValue< 0.05){
      d<-data.frame(n,nplus, correctedPValue)
      colnames(d)<- c("sample size", "postive cases", "corrected P-Value")
      return (d)
    }
    n<- n+1
    
    nplus<- round(n*ratio)
  }
}

calculate.sampleSizeR<-function(wAUC, ratio, no.DataFeatures){

  if(is.null(wAUC)){
    return(NULL)
  }

  if(wAUC<0 | wAUC>1){
    progress$set(message = "Calculating", value = 1)
    stop("The wAUC is not between 0 and 1.")
  }

  if(wAUC< 0.5){
    wAUC<- 1-wAUC
  }

  n<-2
  nplus<- n*ratio
  result <- sampleSize(wAUC, n,nplus,no.DataFeatures)

  return(result)
  
}

calculate.sampleSize<-function(wAUC, n, nplus, no.DataFeatures){
  
  nplus <- as.numeric(nplus)
  n <- as.numeric(n)
  
  if(is.null(wAUC)){
    return(NULL)
  }
  if(wAUC<0 | wAUC>1){
    stop("wAUC is not between 0 and 1 ")
  }
  if(wAUC< 0.5){
    wAUC<- 1-wAUC
  }
  
  
  return(sampleSize(wAUC, n, nplus, no.DataFeatures))
  
}

sampleSize_differentAUCs<-function(auc.values= c(.7,.8,.9,1), n,nplus ,no.DataFeatures){
  
  sample.sizes <- sapply(auc.values, function(x) calculate.sampleSize(x, n,nplus,no.DataFeatures))
  
  sam.size <- sample.sizes[1,]
  return(unlist(sam.size))
  
}

sampleSizeR_differentAUCs<-function(auc.values= c(.7,.8,.9,1), prevalence ,no.DataFeatures){
  
  sample.sizes <- sapply(auc.values, function(x) calculate.sampleSizeR(x, prevalence,no.DataFeatures))
  
  sam.size <- sample.sizes[1,]
  return(unlist(sam.size))
  
}
