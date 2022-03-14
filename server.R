source("AmanisPlots2.R")

## Libraries
if("plotly" %in% rownames(installed.packages())){
  library(plotly)} else{
    install.packages("plotly",dependencies = T)
    library(plotly)}

if("DT" %in% rownames(installed.packages())){
  library(DT)} else{
    install.packages("DT", dependencies = T)
    library(DT)}

if("ROCR" %in% rownames(installed.packages())){
  library(ROCR)} else{
    install.packages("ROCR", dependencies = T)
    library(ROCR)}

if("robust" %in% rownames(installed.packages())){
  library(robust)} else{
    install.packages("robust", dependencies = T)
    library(robust)}

if("corpcor" %in% rownames(installed.packages())){
  library(corpcor)} else{
    install.packages("corpcor", dependencies = T)
    library(corpcor)}

if("foreach" %in% rownames(installed.packages())){
  library(foreach)} else{
    install.packages("foreach",dependencies = T)
    library(foreach)}

if("doParallel" %in% rownames(installed.packages())){
  library(doParallel)} else{
    install.packages("doParallel",dependencies = T)
    library(doParallel)}

server <- function(input, output){
## HAUCA table  
  output$table.output <- DT::renderDataTable({
    data.input <- input$data
    if(!is.null(data.input)){
      HAUCA.table(read.csv2(data.input$datapath))
    }
  })
  
## HAUCA curve  
  output$hauca.curve.output <- renderPlotly({
    data.input <- input$data 
    if(!is.null(data.input)){
      
      HAUCA.curve(read.csv2(data.input$datapath))
      
    } 
  })

## Biomarker Combinations (Theoretically) 
  curve<- eventReactive(input$theoretical.output.button, {
    data.input <- input$data
    if(!is.null(data.input)){
      isolate(plot_Bootsrapping_th(read.csv2(data.input$datapath), 
               input$no.features, 
               input$no.simulations,input$Method.of.mean, 
               input$Method.of.sd,input$correlation, input$low.CI, 
               xlab="no.of combining features", ylab="AUC value", 
               main="Biomarker Combinations with Confidence Interval by Bootstrapping 
               (Theoretically)",
               titles=NULL, shape="l",opacity=NULL))

    }
  })
  
  output$theoretical.output <- renderPlotly({curve()})
    
## Biomarker Combinations (Based on Real Data)    
  
  output$normalize.output <- renderPlotly({
    input$normalize.output.button
    if(input$normalize.output.button != 0){
      data.input <- input$data

    plot_Bootsrapping_rd(read.csv2(data.input$datapath), input$nofeatures,
          input$nosimulations,input$neg.class,input$lowCI, 
          xlab="no.of combining features", ylab="AUC value", 
          main="Biomarker Combinations with Confidence Interval by Bootstrapping
           (Based on Real Data)",
          titles=NULL, shape="l",opacity=NULL)     
    }
    
  })
  
## calculate the sample size  
  
  output$samplesize.numbers.ratio <- renderTable({
    input$samplesize.r.button
    if(input$samplesize.r.button != 0){
      isolate(calculate.sampleSizeR(input$wAUC.value, input$Prevalence, input$no.of.multiple.testing))
    }
  })
  
  output$samplesize.graph.ratio <- renderPlotly({
    
    input$samplesize.r.graph.button
    if(input$samplesize.r.graph.button != 0){
      isolate(auc.values <- seq(input$no1, input$no2, 0.01))
      isolate(sam_sizeR <- sampleSizeR_differentAUCs(auc.values, input$Prevalence ,input$no.of.multiple.testing))
      isolate(plotLines(auc.values, data.frame(sam_sizeR), xlab="AUC value", ylab="Sample Size", main= "sample size estimation",  shape="l", opacity=NULL))
      
    }
  })
  
  output$samplesize.numbers.nnplus <- renderTable({
    
      input$samplesize.nnplus.button
    if(input$samplesize.nnplus.button != 0){
      isolate(calculate.sampleSize(input$wAUC.value, input$n, input$nplus, input$no.of.multiple.testing))
    }
  })
  
   output$samplesize.graph.nnplus <- renderPlotly({
    input$samplesize.nnplus.graph.button
    if(input$samplesize.nnplus.graph.button != 0){
      isolate(auc.values <- seq(input$no1, input$no2, 0.01))
      isolate(sam_size <- sampleSize_differentAUCs(auc.values, input$n, input$nplus ,input$no.of.multiple.testing))
      isolate(plotLines(auc.values, data.frame(sam_size), xlab="AUC value", ylab="Sample Size", main= "sample size estimation",  shape="l", opacity=NULL))
 }
  } )
   
}
#shinyApp(ui, server)
