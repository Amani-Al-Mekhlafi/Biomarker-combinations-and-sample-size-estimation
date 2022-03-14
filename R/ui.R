
ui <- fluidPage(
  fileInput(inputId = "data", label = "Choose a file.", accept = c("text/csv", "text/comma-separated-values, text/plain", ".csv")), 
  tabsetPanel(
    tabPanel("HAUCA",
             tabsetPanel(
               tabPanel("Table", dataTableOutput("table.output")),
               tabPanel("Curve",
                        plotlyOutput("hauca.curve.output")
                        
               )
             )
    ),
    
    tabPanel("Combination",
tabsetPanel(
  
tabPanel("Theoretical",

     selectInput(inputId= "correlation", label= "Type of Correlation", choices= c("No Correlation","kendall","pearson", "corrected.pearson")),
     selectInput(inputId= "Method.of.mean", label= "Method of Average", choices= c("mean", "median")),
     selectInput(inputId= "Method.of.sd", label= "Measure of Dispersion", choices= c("standard", "covRob")),
     splitLayout(
     numericInput(inputId = "no.features", label = "Number of Features", value = 5, min = 2, max = 100, step =1),
     numericInput(inputId = "low.CI", label = "Lower Confidence Interval", value = 0.05, min = 0, max = 1, step = 0.001), 
     numericInput(inputId = "no.simulations", label = "Number of Simulations.", value = 1000, min = 100, max = 1000000, step = 100)
     ),
       
     actionButton("theoretical.output.button", "Visualize"),
     br(),
     plotlyOutput("theoretical.output"),
     br()
     
     
),

tabPanel("Real Data",
numericInput(inputId = "nofeatures", label = "Number of Features", value = 5, min = 2, max = 100, step =1),
textInput(inputId = "neg.class", label = "Enter the negative class."),
numericInput(inputId = "lowCI", label = "Lower Confidence Interval", value = 0.05, min = 0, max = 1, step = 0.001), 
numericInput(inputId = "nosimulations", label = "Number of Simulations.", value = 1000, min = 100, max = 1000000, step = 100),


actionButton("normalize.output.button", "Visualize"),
br(),
plotlyOutput("normalize.output"),
br()
)



)
)   ,
    
    

tabPanel("Sample Size",
numericInput(inputId = "no.of.multiple.testing", label = "Choose the number of multiple testing", value = 2000, min = 1, max = 100000, step = 1),
tabsetPanel(
tabPanel("By Prevalence",
numericInput(inputId = "Prevalence", label = "Choose a Prevalence.", value = 0.5, min = 0, max = 1, step = 0.01),
tabsetPanel(
tabPanel("Calculate",
  numericInput(inputId = "wAUC.value", label = "Choose the AUC value", value = 0.8, min = 0.5, max = 1, step = 0.01),
         
actionButton("samplesize.r.button", "Calculate"),
tableOutput("samplesize.numbers.ratio")
),
tabPanel("Visualize",


splitLayout(
numericInput(inputId = "no1", label = "From auc.val",value = 0.5, min = 0.5, max = 1, step = 0.01),
numericInput(inputId = "no2", label = "To auc.val",value = 0.5, min = 0.5, max = 1, step = 0.01),

br(), br(), br()
),



actionButton("samplesize.r.graph.button", "Visualize"),
br(),
plotlyOutput("samplesize.graph.ratio"),
br()

)
)
),
tabPanel("By Numbers",
numericInput(inputId = "n", label = "Number of Observations", value = 100, min = 1, max = 1000, step = 1),
numericInput(inputId = "nplus", label = "Number Positive Observations", value = 50, min = 1, max = 1000, step = 1),
tabsetPanel(
tabPanel("Calculation",
  numericInput(inputId = "wAUC.value", label = "Choose the AUC value", value = 0.8, min = 0.5, max = 1, step = 0.01),
         
actionButton("samplesize.nnplus.button", "Calculate"),
tableOutput("samplesize.numbers.nnplus")
),
tabPanel("Visualization",


splitLayout(
numericInput(inputId = "no1", label = "From auc.val",value = 0.5, min = 0.5, max = 1, step = 0.01),
numericInput(inputId = "no2", label = "To auc.val",value = 0.5, min = 0.5, max = 1, step = 0.01),
br(), br(), br()
),

actionButton("samplesize.nnplus.graph.button", "Visualize"),
br(),
plotlyOutput("samplesize.graph.nnplus"),
br()
    )
                        )         
               )
               
             )
    )
  )
)