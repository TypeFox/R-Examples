library(shiny)

shinyUI(pageWithSidebar(
  
  headerPanel("Amplification curve simulation"),
  
  sidebarPanel(
    HTML('<p><img src="https://raw.githubusercontent.com/michbur/chipPCR/master/vignettes/logo.png" width="55%" height="55%"/></p>'),
    numericInput("cycles", "Cycles", 35, min = 10, max = 60),
    numericInput("b.eff", "Efficiency", -25), 
    numericInput("bl", "Baseline", 0.05), 
    numericInput("ampl", "Amplitude", 1), 
    numericInput("Cq", "Cq value", 20),
    checkboxInput("noise", "Use noise in simulation", FALSE),
    numericInput("nnl", "Noise level", 0.025, min = 0, max = 10, step = 0.025),
    selectInput("nnl.method", "Variable:",
                list("Constant" = "constant", 
                     "Decreasing" = "decrease", 
                     "Increasing" = "increase")),
    numericInput("th.r", "Fluorescence threshold value", 0.7, step = 0.05),
    checkboxInput("th.auto", "Automatic estimation of the threshold (experimental)", FALSE),
    checkboxInput("th.lin", "Linear regression instead of quadratic", TRUE)
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Amplification plots", plotOutput("AmpSimPlot"), 
               verbatimTextOutput("bgSummary"),
               plotOutput("inderPlot")),
      tabPanel("Amplification data", tableOutput("bgTable"))
    )
  )
)
)
