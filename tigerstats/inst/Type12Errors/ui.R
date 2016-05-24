library(shiny)




# Define UI for testing application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("Types of Error in Significance Tests"),
  
  # Sidebar
  sidebarPanel(
    helpText("You will work with a population that is normally distributed, with a",
             "mean of 70 and a standard deviation of 3.  You can set mu0, the value",
              "that the Null Hypothesis believes mu to be."),
    sliderInput(inputId="mu0", "mu0, the Null's belief about mu",
                                    value=70,min=65,max=75),
    br(),
    
    helpText("Choose the level of significance:  the cut-off value alpha.  If the two-sided",
             "P-value is below alpha, the Null Hypotheses will be rejected."),
    sliderInput(inputId="alpha","Significance level alpha",value=0.05,min=0.01,max=0.20,step=0.01),
    br(),
    
    helpText("One sample at a time, or 5000 at a time?  (Your choice affects the sample",
             "sizes you may select from."),
    selectInput(inputId="actionType",label="One Sample or 5000:",
                choices=list("Just one sample at a time, please."="one",
                             "Give me 5000 samples!"="fiveThousand")),
    
    uiOutput("sizeSelect"),
    br(),
    
    
    
    helpText("When you are ready to sample, push this button:"),
    actionButton("go","Take the Sample")
  ),
  
  
  # Here comes the main panel
  
  mainPanel(
    
    conditionalPanel(
      condition="output.go == 0",
      plotOutput("initialGraph"),
      HTML("<p> </p>"),
      HTML("<ul>
                <li>The population density curve is in red.</li>
                <li>The red vertical line marks the population mean mu.</li>
                <li>The black vertical line marks mu0, the value the Null believes in.</li>
                <li>(When the red and black lines coincide, the Null is true.)</li>
          </ul>")
    ),
    
    conditionalPanel(
      condition="(output.go > 0) && (output.actionType == 'one')",
      plotOutput("graphSample"),
      HTML("<p> </p>"),
      HTML("<ul>
                <li>The population density curve is in red.</li>
                <li>The red line marks the population mean mu.</li>
                <li>The black vertical line marks mu0, the value the Null believes in.</li>
                <li>(When the red and black lines coincide, the Null is true.)</li>
                <li>The histogram of the sample is in light blue.</li>
                <li>The sample mean is the big blue dot.</li>
                <li>The confidence interval is in green.</li>
          </ul>"),
      tableOutput("results")
    ),
    
    conditionalPanel(
      condition="(output.go > 0) && (output.actionType == 'fiveThousand')",
      plotOutput("initialGraph2"),
      tableOutput("summary"),
      dataTableOutput("intervalFrame")      
    )

    
  ) # end main panel
  
))
