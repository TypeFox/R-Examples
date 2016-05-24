library(shiny)




# Define UI for testing application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("The Central Limit Theorem"),
  
  # Sidebar
  sidebarPanel(
    conditionalPanel(
      condition="input.resample == 0 || output.totalPrev == output.total",
      selectInput(inputId="popDist",label="Population Shape",
                choices=list("Normal"="normal",
                             "Skewy"="skew",
                             "REALLY Skewed"="superskew",
                             "Way-Out Outlier Group"="outliers")),
      br(),
    
      helpText("Choose the sample size."),
    
      sliderInput(inputId="n","Sample Size n",value=2,min=2,max=50,step=1),
      br()
    ),
      helpText("How many samples would you like to take at one time?  Limit is 10000. We will ",
             "compute the mean of each sample."),
      numericInput("sims","Number of Samples at Once",1,min=0,step=1),
      actionButton("resample","Sample Now"),
      conditionalPanel(
        condition="(input.resample > 0 && input.reset == 0) || output.total > output.totalPrev",
        actionButton("reset","Start Over")
      )
    
  ),
  
  
  # Here comes the main panel
  
  mainPanel(
    
    conditionalPanel(
      condition="input.resample == 0 || output.totalPrev == output.total",
      plotOutput("initialGraph")
    ),

    
    conditionalPanel(
      condition="(input.resample > 0 && input.reset == 0) || output.total > output.totalPrev",
      
      tabsetPanel(
        tabPanel("Latest Sample",
                 plotOutput("graphSample"),
                 HTML("<p> </p>"),
                 HTML("<ul>
                        <li>The population density curve is in red.</li>
                        <li>The vertical line marks the population mean.</li>
                        <li>The histogram of the sample is in light blue.</li>
                        <li>The sample mean is the big blue dot.</li>
                      </ul>")),
        tabPanel("Sample Mean Distribution",
                 plotOutput("xbar"),
                 HTML(
                   "<p>The plots above compare the distribution of the sample means so far to a normal curve.</p>
                     <p></p>
                     <ul>
                        <li>The normal curve is in red.  If the population is exactly normal, then this curve represents the exact distribution of the sample mean.</li>
                        <li>A density plot of the sample means is in blue.  This plot gives a pretty good estimate of the distribution of the sample mean, for the population and sample size that you have selected.</li>
                        <li>The rug at the bottom of the plot gives individual sample means found so far (if there are no more than 50).  The latest one is marked with a large blue dot.</li>
                    </ul>"))
      ) # end tabset panel
    ) # end conditonal panel
    
  ) # end main panel
  
))
