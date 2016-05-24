## Use navbarPage() to add a related app

library(shiny)

simLimit <- 10000

####################################3
## ui
####################################

navbarPage(
  title = 'Exploring Confidence Intervals',
  tabPanel(title = "Coverage Properties",
           sidebarPanel(
           conditionalPanel(
             condition = "input.takeSample == 0 || output.beginning == true",
             selectInput(
               inputId="popDist",label="Population Shape",
               choices=list("Normal"="normal",
                            "Skewy"="skew",
                            "REALLY Skewed"="superskew",
                            "Way-Out Outlier Group"="outliers")
               ),
             br(),
             sliderInput(inputId="n","Sample Size n",value=2,min=2,max=50,step=1),
             helpText("How confident do you want to be that the population mean is contained",
                      "within the confidence interval?   Use the slider to select a desired",
                      "percent-confidence level."),
             sliderInput(inputId="confLevel","Confidence Level",value=80,min=50,max=99,step=1)
             ),
             helpText("How many samples would you like to take at one time?  Limit is 10000. With each ",
                      "sample, we'll make a confidence interval for the population mean."),
             numericInput("sims", "Number of Samples at Once", 1, min=0, max = simLimit, step=1),
             actionButton("takeSample","Sample Now"),
             conditionalPanel(
               condition = 'output.beginning == false',
               actionButton("reset","Start Over")
               )
           ),
           mainPanel(
             conditionalPanel(
               condition = "input.takeSample == 0 || output.beginning == true",
               plotOutput("initialGraph"),
               HTML("<ul>
                       <li>The population density curve is in red.</li>
                       <li>The vertical line marks the population mean.</li>
                    </ul>")
               ),
             conditionalPanel(
               condition = "output.beginning == false",
               tabsetPanel(
                 id = "coverageTabsetPanel",
                 tabPanel(
                   title = "Latest Interval",
                   plotOutput("plotSample"),
                   HTML("<p> </p>"),
                   HTML("<ul>
                          <li>The population density curve is in red.</li>
                          <li>The vertical line marks the population mean.</li>
                          <li>A density curve for the most recent sample is in light blue.</li>
                          <li>The sample mean is the big blue dot.</li>
                          <li>The confidence interval is in green.</li>
                        </ul>"),
                   br(''),
                   tableOutput("summary")
                   ),
                 tabPanel(
                   title = "t-statistic",
                   plotOutput("tstatistic"),
                   HTML("<p>The plots above compare the actual distribution of the t-statistic to the t-curve with n-1 degrees of freedom.</p>
                        <p></p>
                        <ul>
                          <li>The t-curve is in red.  If the population is exactly normal, then this curve represents the exact distribution of the t-statistic.</li>
                          <li>The density plot of the t-statistics found so far is shown in blue.  This plot gives a pretty good estimate of the actual distribution of the t-statistic, for the population and sample size that you have selected.</li>
                        </ul>")
                   )
                 )
               )
             )
           ),
  tabPanel(
    title = "Fifty at a Time",
    sidebarPanel(
      conditionalPanel(
        condition = "input.takeSample2 == 0 || output.beginning2 == true",
        selectInput(
          inputId="popDist2",
          label="Population Shape",
          choices=list("Normal"="normal",
                       "Skewy"="skew",
                       "REALLY Skewed"="superskew",
                       "Way-Out Outlier Group"="outliers")
          ),
        br()
        ),
      sliderInput(inputId="n2","Sample Size n",value=2,min=2,max=50,step=1),
      helpText("How confident do you want to be that the population mean is contained",
               "within the confidence interval?   Use the slider to select a desired",
               "percent-confidence level."),
      sliderInput(inputId="confLevel2","Confidence Level",value=80,min=50,max=99,step=1),
      actionButton("takeSample2","Fifty Samples Now"),
      conditionalPanel(
        condition = 'output.beginning2 == false',
        actionButton("reset2","Start Over")
        )
      ),
    mainPanel(
      conditionalPanel(
        condition = "input.takeSample2 == 0 || output.beginning2 == true",
        plotOutput("initialGraph2"),
        HTML("<ul>
                <li>The population density curve is in red.</li>
                <li>The vertical line marks the population mean.</li>
             </ul>")
        ),
      conditionalPanel(
        condition = "output.beginning2 == false",
          plotOutput("plotSample2"),
          HTML("<ul>
                  <li>The population density curve is in red.</li>
                  <li>The vertical line marks the population mean.</li>
                  <li>Intervals covering the mean are green.</li>
                  <li>Intervals NOT covering the mean are in burlywood.</li>
               </ul>"),
          tableOutput("summary2")
        )
      )
    )
  )