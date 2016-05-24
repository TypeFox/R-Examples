################################################################################
# Define UI ####################################################################
shinyUI(pageWithSidebar(


  # Application title
  headerPanel(windowTitle="Density based peak picking",title="enviPick 1.0"),

  sidebarPanel(
    fileInput("file1", "Choose centroided .mzXML input file"),
    radioButtons("MSlevel", "MS level", c("1"="1","2"="2"), selected = NULL),
    actionButton("Calculate","Calculate"),
    tags$h6("Progress:"),
    verbatimTextOutput("text2"),
    verbatimTextOutput("error"), 
    tags$h6("Job status:"),
    verbatimTextOutput("text3"),
    verbatimTextOutput("text4"),
    verbatimTextOutput("text1"),
    verbatimTextOutput("text5")       
  ),

  mainPanel(
  
  tags$h4("Parameter settings"),
     tabsetPanel(
      tabPanel("(1) EIC partitioning & clustering",
          sliderInput("drtgap", "Maximum retention time gap in an EIC", min = 20, max = 1500, value = 300, step= 1),
          #sliderInput("drtdens1", "Retention time (RT) tolerance for clustering; defined as (+/-) seconds relative to the lowest and highest RT value in each cluster ", min = 10, max = 1500, value = 60, step= 1),
          sliderInput("dmzdens", "Maximum deviation (+/-) of m/z of a measurement from its EIC mean [ppm]", min = 1, max = 20, value = 3.5, step= 0.1)       
       ),
       tabPanel("(2) Peak definition",
          numericInput("minpeak", "Minimum number of measurements per peak ...", 4),
          sliderInput("drtsmall2", "... within a given RT window [s]", min = 5, max = 200, value = 20, step= 0.1),
          sliderInput("drtfill", "Maximum RT gap length to be interpolated [s]", min = 0, max = 60, value = 10, step= 0.1),
          sliderInput("drtdens2", "Peak definition - Maximum RT length of a single peak", min = 10, max = 1500, value = 120, step= 0.1),
          sliderInput("minint", "Minimum log10(intensity) threshold", min = 0, max = 10, value = 4, step= .1),
          numericInput("SN", "Minimum Signal/Noise", 5),
          numericInput("SB", "Minimum Signal/Base", 2),
          numericInput("recurs", "Maximum possible number of peaks within a single EIC", 3)
       ),
       tabPanel("(3) Advanced",
          numericInput("ended", "How often can a peak detection fail to end the recursion? - peak picking", 1),
          numericInput("weight", "Weight for assigning measurements to a peak - peak picking", 1),
          sliderInput("maxint", "Upper log10(intensity) safety threshold", min = 0, max = 15, value = 6.7, step= .1)
       )
     ),
     
     conditionalPanel(
      condition = "output.text4 == 'calculation completed'",
      tags$h4("Results"),
      tabsetPanel(
      tabPanel("Summary",
        tags$h6("Total number of measurements:"),      
          verbatimTextOutput("sum_1"),
        tags$h6("Number of partitions:"),     
          verbatimTextOutput("sum_2"),      
        tags$h6("Number of measurements in partitions:"), 
          verbatimTextOutput("sum_3"),      
        tags$h6("Number of EICs:"), 
          verbatimTextOutput("sum_4"),
        tags$h6("Number of measurements in EICs:"), 
          verbatimTextOutput("sum_5"),        
        tags$h6("Number of peaks:"), 
          verbatimTextOutput("sum_6"),
        tags$h6("Number of measurements in peaks:"), 
          verbatimTextOutput("sum_7")                  
      ),
      tabPanel("Plot",
        plotOutput("plotit1"),
        plotOutput("plotit2"),
        plotOutput("plotit3"),
        plotOutput("plotit4"),
        actionButton("iRplot","Start interactive R plot")                
      ),
      tabPanel("Peak table & download",
        downloadButton('downloadData', 'Download'),      
        tableOutput("peaks")
      )
     )
      
   )
  )

))
################################################################################
