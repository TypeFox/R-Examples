#library(shiny)
#library(diveRsity)
if(!require("shinyIncubator")){
  stop("The package 'shinyIncubator' is required to use this application.")
}
# Define UI
shinyUI(
  
  pageWithSidebar(
    
    headerPanel("divMigrate-online: Visualise and test gene flow patterns among populations"),
    
    sidebarPanel(      
      fileInput("file", h5("1. Input file"), multiple = FALSE, accept = NULL),      
      numericInput("nbs", h5("2. Number of bootstraps"), value = 0, min = 0, 
                   max = 5000, step = 1),      
      radioButtons("stat", h5("3. Migration Statistic"),
                   c("D", "Gst", "Nm"), selected = "D"),      
      sliderInput("filter_threshold", h5("4. Filter Threshold"), min = 0, 
                  max = 1, value = 0, step = 0.05),
      conditionalPanel(
        condition = "input.tabs=='Network Plots'",
        uiOutput("popnames")
      ),
      conditionalPanel(
        condition = "input.tabs=='Network Plots'",
        uiOutput("stdplt")
      ),
      conditionalPanel(
        condition = "input.tabs=='Network Plots'",
        uiOutput("pltFormat")
      ),
      HTML("<br>"),
      conditionalPanel(
        condition = "input.tabs=='Network Plots'",
        uiOutput("pltDL")
      ),
      HTML("<br>"),
      #actionButton("goButton", h5("Calculate")),
      helpText(""),      
      helpText("Written and designed by Kevin Keenan, using shiny",
               "from RStudio and Inc. (2012).")
    ),
    
    mainPanel(progressInit(),
              tabsetPanel(id = "tabs",
                          tabPanel(
                            HTML("<h5><font color = #215F9C>Usage Instructions</font></h5>"),
                            HTML("<h3><font color = #215F9C>Welcome to the divMigrate-Online       web app!</font></h3><br>"),
                            h4("Getting started:"),
                            HTML("<ol>",
                                 "<li>Read your genepop file (2-digit or 3-digit)</li>",
                                 "<li>Choose which statistic you want to use to calculate",
                                 "relative migration</li>",
                                 "<li>Statistical significance of directional migration ",
                                 "can be calculated by setting the number of bootstraps ",
                                 " to a value > 0</li>",
                                 "<li>The number of links in the network plot can be ",
                                 "manipulated using the 'filter threshold' slider and by excluding populations.</li>",
                                 "<li>Download your result matrix (as .txt) for reporting ",
                                 "puposes (make sure you have cliked calculate before trying to download your file).</li><br>"),
                            HTML("<font size = 1>Any issues or queries regarding this application can be directed to kkeenan02[at]qub.ac.uk.<br> The methods used in this application were originally presented in the publications below. Please cite them.<br><b>Sundqvist et al, (in prep.)<br>Keenan et al, (in prep.)</b></font>")
                          ),
                          
                          tabPanel(
                            HTML("<h5><font color = #215F9C>Network Plots</font></h5>"),
                            plotOutput("plt", "auto", "650px")
                          ),
                          
                          tabPanel(
                            HTML("<h5><font color = #215F9C>Results Matrix</font></h5>"),
                            h6("To download a file containing your results as a tab",
                               "delimited file, please click the button below:"),
                            downloadButton("mat", "Download your file!")
                          )
              )   
    )
  )
)