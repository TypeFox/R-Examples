#Preamble

# non-standard plot dimensions
verName<-'EAGLE.1' 
width <- "90%"          # narrower
height <- "500px"       # more than 400

#check for new version
#If need2update = TRUE, a message will come up in the main panel
need2update<-'false'
try({
  library(RCurl)
  newestVerName<-strsplit(getURL('https://raw.githubusercontent.com/aaronjfisher/interAdapt/master/version.txt'),split='\n')[[1]][1]
  if(newestVerName!=verName) need2update<-'true'
})

my_plotOutput <- function(...)
	plotOutput(..., width=width, height=height)

# NOTE: header element defaults to h3 instead of h1
my_headerPanel <- function (title, windowTitle = title, h=h3)
{
    tagList(tags$head(tags$title(windowTitle)), div(class = "span12", 
        style = "padding: 10px 0px;", h(title)))
}

pbreak<-HTML('<P CLASS=breakhere>')



# bt and st are called in server.R, and so 
# we don't need to call them again here
# Adjustments for whether we're on 
# RStudio or not (with the onRStudioServer variable)
# are also done in server.R
readHelpTabHTML<- paste0(readLines('help_tab.html'),collapse='') #will be converted to a the welcome (help) page for interAdapt



animationOptions(interval = 5000, loop = FALSE, playButton = NULL, pauseButton = NULL)






  #  _____ _     _               _____           _      
  # /  ___| |   (_)             /  __ \         | |     
  # \ `--.| |__  _ _ __  _   _  | /  \/ ___   __| | ___ 
  #  `--. \ '_ \| | '_ \| | | | | |    / _ \ / _` |/ _ \
  # /\__/ / | | | | | | | |_| | | \__/\ (_) | (_| |  __/
  # \____/|_| |_|_|_| |_|\__, |  \____/\___/ \__,_|\___|
  #                       __/ |                         
  #                      |___/                          



cat('\n Running interAdapt...\n To stop interAdapt, click on the R window and press escape or press the "STOP" button in the app. Then close the interAdapt browser tab.')

shinyUI(pageWithSidebar(


  my_headerPanel("interAdapt - An Interactive Planner for Group Sequential, Adaptive Enrichment Designs"),  



    #      _     _      
    #     (_)   | |     
    #  ___ _  __| | ___ 
    # / __| |/ _` |/ _ \
    # \__ \ | (_| |  __/
    # |___/_|\__,_|\___|
                      
  #Side panel for inputs

  # "Advanced" forces batch mode
  sidebarPanel(

	# Stop button
	actionButton("stopButton", "STOP"),
	br(), br(),

        #TOP PANEL
        selectInput("Which_params", "", c("Show Basic Parameters" = "1",
                "Show Advanced Parameters" = "2", "Show All Parameters and Save/Load Option" = "3") ),

        #SAVE & LOAD

        conditionalPanel(condition="input.Which_params == '3'",
            br(),
            strong('Save current parameters to file:'),br(),
            downloadButton('downloadInputs', 'Save Inputs'),
            br(),br(),
            strong('Load previous parameters from file:'),br(),
            HTML('(Press Apply after loading)'),
            fileInput('uploadCsvInput', '',
                    accept=c('text/csv', 'text/comma-separated-values,text/plain')),
            strong('Upload trial data from file:'),
            fileInput('uploadDataTable', '',
                    accept=c('text/csv', 'text/comma-separated-values,text/plain')),
            strong("Current Parameters:")
          ),

        #INTERACTIVE MODE ONLY IF JUST SLIDERS
          #Note interactive mode is auto-disabled if you're not in the basic parameters view
        conditionalPanel(condition= "input.Which_params== '1'",
          selectInput("Batch", "", c("Batch mode" = "1",
                "Interactive mode" = "2"))
        ),
        #BASIC SLIDERS
        conditionalPanel(condition = "input.Which_params == '1' || input.Which_params == '3'",
          #show apply button if you're in batch mode, or if we're showing all inputs and batch mode is enforced
          conditionalPanel(condition='input.Batch== "1" || input.Which_params== "3" ',
            actionButton("Parameters1", "Apply")),
          br(), br(),
          uiOutput('fullSliders')
        ),
        #ADVANCED BOXES
        conditionalPanel(condition = "input.Which_params == '2' || input.Which_params == '3'",
          #always show apply button
          actionButton("Parameters2", "Apply"),
          br(), br(),
          uiOutput('fullBoxes'))
  ),





  #                  _       
  #                 (_)      
  #  _ __ ___   __ _ _ _ __  
  # | '_ ` _ \ / _` | | '_ \ 
  # | | | | | | (_| | | | | |
  # |_| |_| |_|\__,_|_|_| |_|
                           
                           

  mainPanel(
  #INITIALIZE PAGE BREAK CODE
  HTML('<STYLE TYPE=text/css> P.breakhere {page-break-before: always} </STYLE>'),

  #STOP alert
  h4(textOutput('stop')),

  #UPDATE MESSAGE
  conditionalPanel(need2update, 
    h4('Updates are available!'),
    'An updated version of this software can been downloaded ',
    a('here',href='https://rawgithub.com/aaronjfisher/interAdapt/master/install_instructions.html'), 
    br(),br()
  ),

  #WARNINGS
  h4(textOutput('warn1')),
  h4(textOutput('warn2')),
  h4(textOutput('warn3')),
  #br(),br(),

  #OUTPUT
  radioButtons("OutputSelection", em(strong("Output selection")),
  c("About interAdapt" = "1", "Designs" = "2", "Performance" = "3"), selected="1"),
  
  br(), pbreak,

  #About interAdapt Tab
  conditionalPanel(condition = "input.OutputSelection == '1'",
    HTML(paste(readHelpTabHTML,collapse=''))),

  conditionalPanel(condition = "input.OutputSelection == '2'",
    em(strong("Designs")),
    tabsetPanel(
    	tabPanel("Adaptive",
    		my_plotOutput("adapt_boundary_plot"),
        br(),pbreak,
        tableOutput("adaptive_design_sample_sizes_and_boundaries_table"),
        downloadButton('downloadDesignAD.1', 'Download table as csv'),br()),
    	tabPanel("Standard, Total Population",
        my_plotOutput("standard_H0C_boundary_plot"),
        br(),pbreak,
    		tableOutput("standard_H0C_design_sample_sizes_and_boundaries_table"),
        downloadButton('downloadDesignSC.1', 'Download table as csv'),br()) ,
    	tabPanel("Standard, Subpop. 1 only",
    		my_plotOutput("standard_H01_boundary_plot"),
        br(),pbreak,
        tableOutput("standard_H01_design_sample_sizes_and_boundaries_table"),
        downloadButton('downloadDesignSS.1', 'Download table as csv'),br()),
    	tabPanel("All designs",
    		tableOutput("adaptive_design_sample_sizes_and_boundaries_table.2"),
        pbreak,
    		tableOutput("standard_H0C_design_sample_sizes_and_boundaries_table.2"),
        pbreak,
        tableOutput("standard_H01_design_sample_sizes_and_boundaries_table.2"),
        br(),br(),downloadButton('downloadDesignAD.2', 'Download AD design table as csv'),
        br(),br(),downloadButton('downloadDesignSS.2', 'Download SS design table as csv'),
        br(),br(),downloadButton('downloadDesignSC.2', 'Download SC design table as csv')),
  selected="Adaptive")
  ),

  conditionalPanel( condition = "input.OutputSelection == '3'",
    tabsetPanel(
    	tabPanel("Power", my_plotOutput("power_curve_plot")),
    	tabPanel("Sample Size", my_plotOutput("expected_sample_size_plot")),
    	tabPanel("Duration", my_plotOutput("expected_duration_plot")),
      selected='Power'),
    pbreak,
    tableOutput("performance_table"),
    downloadButton('downloadPerformance.1', 'Download table as csv')
  ),
  conditionalPanel(condition = "input.OutputSelection != '1'",
    br(),
    downloadButton('knitr', 'Generate report')
  ),
  br(),br() )
))
