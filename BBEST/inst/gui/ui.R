############################################################################
##
##
##  TITLE PANEL
##  SIDEBAR PANEL
##    MAIN MENU
##    LOCAL MENU
##      LOAD DATA 
##      PREPARE DATA 
##        TRUNCATE DATA
##        LAMBDA
##        BASELINE
##        NOISE
##      G(r) INFORMATION
##      DIFFERENTIAL EVOLUTION
##      FIT RESULTS
##    Plot Bounds
##  MAIN PANEL

# CB Check Box
# RB Radio Button
# TI Text Input
#

# mainRB: RadioButton Main Menu
# data file: upload data
# headerCB: datafile with header
#

shinyUI(fluidPage(
############################################################################
##   === TITLE PANEL ===

 # titlePanel(title="Bayesian Background Estimation", windowTitle="BBEST"),
  tags$head(
    tags$title('BBEST'),
    h2(span("B", style="letter-spacing: -0.3em; color:#0055FF"), span("ayesian ", style = "color:#0099FF;"),
       span("B", style="letter-spacing: -0.3em; color:#0055FF"), span("ackground", style = "color:#0099FF;"),
       span("E", style="letter-spacing: -0.3em; color:#0055FF"), span("S", style="letter-spacing: -0.3em; color:#0055FF"),
       span("T", style="letter-spacing: -0.3em; color:#0055FF"), span("imation", style = "color:#0099FF;") )
  ),


############################################################################
######################################
##   === SIDEBAR PANEL ===
  sidebarLayout(
    sidebarPanel(
######################################
##   == MAIN MENU ==
      tags$table(style="border: 3px solid #E8E8E8; border-style:solid; width: 100%;" ,
	      tags$tr(
		      tags$td(
            h3("Main Menu"),
            p(),
            radioButtons('mainRB', '',
              choices=c("Load data"='load',
                  "Set additional parameters"='prepare',
                  "Set real-space condition"='gr',
                  "Optimize background with DifEv"='difev',
                  "Fit results"='save',
                  "Plot options"='plot'),
              selected='load'
            )
          )
	      )
	    ),

######################################
##   == LOCAL MENU ==
      tags$hr(),
      tags$table(style="border: 3px solid #E8E8E8; padding:15px;border-style:solid; width: 100%;" ,
        tags$tr(
          tags$td(
##   = LOAD DATA =
            conditionalPanel(
              condition = "input.mainRB == 'load'",
              h3("Load Data"),
              p(),
              fileInput('datafile', strong('Choose RData, CSV, text, .sqa or .sqb file'),
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.sqa', '.sqb', '.sq', '.RData')),
              p(),
              checkboxInput('headerCB', strong('Data include header'), TRUE),
              p(),
              radioButtons('separatorRB', strong('Separator'),
                c(Comma=',',
                  Semicolon=';',
                  Tab='\t',
                  Space=' '),
                selected='\t'
              ),
              p(),
              uiOutput("sqaSplit")
            ),
##   = PREPARE DATA =
            conditionalPanel(
              condition = "input.mainRB == 'prepare'",
              h3("Prepare data"),
              p(),
##     TRUNCATE DATA
              strong("Truncate data"),
              uiOutput("truncLimitsR"),
              p(),
##     LAMBDA
              strong("Useful signal level"),
              textInput("lambda", label = c("Type x_1, lambda_1,   x_2, lambda_2, lambda_0"), value = ""),
              p(),
##     BASELINE
              strong("Baseline"),
              checkboxInput("setSB", "Set/Recalculate baseline"),
              conditionalPanel(
                condition = "input.setSB == true",
                textInput("SBNAtoms", label = c("Type number of atoms of each type per unit cell"), value = ""),
                textInput("SBScLen", label = c("Type neutron scattering lengths"), value = ""),
                checkboxInput("oneADP", "use one ADP for all atoms"),
                checkboxInput("fitADP", "fit ADP(s)"),
                conditionalPanel(
                  condition = "input.fitADP == false",
                  textInput("ADP", label = c("Type ADP(s)"), value = "")
                )
              ),
     #          p(),
               tags$hr(),
##     NOISE
              strong("Noise level"),
              textInput("sigma", label = c("Type number of regions or bounds for a signal-free region"),  value = ""),
              textInput("sigmaTS", label = c("Type threshold scale (degree of smoothing)"),  value = "1"),
              p(style="margin:0; padding:0;"),
              actionButton("calcSigmaButton", label = "Estimate noise"),
              p(),              
##     P.BKG    
              strong("P(bkg)"),
              numericInput("pbkg", min=-1, max=1, step=0.01,
                label = "Type probability that data points contain no signal contribution (only background). Type '-1' to estimate P(bkg) iteratively",
                value = 0.5),
              p(),
##     G(R) PLOT    
              strong("G(r)"),
              checkboxInput("estGr", "Estimate PDF"),
              conditionalPanel(
                condition = "input.estGr == true",
                textInput("rGrid", label = c("Type minimum r, maximum r, grid spacing dr"), value = ""),
                p(style="margin:0; padding:0;"),
                actionButton("plotPrelimGr", label = "Plot G(r)")                
              )              
              
              
            ),


##   = G(r) INFORMATION =
            conditionalPanel(
              condition = "input.mainRB == 'gr'",
              p(),
              uiOutput("GrNoteForBanks"),
              h3("Low-r Correction"),
              p(),
              checkboxInput("inclGr",
                strong(HTML(paste("Use low-", em("r"), " conditions in ", em("G(r) "),
                  "to construct the ", em("r"), "-space likelihood", sep = "")))),
              tags$hr(),
              conditionalPanel(
                condition = "input.inclGr == true",
                selectInput("GrNoiseType", label = strong("Condition type"),
                  choices = list("Gaussian noise" = "gauss",
                      "Correlated noise" = "correlated"),
                  selected = "gauss"),                          
                numericInput("rhoInclGr", value=NA, min=0, max=1,  step=0.0000001,
                  label = strong(HTML(paste("Number density of the material &rho;", tags$sub(0), sep = "")))),
                numericInput("rminInclGr", min=0, max=2.5, step=0.001,
                  label = strong("minimum r"),
                  value = 0),
                numericInput("rmaxInclGr", min=0, max=2.5, step=0.001,
                  label = HTML(paste(strong("maximum r"), "(should be below the shortest possible interatomic distance)")),
                  value = 2),
                numericInput("drInclGr", min=0.001, max=0.05, step=0.001,
                  label = div( span(strong("grid spacing")), span(strong(em("dr"))) ),
                  value = 0.005),
                p(),  
                actionButton("setGrButton", label = strong("Submit")),
                p(),
                HTML(paste("(make sure", strong(HTML("&epsilon;")), "(noise level) has been estimated)"))
              )
            ),
         #   span("B", style="letter-spacing: -0.3em;")


##   = DIFFERENTIAL EVOLUTION =
            conditionalPanel(
              condition = "input.mainRB == 'difev'",
              h3("DifEv Parameters"),
              p(),
              numericInput("fitNP", min=10, max=1000, step=10,
                label = strong("Number of population members"),
                value = 100),
              p(),
              numericInput("fitItermax", min=10, max=10000, step=10,
                label = strong("Number of iterations"),
                value = 500),
              p(),
              numericInput("fitCR",
                label = strong("Crossover probability (CR)"),
                min=0, max=1, step=0.01, value = .85),
              p(),
              numericInput("fitF",
                label = strong("Differential weighting factor (F)"),
                min=0, max=2, step=0.01, value = .7),
              p(),
              textInput("fitScale",
                label = strong("Lower and upper bounds for scale factor fit"),
                value = "1, 1"),
              p(),
              uiOutput("bkgBoundsR"),
              p(),  
              uiOutput("fitWithR"),              
              p(),
              conditionalPanel(
                condition = "input.fitWith == 'fitWith.splines'",
                textInput("fitKnots",
                  label = strong("Number of splines or spline knot positions"),
                  value = 20
                )            
              ),
              p(),
              actionButton("doFit", label = strong("Start fit"))
            ),

##   = FIT RESULTS =
            conditionalPanel(
              condition = "input.mainRB == 'save'",
              h3("Fit Results"),
              p(),
              uiOutput("downloadRDataR"),
              p(),
              uiOutput("downloadFitResAsTxtR"),
              p(),
              uiOutput("downloadFixR"),
              p(),
              br(),
#              uiOutput("messageFixR"),
#              p(),
              uiOutput("selectFixR"),
              uiOutput("appendFixR"),
              p(),              
              br(),
              uiOutput("outHeaderGr"),
              p(),
              uiOutput("rminCalcGrR"),
              uiOutput("rmaxCalcGrR"),
              uiOutput("drCalcGrR"),
              uiOutput("calcGrButtonR"),
              p(),
              uiOutput("downloadGrR"),
              p(),
              br(),
              uiOutput("iterHeader"),
              p(),
              uiOutput("iterTechniqueR"),
              uiOutput("iterEpsR"),
              uiOutput("iterNIterR"),
              uiOutput("doIterationR")
            )         
          )
        )
	    ),
      ##   = Plot Options =
      conditionalPanel(
        condition = "input.mainRB == 'plot'",               
        tags$table(style="border: 3px solid #E8E8E8; padding:2px;border-style:solid; width:100%; height: 450px;" ,
          tags$tr(valign="top", 
            tags$td(
              p(),
              uiOutput("selectPlotR"),
              p(),
              uiOutput("youCanSeePlot"),
              p(),
              p(),
              uiOutput("axisLimsTxt"),
              p(),
              uiOutput("plotLimXR"),
              p(),
              uiOutput("plotLimYR"),
              p(),
              actionButton("rescaleY", label = strong("rescale Y slider")),
              actionButton("resetY", label = strong("reset Y slider"))
              
            )
          )
        )
      ) # end of the last conditional panel     
    ), ## end of == sidebar panel ==

############################################################################
##   === MAIN PANEL ===
    mainPanel(    
      tabsetPanel(
        tabPanel("Data Plot",  
                 progressInit(),
                 tags$table(style="border-style:none;  width:100%",
                   tags$tr(
                     tags$td(uiOutput("progress"), align="center", style="width:70%"),
                     tags$td(uiOutput("selectBank"), align="left", style="width:30%")
                      )
                    ), 
                 plotOutput("dataPlot", clickId="mainClick", hoverId="mainHover", hoverDelay=50, hoverDelayType="debounce", width='750px'),
                 plotOutput("legendPlot", height="20px", width='750px'),
                 uiOutput('downloadMainPlotR'),
                 plotOutput("prelimGrPlot"),
                 uiOutput('downloadestGrPlotR'),
                 uiOutput('downloadestGrDataR')),
        tabPanel("Data Table",  
                 downloadButton('downloadData', 'Download'), tags$hr(), 
                 tableOutput('datatable')),
        tabPanel("Help", includeHTML("./help/help.html")),
        tabPanel("Fit Results Plot",  
                 plotOutput("fitResPlot"), 
                 uiOutput('downloadFitResPlotR'), p(), 
                 plotOutput("GrPlot"),
                 uiOutput('downloadGrPlotR'))
      )
    )
  )  ## end of === layout ===

)
)


