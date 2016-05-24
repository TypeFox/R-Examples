options(rgl.useNULL=TRUE)
## CHECKS ##
if(!require("scatterD3")) stop("package scatterD3 is required")
if(!require("shiny")) stop("package shiny is required")
if(!require("rglwidget")) stop("package rglwidget is required")
if(!require("rgl")) stop("package rgl is required")
if(!require("RLumShiny")) stop("package RLumShiny is required")
if(!require("shinyBS")) stop("package shinyBS is required")

## DEFINE UI ##
shinyUI(
  navbarPage("",position="fixed-top", collapsible = TRUE,
             theme = "bootstrap.simplex.css", 
             
             tabPanel("Tree landscape explorer",
                      tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                      tags$style(type="text/css", "body {padding-top: 40px;}"),
                      pageWithSidebar(
                        ##  TITLE ##
                        headerPanel(
                          img(src="img/logo.png", height="160")
                        ),
                        
                        ## SIDE PANEL CONTENT ##
                        sidebarPanel(
                          tags$head(tags$style(
                            type = 'text/css',
                            'form.well { max-height: 1600px; overflow-y: auto; }'
                          )),
                          
                          ## SPECIFIC TO TREE LANDSCAPE EXPLORER ##
                          conditionalPanel(condition = "$('li.active a').first().html()== 'Tree landscape explorer'",
                                           ## INPUT
                                           ## choice of type of data source
                                           img(src="img/line.png", width="100%"),
                                           h2(HTML('<font color="#6C6CC4" size="6"> > Input </font>')),
                                           radioButtons("datatype", HTML('<font size="4"> Choose data source:</font>'),
                                                        list("Example: Dengue fever"="exDengue", 
                                                             "Example: Woodmice"="exWoodmice",
                                                             "Input file"="file")),
                                           
                                           ## choice of dataset if source is a file
                                           conditionalPanel(condition = "input.datatype=='file'",
                                                            fileInput("datafile", p(HTML(' <font size="4"> Choose input file:</font>'), br(),
                                                                                    strong("accepted formats:"), br(),
                                                                                    em("- multiphylo"), "saved from R (.RData/.rda/.rds)", br(),
                                                                                    em("- nexus"), "file (.nex/.nexus)")
                                                            ),
                                                            checkboxInput("randSamp","Randomly sample from the trees?", value=TRUE),
                                                            bsTooltip("randSamp", "For large sets of trees and/or trees with many tips the app may be slow, so beginning the analysis with a small random sample is recommended.", 
                                                                      placement = "right", trigger = "hover", options = NULL),
                                                            conditionalPanel(
                                                              condition="input.randSamp",
                                                              sliderInput("sampleSize", "Size of random sample:", value=10, min=10, max=300, step=10)
                                                            )
                                           ),
                                           
                                           
                                           ## ANALYSIS
                                           img(src="img/line.png", width="100%"),
                                           h4(HTML('<font color="#6C6CC4" size="6"> > Analysis </font>')),
                                           
                                           ## choose metric
                                           selectInput("treemethod", "Choose a tree summary:",
                                                       choices=c(
                                                         "Kendall Colijn" = "metric",
                                                         "Billera, Holmes, Vogtmann" = "BHV",
                                                         "Robinson Foulds (unrooted)" = "RF",
                                                         "Tip-tip path distance (unrooted)" = "nNodes",
                                                         "Tip-tip branch-length distance (unrooted)" = "patristic",
                                                         "Abouheif test" = "Abouheif",
                                                         "Sum of direct descendents" = "sumDD")),
                                           
                                           ## lambda, axes
                                           uiOutput("lambda"),
                                           bsTooltip("lambda","When lambda=0 trees are compared topologically; increasing lambda gives more focus to branch lengths"),
                                           
                                           conditionalPanel(
                                             condition="input.plotType==1",
                                             uiOutput("naxes")
                                             
                                             ## Future: highlight median trees (if plotType==1)
                                             #checkboxInput("showMedians", label=strong("Highlight median tree(s)?"), value=FALSE)
                                           ),
                                           
                                           ## show Shepard plot?
                                           checkboxInput("quality", label=strong("Assess quality of projection (Shepard plot)?"), value=FALSE),
                                           bsTooltip("quality","A Shepard plot gives an indication of the quality of the MDS projection. It will be displayed below the main plot.", placement="right"),
                                           
                                           ## show screeplot?
                                           conditionalPanel(
                                             condition="input.graphics==1",
                                             checkboxInput("scree", label=strong("Show screeplot?"), value=FALSE),
                                             bsTooltip("scree","Display screeplot of the eigenvalues associated with each componenet? It will be displayed below the main plot.", placement="right")
                                           ),
                                           
                                           ## find clusters?
                                           checkboxInput("findClusters", label=strong("Identify clusters?"), value=FALSE),
                                           bsTooltip("findClusters","Statistical tools for choosing an appropriate clustering method and number of clusters will be added to treescape soon.", placement="right"),
                                           
                                           conditionalPanel(condition ="input.findClusters", 
                                                            radioButtons("clusterType", label="Method:", 
                                                                         choices=c("statistically"="stat","by metadata"="meta"), selected="stat"),
                                                            conditionalPanel(
                                                              condition="input.clusterType=='stat'",
                                                              
                                                              ## clustering method
                                                              selectInput("clustmethod", "Clustering method:",
                                                                          choices=c(
                                                                            "Ward" = "ward.D2",
                                                                            "Single" = "single",
                                                                            "Complete" = "complete",
                                                                            "UPGMA" = "average")),
                                                              
                                                              ## number of clusters
                                                              uiOutput("nclust")
                                                            ),
                                                            conditionalPanel(
                                                              condition="input.clusterType=='meta'",
                                                              fileInput("metadatafile", p(HTML(' <font size="4"> Choose input file:</font>'), br(),
                                                                                          strong("accepted formats:"), br(),
                                                                                          em("- object of class factor/numeric/character/list"), "saved from R (.RData/.rda)", br(),
                                                                                          em("- csv file"), "(.csv) (first column will be used)")
                                                              )
                                                            )
                                                            
                                                            
                                           ),
                                           
                                           
                                           
                                           
                                           ## relevant if method = KC metric, allow tip emphasis
                                           conditionalPanel(
                                             condition="input.treemethod=='metric'",
                                             ## Emphasise tips
                                             checkboxInput("emphTips", label=strong("Emphasise tips?"), value=FALSE),
                                             bsTooltip("emphTips","Choose tips to emphasise or de-emphasise: the vector elements corresponding to these tips are multiplied by the weight chosen below.", placement="right"),  
                                             ## if tip emphasis is chosen, provide options:
                                             conditionalPanel(
                                               condition="input.emphTips",
                                               
                                               uiOutput("whichTips"),
                                               
                                               sliderInput("emphWeight", "Weight of emphasis", value=2,min=0.1,max=100)
                                             )
                                           ),
                                           
                                           ## AESTHETICS
                                           img(src="img/line.png", width="100%"),
                                           
                                           h2(HTML('<font color="#6C6CC4" size="6"> > Aesthetics </font>')),
                                           
                                           ## tree landscape or compare to single reference tree
                                           conditionalPanel(
                                             condition="input.plot3D==2",
                                             radioButtons("plotType", "View",
                                                          choices=c("Full tree landscape"=1,"Distances from a reference tree"=2),
                                                          selected=1),
                                             bsTooltip("plotType", "Choose whether to view the relative distances between all trees, or a 1-dimensional plot of their distances from a fixed reference tree")
                                           ),
                                           
                                           
                                           ## Dimensions (3D possible if 3 or more axes retained, and full tree landscape)
                                           conditionalPanel(condition="input.naxes>2",
                                                            conditionalPanel(
                                                              condition="input.plotType==1",
                                                              radioButtons("plot3D", "Plot dimensions",
                                                                           choices=c("2D"=2,"3D"=3),
                                                                           selected=2)
                                                            )
                                           ),
                                           
                                           conditionalPanel(
                                             condition="(input.plot3D==2)&&(input.plotType==1)",
                                             radioButtons("graphics", "Display using",
                                                          choices=c("plotGrovesD3"=1,"plotGroves"=2),
                                                          selected=1),
                                             bsTooltip("graphics", "Choose whether to view the tree landscape using plotGrovesD3 which uses scatterD3 (interactive html) or plotGroves which uses adegraphics")
                                           ),
                                           
                                           # if plotType=1, pick the axes to view:
                                           conditionalPanel(    
                                             condition="input.plotType==1",
                                             ## select first axis to plot
                                             numericInput("xax", "Indicate the x axis", value=1, min=1, max=3),
                                             
                                             ## select second axis to plot
                                             numericInput("yax", "Indicate the y axis", value=2, min=1, max=3),
                                             bsTooltip("yax", "If multiple MDS axes have been retained, any combination of axes can be viewed"),
                                             
                                             ## if in 3D, need a z axis:
                                             conditionalPanel(condition="input.plot3D==3",
                                                              numericInput("zax", "Indicate the z axis", value=3, min=1, max=3)
                                             )
                                           ),
                                           
                                           ## aesthetics for tree landscape view
                                           conditionalPanel(
                                             condition="input.plotType==1",
                                              conditionalPanel(
                                                condition="(input.plot3D==2)&&(input.graphics==1)",
                                             
                                                   ## Animate transitions?
                                                  checkboxInput("transitions", label="Animate transitions?", value=TRUE)
                                              ),
                                           
                                              conditionalPanel(
                                                condition="(input.plot3D==2)&&(input.graphics==2)",
                                             
                                                 ## convex hulls or ellipses when clusters identified
                                                    conditionalPanel(
                                                      condition="input.findClusters",
                                                      radioButtons("scattertype", "Type of scatterplot",
                                                                   choices=c("chull","ellipse"),
                                                                   selected="chull")
                                                    ),
                                                
                                                selectInput("screemds", "Position of the MDS screeplot:",
                                                         choices=c("None" = "none",
                                                                   "Bottom right" = "bottomright",
                                                                   "Bottom left" = "bottomleft",
                                                                   "Top right" = "topright",
                                                                   "Top left" = "topleft"),
                                                         selected="bottomleft")
                                              ),
                                          
                                             ## symbol size
                                             sliderInput("pointsize", "Size of the points", value=2, min=0, max=10, step=0.2),
                                             
                                             conditionalPanel(
                                               condition="(input.plot3D==2)&&(input.graphics==1)",
                                               ## symbol size
                                               sliderInput("pointopacity", "Opacity of the points", value=0.6, min=0, max=1, step=0.05)
                                             ),
                                             
                                             conditionalPanel(
                                               condition="input.plot3D==2",
                                                  ## display labels
                                                  checkboxInput("showlabels", label="Display tree labels?", value=FALSE),
                                             
                                                  conditionalPanel(
                                                      condition="input.showlabels",
                                                      ## label size
                                                      sliderInput("labelsize", "Size of the labels", value=1, min=0, max=10, step=1),
                                               
                                                      conditionalPanel(
                                                          condition="input.graphics==2",
                                                          checkboxInput("optimlabels", label="Optimize label position?", value=FALSE)
                                                      )  
                                                  )
                                             )
                                             
                                          ),
                                           
                              
                                              
                                           
                                          # if plotType=2, option to stretch
                                           conditionalPanel(
                                             condition="input.plotType==2",
                                             
                                             uiOutput("selectedRefTree"),
                                             
                                             sliderInput("stretch", "Height of plot (pixels)", value=1600, min=800, max=12800, step=200)
                                             
                                           ),
                                           
                                           
                                           ## choose color palette (if clusters detected)
                                           conditionalPanel(
                                             ## condition
                                             condition="input.findClusters",
                                             
                                             selectInput("palette", "Palette for the clusters",
                                                         choices=c("funky", "spectral",
                                                                   "seasun", "lightseasun", "deepseasun",
                                                                   "rainbow", "azur", "wasp"),
                                                         selected="funky")
                                           ),
                                           
                                           conditionalPanel(
                                             condition="input.plotType==1",

                                             ## choose label colors
                                             jscolorInput("labcol", "Label / point color", value="#1B2266", close=TRUE)
                                             
                                           ),
                                           
                                           
                                           
                                           
                                           br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                           br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                           width=4)), # end conditional panel and sidebarPanel; width is out of 12
                        ## MAIN PANEL
                        mainPanel("",
                                  
                                  # TITLE #
                                  h2(HTML('<font color="#6C6CC4" size="6"> Tree landscape explorer </font>')),
                                  br(),br(),
                                  
                                  ## function I was using for testing:
                                  #verbatimTextOutput("plot_click"),
                                  
                                  # Removed:
                                  #verbatimTextOutput("caption"),
                                  
                                  conditionalPanel(
                                    condition="input.plot3D==2",
                                    uiOutput("treescapePlot")
                                  ),
                                  conditionalPanel(
                                    condition="input.plot3D==3",
                                    rglwidgetOutput("treescapePlot3D", width="800px")
                                  ),
                                  
                                  conditionalPanel(
                                    condition="(input.plotType==1)&&(input.plot3D==2)&&(input.graphics==1)",
                                      tags$p(actionButton("scatterD3-reset-zoom", HTML("<span class='glyphicon glyphicon-search' aria-hidden='true'></span> Reset Zoom")))
                                  ),
                                  
                                  conditionalPanel(
                                    condition="input.quality",
                                    uiOutput("shep")
                                  ),
                                  
                                  conditionalPanel(
                                    condition="(input.scree)&&(input.graphics==1)",
                                    uiOutput("scree")
                                  ),
                                  
                                  br(), br(),
                                  
                                  ## OUTPUT (save)
                                  img(src="img/line.png", width="400px"),
                                  h2(HTML('<font color="#6C6CC4" size="6"> > Output </font>')),
                                  
                                  ## save MDS plot
                                  conditionalPanel(
                                    condition="(input.plotType==1)&&(input.plot3D==2)&&(input.graphics==1)",
                                        tags$p(tags$a(id = "scatterD3-svg-export", href = "#",
                                              class = "btn btn-default", HTML("<span class='glyphicon glyphicon-save' aria-hidden='true'></span> Save treescape plot as svg"))),
                                        downloadButton("downloadMDS2Dhtml", "Save treescape plot as interactive html")
                                    ),
                                  
                                  conditionalPanel(
                                    condition="(input.plotType==1)&&(input.plot3D==2)&&(input.graphics==2)",
                                    downloadButton("downloadMDS", "Save treescape image as png file")
                                  ),
                                  
                                  conditionalPanel(
                                    condition="input.plot3D==3",
                                    downloadButton("downloadMDS3Dhtml", "Save treescape 3D plot as interactive html")
                                  ),
                                  
                                  conditionalPanel(
                                    condition="input.quality",
                                    downloadButton("downloadShep", "Save Shepard plot as png file")
                                  ),
                                  
                                  conditionalPanel(
                                    condition="input.scree",
                                    downloadButton("downloadScree", "Save screeplot as png file")
                                  ),
                                  
                                  ## save trees to nexus file
                                  downloadButton('exporttrees', "Save trees to nexus file"),
                                  
                                  ## save results to csv
                                  downloadButton('exportrestocsv', "Save results (MDS+clusters) to csv file"),
                                  
                                  ## save results to RData
                                  downloadButton('exportrestordata', "Save results (MDS+clusters) to R object")
                        ) # end mainPanel
                      ) # end page with sidebar
             ), # end tabPanel
             tabPanel("Tree viewer",
                      tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                      tags$style(type="text/css", "body {padding-top: 40px;}"),
                      pageWithSidebar(
                        ##  TITLE ##
                        headerPanel(
                          img(src="img/logo.png", height="160")
                        ),
                        
                        ## SIDE PANEL CONTENT ##
                        sidebarPanel(
                          tags$head(tags$style(
                            type = 'text/css',
                            'form.well { max-height: 1600px; overflow-y: auto; }'
                          )),
                          
                          ## INPUT
                          ## choice of tree type
                          img(src="img/line.png", width="100%"),
                          h2(HTML('<font color="#6C6CC4" size="6"> > Input </font>')),
                          
                          radioButtons("treeChoice", "Tree selection",
                                       choices=c("Median tree"="med","General tree selection"="gen"),
                                       selected="med", width="100%"),
                          bsTooltip("treeChoice", "A geometric median tree is plotted by default. If clusters have been identified, the median for each can be viewed. Alternatively, any individual tree can be plotted."),
                          
                          conditionalPanel(condition = "input.treeChoice=='med'",
                                           selectInput("selectedMedTree", "Median tree from:", 
                                                       choices=c("All trees"="all"))
                          ),
                          
                          conditionalPanel(condition = "input.treeChoice=='gen'",
                                           uiOutput("selectedGenTree")
                          ),
                          
                          
                          ## TREE AESTHETICS
                          img(src="img/line.png", width="100%"),
                          h2(HTML('<font color="#6C6CC4" size="6"> > Aesthetics </font>')),
                          
                          ## condition on tree being displayed
                          conditionalPanel(condition = "input.selectedTree!=''",
                                           ## use edge lengths?
                                           checkboxInput("edgelengths", label="Use original branch lengths?", value=TRUE),
                                           
                                           ## ladderize
                                           checkboxInput("ladderize", label="Ladderize the tree?", value=TRUE),
                                           
                                           ## type of tree
                                           radioButtons("treetype", "Type of tree",
                                                        choices=c("phylogram","cladogram", "fan", "unrooted", "radial"),
                                                        selected="phylogram", width="100%"),
                                           
                                           ## tree direction
                                           radioButtons("treedirection", "Direction of the tree",
                                                        choices=c("rightwards", "leftwards", "upwards", "downwards"),
                                                        selected="rightwards", width="100%"),
                                           
                                           ## tip labels
                                           checkboxInput("showtiplabels", label="Display tip labels?", value=TRUE),
                                           
                                           conditionalPanel(condition="input.showtiplabels",
                                                            ## tip label font
                                                            selectInput("tiplabelfont", "Tip label font", 
                                                                        choices=c("Plain"="1","Bold"="2","Italic"="3","Bold italic"="4")),
                                                            
                                                            ## tip label size
                                                            sliderInput("tiplabelsize", "Size of the tip labels", value=1, min=0, max=5, step=0.1)
                                                            
                                           ),
                                           
                                           
                                           ## edge width
                                           sliderInput("edgewidth", "Width of the edges", value=2, min=1, max=20, step=0.2),
                                           
                                           ## edge colour
                                           jscolorInput("edgecolor", "Edge colour", value="#000000", close=TRUE),
                                           
                                           
                                           br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                           br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                           width=4)), # end conditional panel and sidebarPanel; width is out of 12
                        
                        ## MAIN PANEL
                        mainPanel("",
                                  # TITLE #
                                  h2(HTML('<font color="#6C6CC4" size="6"> Tree viewer </font>')),
                                  br(),br(),
                                  
                                  ## conditional panel: plot tree if needed
                                  conditionalPanel(condition = "input.selectedTree!=''",
                                                   plotOutput("tree", height = "800px"),
                                                   
                                                   br(), br(),
                                                   
                                                   ## OUTPUT (save)
                                                   img(src="img/line.png", width="400px"),
                                                   h2(HTML('<font color="#6C6CC4" size="6"> > Output </font>')),
                                                   downloadButton("downloadTree", "Save tree image"),
                                                   
                                                   br(), br(), br(), br(), br(), br()
                                  ) # end single tree conditional panel  
                        ) # end mainPanel
                      ) # end page with sidebar
             ), # end tabPanel "Tree Viewer"
             tabPanel("densiTree viewer",
                      tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                      tags$style(type="text/css", "body {padding-top: 40px;}"),
                      pageWithSidebar(
                        ##  TITLE ##
                        headerPanel(
                          img(src="img/logo.png", height="160")
                        ),
                        
                        ## SIDE PANEL CONTENT ##
                        sidebarPanel(
                          tags$head(tags$style(
                            type = 'text/css',
                            'form.well { max-height: 1600px; overflow-y: auto; }'
                          )),
                          
                          ## INPUT
                          ## choice of tree type
                          img(src="img/line.png", width="100%"),
                          h2(HTML('<font color="#6C6CC4" size="6"> > Input </font>')),
                          
                          ## add densiTree selector (gets updated to number of clusters by )
                          selectInput("selectedDensiTree", "Choose collection of trees to view in densiTree plot", 
                                      choices=c("Choose one"="","All trees"="all"), width="100%"),
                          #h2(HTML('<font color="#6C6CC4" size="2"> Note: this can be slow for large sets of trees </font>')),
                          
                          bsTooltip("selectedDensiTree", "View all trees together in a densiTree plot. If clusters have been identified, the set of trees from a single cluster can be plotted. Note this function can be slow if many trees are included.", placement="bottom"),
                          
                          
                          ## DENSITREE AESTHETICS
                          img(src="img/line.png", width="100%"),
                          h2(HTML('<font color="#6C6CC4" size="6"> > Aesthetics </font>')),
                          
                          conditionalPanel(condition = "input.selectedDensiTree!=''",
                                           
                                           ## alpha (semitransparency of edges)
                                           sliderInput("alpha", "Transparency of edges", value=0.5, min=0, max=1, step=0.05),
                                           
                                           
                                           checkboxInput("scaleX", label="Scale trees to equal heights?", value=FALSE),                 
                                           
                                           br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                           br(),br(),br(),br(),br(),br(),br(), # add some blank space at the end of side panel
                                           width=4)), # end conditional panel and sidebarPanel; width is out of 12
                        
                        ## MAIN PANEL
                        mainPanel("",
                                  
                                  # TITLE #
                                  h2(HTML('<font color="#6C6CC4" size="6"> densiTree viewer </font>')),
                                  br(),br(),
                                  
                                  ## conditional panel: plot tree if needed
                                  conditionalPanel(condition = "input.selectedDensiTree!=''",
                                                   
                                                   plotOutput("densiTree", height = "800px"),
                                                   
                                                   br(), br(),
                                                   
                                                   ## OUTPUT (save)
                                                   img(src="img/line.png", width="400px"),
                                                   h2(HTML('<font color="#6C6CC4" size="6"> > Output </font>')),
                                                   downloadButton("downloadDensiTree", "Save densiTree image"),
                                                   
                                                   br(), br(), br(), br(), br(), br()
                                  ) # end densiTree conditional panel
                        ) # end of main panel "Multi-tree viewer"
                      ) # end of page with sidebar
             ),# end of tabPanel "densiTree viewer"
             ## HELP SECTION
             tabPanel("Help",
                      tags$style(type="text/css", "body {padding-top: 40px;}"),
                      HTML(paste(readLines("www/html/help.html"), collapse=" "))
             ),
             
             ## SERVER INFO ##
             tabPanel("System info", 
                      tags$style(type="text/css", "body {padding-top: 40px;}"),
                      verbatimTextOutput("systeminfo"))
  ) # end of tabsetPanel
) # end of Shiny UI