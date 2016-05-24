
library(statnetWeb)

# Everything that gets displayed inside the app is enclosed in a call to `shinyUI`.
# The first thing to be specified is the type of page to display. The `navbarPage`
# includes a navigation bar at the top of the page and each tab leads to different
# pages of content.

shinyUI(
  navbarPage(
    #theme="mycosmo.css",
    title=NULL,
    id= 'navbar', windowTitle = 'statnetWeb', collapsible=TRUE,



# Front Page (About) ------------------------------------------------------


tabPanel(title=span('statnetWeb', id="sWtitle"),
         value='tab1',
         fluidRow(
          column(2,
                 actionButton("aboutButton", label = "About statnetWeb",
                              class = "btn active"),
                 actionButton("citeButton", label = "Citing statnetWeb",
                              class = "btn"),
                 actionButton('startButton', label='Get Started',
                              class="btn btn-primary")
          ),
   column(6, style="padding: 0 30px 0 0;",
          div(id="aboutbox",
            p("Welcome to the first interactive interface for the", strong("ergm"),
              "package!", strong("ergm"), "is part of the statnet network analysis software -",
              "a suite of packages written in R - and this GUI also includes some of the functionality",
              "from the associated packages", strong("network"), " and ", strong("sna"), ".  This web application",
              "is written with the Shiny framework from RStudio and development is via GitHub.  More information",
              "on the statnet software, the ergm package, Shiny and our GitHub repository can be found in the",
              "resource links on the right."),

            p("This interface is useful for teachers and students of introductory network analysis,",
              "for newcomers to exponential random graphs models, and for experienced network modelers",
              "who want easier access to analysis results. If you are new to ergm, you may find it helpful",
              "to work through the", a("ergm tutorial",
                                       href="http://statnet.csde.washington.edu/workshops/SUNBELT/EUSN/ergm/ergm_tutorial.html",
                                       target="_blank"), "using this interface. Advanced users will still want to interact",
              "via the command line in order to access the full functionality of ergm."),
            p("A typical network analysis will move sequentially through the tabs at the top of the page.",
              "Click on the help icon at the top of any page for guidance."),
            p("Do you have comments/suggestions/complaints on this prototype app? Please share them with us.",
              "They are best submitted through our", a('GitHub site,',
                                                       href='https://github.com/statnet/statnetWeb',
                                                       target='_blank'),
              "or by email to the statnet_help mailing list (see", actionLink("helpLink", "Help"), "tab).")
          ),
          div(id="citebox",
            tabsetPanel(
              tabPanel("BibTeX",
p(strong("statnet")),
tags$pre(id='scitation','@Manual{handcock:statnet,
  title = {statnet: Software tools for the Statistical Modeling of Network Data},
  author = {Mark S. Handcock and David R. Hunter and Carter T. Butts and Steven M. Goodreau and Martina Morris},
  year = {2003},
  address = {Seattle, WA},
  url = {http://statnetproject.org}
}'),

p(strong("statnetWeb")),
tags$pre(id='swcitation',"@Manual{beylerian:statnetWeb,
  title = {\\pkg{statnetWeb}: A Graphical User Interface for Network Modeling with 'Statnet'},
  author = {Emily N. Beylerian and Samuel Jenness and Kirk Li and Martina Morris},
  year = {2015},
  note = {\\proglang{R}~package version~0.3.4},
  address = {Seattle, WA},
  url = {https://cran.r-project.org/web/packages/statnetWeb/}
}")
                       ),
              tabPanel("Other",
p(strong("statnet")),
tags$pre("Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau, and
Martina Morris (2003). statnet: Software tools for the Statistical Modeling
of Network Data. URL http://statnetproject.org"),

p(strong("statnetWeb")),
tags$pre("Emily N. Beylerian, Samuel Jenness, Kirk Li, and Martina Morris (2014).
statnetWeb: A Graphical User Interface for Network Modeling with 'Statnet'.")
                       )
            ),

            p('If you use statnet or statnetWeb, please cite them.',
              'Additional citation information for statnet',
              'and the component packages can be found here:'),
            tags$ul(
              tags$li(a('Citing statnet',
                        href='https://statnet.csde.washington.edu/trac/wiki/citation%20information',
                        target='_blank')),
              tags$li(a('License and source code attribution requirements',
                        href = 'http://statnet.csde.washington.edu/attribution.shtml',
                        target = '_blank')),
              tags$li(a('statnet Development Team',
                        href = 'http://statnet.csde.washington.edu/about_us.shtml',
                        target = '_blank'))
            )
            )
          ),
   column(4,
          wellPanel(
              h5(tags$u('Resources')),
              div(title = "Wiki page for statnetWeb",
                a("About statnetWeb",
                  href = "https://statnet.csde.washington.edu/trac/wiki/statnetWeb",
                  target = "_blank")),
              div(title=paste("Homepage of the statnet project with tutorials,",
                              "publications and recent news."),
                  a("About statnet",
                    href = "https://statnet.csde.washington.edu/trac", target = "_blank")
              ),

              column(11, offset = 1,
                    span(id="linktitle1",'Key background papers',icon('angle-double-left')),br(),
                    div(id="linkbox1",
                      a("ergm: Journal of Statistical Software",
                        href = "http://www.jstatsoft.org/v24/i03/", target = "_blank"),
                      br(),
                      a("Using ergm: Journal of Statistical Software",
                        href = "http://www.jstatsoft.org/v24/i04/", target = "_blank")),

                    span(id="linktitle2",'Tutorials and documentation',icon('angle-double-left')),br(),
                    div(id="linkbox2",
                        a("ergm tutorial from Sunbelt EUSN 2014 Workshop",
                        href = "http://statnet.csde.washington.edu/workshops/SUNBELT/EUSN/ergm/ergm_tutorial.html",
                        target= "_blank"),
                      br(),
                      a("ergm documentation on CRAN",
                        href = "http://cran.r-project.org/web/packages/ergm/ergm.pdf",
                        target = "_blank")),
                    style="margin-bottom:10px;"),
              br(),
              div(a("statnetWeb on GitHub", href="https://github.com/statnet/statnetWeb",
                    target="_blank")),
              div(a("Shiny: a web application framework for R", href="http://shiny.rstudio.com/",
                    target="_blank"))
   ),
   fluidRow(img(src= 'UW.Wordmark_ctr_K.jpg', width=200), style="margin-left:15px;"),
   fluidRow(a(img(src = 'csdelogo_crop.png', height = 40, width = 40),
             href = 'https://csde.washington.edu/', target = '_blank'),
            a(img(src = 'csde_goudy.fw.png', width=150), href = 'https://csde.washington.edu/',
             target = '_blank'), style="margin-left:15px;")
   )
   )
 ),


# Data Upload -------------------------------------------------------------


# Before the code for what is displayed on the Data Upload page,
# various javaScript and CSS files that will be useful later in the
# script are linked. For example, since network plotting and model
# fitting do not happen instantly (especially for large networks),
# a loading icon will help to assure users that the app is still working
# on producing output. The file busy.js controls the behavior of the
# loading message and style.css controls the appearance. To display
# the loading message on subsequent tabs, we only need to include the
# div statement within those tabs.

tabPanel(title='Data', value='tab2',
         #busy.js is for calculation in progress boxes
         #alert.js is for popup boxes,
         #jquery libraries are loaded from google cdn, needed for autocomplete
         #this tagList command has to go inside a tabPanel
         tagList(
           tags$head(
             tags$link(rel="stylesheet", type="text/css",href="style.css"),
             #tags$link(rel="stylesheet", type="text/css",href="autocomplete.css"),
             tags$link(rel="stylesheet", type="text/css",
                       href="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.1/themes/smoothness/jquery-ui.css"),
             tags$script(type="text/javascript", src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"),
             #tags$script(type="text/javascript", src="autocomplete.js"),
             tags$script(type="text/javascript", src="//ajax.googleapis.com/ajax/libs/jqueryui/1.11.1/jquery-ui.min.js"),
             tags$script(type="text/javascript", src="busy.js"),
             tags$script(type="text/javascript", src="alert.js"),
             tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",
                              function(message) {
                              console.log(message)
                              eval(message.code);
                              });'))
           )
         ),

# Conditional panels are only displayed when a specified condition is true.
# The condition is a javascript expression that can refer to the current
# values of input or output objects. When the condition is false, the panel
# does not take up any space in the UI.


fluidRow(
  column(7,
    tabsetPanel(id='datatabs',
      tabPanel('Upload Network', br(),
         wellPanel(
           fluidRow(
             column(6,
                    selectInput('filetype',label='File type',
                                 choices=c('built-in network'= 5,
                                           'statnet network object (*.rds)' = 1,
                                           'matrix of relational data (*.csv or *.rds)' = 4,
                                           'Pajek network (*.net)' = 2,
                                           'Pajek project (*.paj)' = 3))
                    ),
             conditionalPanel(condition = 'input.filetype < 5',
               column(6,
                    br(),
                    fileInput(inputId='rawdatafile', label=NULL, accept='text'),
                    verbatimTextOutput('rawdatafile'))
                ),
             conditionalPanel(condition = 'input.filetype == 5',
                column(6,
                    br(style="line-height:26px;"),
                    selectizeInput('samplenet', label=NULL,
                                choices=c("Choose a network" = '', 'ecoli1', 'ecoli2',
                                          'faux.mesa.high','flobusiness',
                                          'flomarriage', 'kapferer', 'kapferer2',
                                          'molecule', 'samplike', 'samplk1',
                                          'samplk2', 'samplk3'))
                )
               )
             ),
           fluidRow(
             conditionalPanel(condition='input.filetype == 4',
                 column(1, align="right",
                        style="margin-top:5px; margin-left:0px;",
                        br(),br(),
                        span(style="line-height:25px;", class="helper",
                            span(id="filetypehelper1",
                                 icon('question-circle'),
                                 div(id="filetypebox1", class="smallhelperbox",
                                     "For adjacency matrices,",
                                     "the first row and column of .csv files",
                                     "should hold vertex labels.")),
                            br(),
                            span(id="filetypehelper2",
                                 icon('question-circle'),
                                 div(id="filetypebox2", class="smallhelperbox",
                                     "For adjacency matrices,",
                                     "the first row and column of .csv files",
                                     "should hold vertex labels.")),
                            br(),
                            span(id="filetypehelper3",
                                 icon('question-circle'),
                                 div(id="filetypebox3", class="smallhelperbox",
                                     "For incidence matrices,",
                                     "the first row of .csv files should hold",
                                     "edge labels, the first column should hold vertex labels")),
                            br(),
                            span(id="filetypehelper4",
                                 icon('question-circle'),
                                 div(id="filetypebox4", class="smallhelperbox",
                                     "For edge lists, .csv files should not have row or column labels.")),

                            br())),
                 column(4,
                        br(),
                        radioButtons('matrixtype', label='Matrix Type',
                                     choices=c('Adjacency matrix'='adjacency',
                                               'Bipartite adjacency matrix'='bipartite',
                                               'Incidence matrix' = 'incidence',
                                               'Edge list' = 'edgelist'))),

                 column(5,
                        br(),
                        span(strong('Network Attributes')),
                        checkboxInput('dir', 'directed?', value=TRUE),
                        checkboxInput('loops', 'loops?', value=FALSE),
                        checkboxInput('multiple', 'multiple?', value=FALSE),
                        checkboxInput('bipartite', 'bipartite?', value=FALSE))
             ),
             conditionalPanel(condition='input.filetype == 3',
                 column(6,
                        uiOutput('pajchooser'))),
             conditionalPanel(condition='input.filetype == 1',
                              p(class="helper", id="Robjhelp", icon("question-circle"), span("What is an .rds file?", style="font-size:0.85em;")),
                              div(class="mischelperbox", id="Robjbox", 'When working in R, an object in your environment',
                                  'can be saved to a .rds file from the command line in the following way:',
                                  code('saveRDS(objectname, file="newfilename.rds")'),br(),'By default the file will be saved',
                                  'into the current working directory. The full path to a new location can be',
                                  'specified in the ', code('file='), 'argument, or set', code('file=file.choose(new=TRUE)'),
                                  'to use a save dialog box.')
                              )
           )),
         conditionalPanel(
           condition="input.filetype == 5 & input.samplenet != ''",
           wellPanel(uiOutput("datadesc"))
           )
         ),
    tabPanel('Edit Network', br(),
         wellPanel(
           fluidRow(

             column(6,strong('Symmetrize'),
                conditionalPanel(condition="output.nwsum != 'NA'",
                    br(),
                    selectizeInput('symmetrize', label=NULL,
                                 choices=c('Do not symmetrize',
                                           'upper: Copy upper triangle over lower'='upper',
                                           'lower: Copy lower triangle over upper'='lower',
                                           'strong: Intersection of ties'='strong',
                                           'weak: Union of ties'='weak')),
                    conditionalPanel(condition="input.symmetrize != 'Do not symmetrize'",
                                     p("After symmetrizing, network should be:"),
                                     actionButton("symmdir", "directed", class="btn-sm"),
                                     actionButton("symmundir", "undirected", class="btn-sm active")
                                     )

                    )),

             column(5,strong('Import new attribute information'),
                    conditionalPanel(condition="output.nwsum != 'NA'",
                       br(),
                       selectizeInput("newattrtype", label=NULL,
                                 choices=c("Vertex attributes" = "vertexattr",
                                           "Vertex names" = "vertexnames",
                                           "Edge attributes" = "edgeattr"),
                                 options=list(placeholder="Select attribute type",
                                              onInitialize = I('function() { this.setValue(""); }')
                                              )),
                       conditionalPanel(condition="input.newattrtype == 'edgeattr'",
                            radioButtons("edgeform", label = NULL,
                                         choices = c("Attributes are in adjacency matrix form" = "matrix",
                                                     "Attributes are in vector form" = "vector")),

                         conditionalPanel(condition="input.edgeform == 'matrix'",
                              span('Upload a file of one of the following types:',br(),
                              tags$ul(
                                tags$li('.rds file',
                                        span(class="helper",id="filetypehelper5",
                                             icon("question-circle"),
                                             div(id="filetypebox5", class="mischelperbox",
                                                 strong("Edge values"), "should be in adjacency matrix form",
                                                 "and saved into a list object in R. For example,",br(),
                                                 code("mylist <- list()"),br(),
                                                 code("mylist$eval1 <- matrix(...)"), br(),
                                                 code("mylist$eval2 <- matrix(...)"), br(),
                                                 "The named elements of the list ('eval1' and 'eval2') will",
                                                 "become the names of the edge attributes in statnetWeb.",br(),br(),
                                                 strong(".rds files"), "can be saved with",
                                                 br(),
                                                 code('saveRDS(objectname, file="newfilename.rds")'),br(),br(),
                                                 "Multiple edge value matrices can be uploaded from one list object."
                                                 ))),
                                tags$li('.csv file',
                                        span(class="helper", id="filetypehelper6",
                                             icon("question-circle"),
                                             div(id="filetypebox6", class="mischelperbox",
                                                 strong(".csv files"), "should include vertex labels in the first",
                                                 "row and column of the matrix. The edge attribute name will be taken",
                                                 "from the filename.",br(),br(),
                                                 "Only one edge value matrix can be uploaded at a time.")))))
                                          ),
                         conditionalPanel(condition="input.edgeform == 'vector'",
                              span('Upload a file of one of the following types:',br(),
                                   tags$ul(
                                     tags$li('.rds file',
                                             span(class="helper",id="filetypehelper9",
                                                  icon("question-circle"),
                                                  div(id="filetypebox9", class="mischelperbox",
                                                      strong("Edge attributes"), "should be in vector form",
                                                      "and saved into a list object in R. For example,",br(),
                                                      code("mylist <- list()"),br(),
                                                      code("mylist$eval1 <- c()"), br(),
                                                      code("mylist$eval2 <- matrix(...)"), br(),
                                                      "The named elements of the list ('eval1' and 'eval2') will",
                                                      "become the names of the edge attributes in statnetWeb.",
                                                      "The values in each vector should be in the same order as the",
                                                      "edge IDs they will be applied to.",
                                                      br(),br(),
                                                      strong(".rds files"), "can be saved with",
                                                      br(),
                                                      code('saveRDS(objectname, file="newfilename.rds")'),br(),br(),
                                                      "Multiple edge attributes can be uploaded from one list object."
                                                  ))),
                                     tags$li('.csv file',
                                             span(class="helper", id="filetypehelper10",
                                                  icon("question-circle"),
                                                  div(id="filetypebox10", class="mischelperbox",
                                                      strong(".csv files"), "should include headers in the first row.",
                                                      "The header of each column will become an edge attribute name.",
                                                      "The values in each column should be in the same order as the",
                                                      "edge IDs they will be applied to.")))))
                                          )
                       ),
                       conditionalPanel(condition="input.newattrtype != '' & input.newattrtype != 'edgeattr'",
                            span('Upload a file of one of the following types:',br(),
                                 tags$ul(
                                   tags$li('.rds file',
                                           span(class="helper",id="filetypehelper7",
                                                icon("question-circle"),
                                                div(id="filetypebox7", class="mischelperbox",
                                                    "Each attribute should be a vector element of an R list object.",
                                                    strong("R lists"), "should be named, for example:",
                                                    br(),
                                                    code("mylist <- list()"),br(),
                                                    code("mylist$age <- c(18,27,20)"),br(),
                                                    code("mylist$sex <- c('M','F','F')"), br(),
                                                    "The names of the elements in the list (e.g. 'age'",
                                                    "and 'sex' in the list above) will become the attribute names.",
                                                    "The values in each vector should be in the same order as the",
                                                    "vertex IDs they will be applied to.",
                                                    br(),br(),
                                                    tags$u("Note:"),'Attributes uploaded as vertex names',
                                                    "will automatically be saved into the", code("vertex.names"),
                                                    "attribute and the names of the list will be ignored.",
                                                    br(),br(),
                                                    strong(".rds files"), "can be saved with",
                                                    br(),
                                                    code('saveRDS(objectname, file="newfilename.rds")')
                                                    ))),
                                   tags$li('.csv file',
                                           span(class="helper", id="filetypehelper8",
                                                icon("question-circle"),
                                                div(id="filetypebox8", class="mischelperbox",
                                                    strong(".csv files"), "should include headers in the first row.",
                                                    "The header of each column will become an attribute name.",
                                                    "The values in each column should be in the same order as the",
                                                    "vertex IDs they will be applied to.",
                                                    br(),br(),tags$u("Note:"),"Attributes uploaded as vertex names",
                                                    "will automatically be saved into the", code("vertex.names"),
                                                    "attribute and the names in the .csv file will be ignored.")))))
                                        ),
                    conditionalPanel(condition="input.newattrtype != ''",
                          fileInput(inputId='newattrvalue', label=NULL),
                          p('New attribute name(s):'),
                          verbatimTextOutput('newattrname'),
                          actionButton('newattrButton', label='Set Attribute', class="btn-sm")
                                        )

                       )

                    )

           )
           )
         )
# tabPanel('Modify Attributes', br(),
#          wellPanel(
#            p('In the future we will build in functions that will ',
#              'allow you to modify the attributes of your network.',
#              'This will include options like:'),
#            tags$ul(
#              tags$li('Applying a function (e.g.', code('sqrt()'), ') to an attribute'),
#              tags$li('Recoding (mapping a set of attributes onto a new set)'),
#              tags$li('Conditional transformations (', code('do...if...'),')'))
#            #uiOutput('modifyattrchooser')
#            )
#          )
  )
),

column(4,
tabsetPanel(
  tabPanel('Network Summary', br(),
           verbatimTextOutput('nwsum')
           ))
  )
),

icon('question-circle', class='fa-2x helper-btn'),
div(class="helper-box", style="display:none",
    p('Upload a file of observed network data (must be of a supported type).',
    'Add custom attributes or symmetrize on the "Edit Network" tab.')),
actionLink('dataleft', icon=icon('arrow-left', class='fa-2x'), label=NULL),
actionLink('dataright', icon=icon('arrow-right', class='fa-2x'), label=NULL)
),


# Network Descriptives ----------------------------------------------------

# There are no calls to selectInput for the options to color code or size the nodes,
# even though they appear in the app. Most widget functions are called in ui.R, but
# this means that all the options passed to them must be static. If the options depend
# on user input (the coloring and sizing menus depend on which network the user
# selects), the widget must be rendered in server.R and output in ui.R with
# iuOutput.

tabPanel(title='Network Descriptives', value='tab3',
 #include progress box when this tab is loading
 div(class = "busy",
     p("Calculation in progress..."),
     img(src="ajax-loader.gif")
 ),

fluidRow(
 column(7,
    tabsetPanel(id='plottabs',
      tabPanel('Network Plot', br(),
                plotOutput('nwplot', click = "plot_click",
                           dblclick = dblclickOpts(id = "plot_dblclick"),
                           hover = hoverOpts(id = "plot_hover", delay = 100,
                                             delayType = "throttle"),
                           brush = brushOpts(id = "plot_brush")
                           )
        ),
      tabPanel('Attributes', br(),
               conditionalPanel('input.attrview == "Large table"',
                                dataTableOutput("attrtbl_lg")
                                ),
               conditionalPanel('input.attrview == "Small tables"',
                                verbatimTextOutput("attrtbl_sm")
                                ),
               conditionalPanel('input.attrview == "Plot summaries"',
                                tags$label("Type of plots"),
                                helpText("Density plots will only be created for",
                                         "numeric attributes with more than nine",
                                         "levels."),
                                selectInput("attrhistaxis",
                                            label = NULL,
                                            choices = c("Barplot: counts" = "count",
                                                        "Barplot: percents" = "percent",
                                                        "Density plot" = "density")),
                                uiOutput("attrhistplotspace"))

               ),
      tabPanel('Degree Distribution',
               p(class='helper', id='ddhelper', icon('question-circle')),
               div(class='mischelperbox', id='ddhelperbox',
                   "Degrees are a node level measure for the number of edges incident on each node.",
                   "The degree distribution shows how the edges in a network are distributed among the",
                   "nodes. The amount (or proportion) of nodes with low, medium",
                   "or high degrees contribute to the overall structure of the",
                   "network. The degree distributions of directed graphs can",
                   "be subset by in-degree or out-degree."),
               plotOutput('degreedist')
               ),
      tabPanel('Geodesic Distribution',
               p(class='helper', id='gdhelper', icon('question-circle')),
               div(class='mischelperbox', id='gdhelperbox',
                   "Geodesics are a dyad level measure for the shortest possible path between a pair of nodes.",
                   'If there is no path between a pair of nodes, the geodesic distance is "inf".',
                   "The geodesic distribution among all possible dyads contributes to the",
                   "structure of network connectivity."),
               plotOutput('geodistplot')
               ),
      tabPanel('More', value='More', br(),
               h5('Conditional uniform graph tests', icon('angle-double-left'),
                  id="cugtitle"),
               wellPanel(id="cugbox",
                 column(4, uiOutput("dynamiccugterm")),
                 column(4, selectInput("ncugsims",
                                       label = "Number of simulations",
                                       choices = c(100, 200, 500))),
                 column(3, actionButton("cugButton", label = "Run",
                                        style="margin-top: 25px;")),
                 br(),
                 plotOutput("cugtest"),
                 br(),
                 downloadButton('cugtestdownload', label = "Download Plot",
                                class="btn-sm")
               ),
               h5('Mixing matrix', icon('angle-double-left'),
                  id="mixmxtitle"),
               wellPanel(id="mixmxbox",
                 fluidRow(
                   column(6, uiOutput('mixmxchooser')),
                   column(6, downloadButton("mixmxdownload",
                                            class = "shiftdown25"))
                 ),
                 fluidRow(
                   verbatimTextOutput('mixingmatrix')
                 )
               ),
               h5('Graph-level descriptive indices',
                  icon('angle-double-left'), id="graphleveltitle"),
               wellPanel(id="graphlevelbox",
                 fluidRow(
                  column(4, offset=7,tags$u('Measure')),
                 fluidRow(
                  column(4, p('Density:', class='stitle')),
                  column(3, p(textOutput('gden'), class='snum'))),

                 fluidRow(
                  column(4, p('Degree:', class='stitle')),
                  column(3, p(textOutput('gdeg'), class='snum')),
                  column(4, selectInput('gdegcmode', label=NULL,
                                        choices=c('indegree', 'outdegree', 'total')
                                        ))),
                 fluidRow(
                  column(4, p('Reciprocity:', class='stitle')),
                  column(3, p(textOutput('grecip'), class='snum')),
                  column(4, selectInput('grecipmeas',label=NULL,
                             choices=c('dyadic','dyadic.nonnull','edgewise',
                                       'edgewise.lrr','correlation')))),
                 fluidRow(
                  column(4, p('Transitivity:'), class='stitle'),
                  column(3, p(textOutput('gtrans'), class='snum')),
                  column(4, selectInput('gtransmeas',label=NULL,
                                choices=c('weak','strong','weakcensus',
                                          'strongcensus','rank','correlation'))
                         )),
                 fluidRow(
                  column(4, p('Betweenness:', class='stitle')),
                  column(3, p(textOutput('gbetw'), class='snum')),
                  column(4, selectInput('gbetwcmode', label=NULL,
                                         choices=c('directed','undirected',
                                                   'endpoints','proximalsrc',
                                                   'proximaltar','proximalsum',
                                                   'lengthscaled', 'linearscaled'))
                         )),
                 fluidRow(
                  column(4, p('Closeness:', class='stitle')),
                  column(3, p(textOutput('gclose'), class='snum')),
                  column(4, selectInput('gclosecmode', label=NULL,
                                        choices=c('directed','undirected',
                                                  'suminvdir','suminvundir')))),
                 fluidRow(
                   column(4, p('Stress Centrality:', class='stitle')),
                   column(3, p(textOutput('gstress'), class='snum')),
                   column(4, selectInput('gstresscmode', label=NULL,
                                         choices=c('directed','undirected')))
                 ),
                 fluidRow(
                   column(4, p('(Harary) Graph Centrality:', class='stitle')),
                   column(3, p(textOutput('ggraphcent'), class='snum')),
                   column(4, selectInput('ggraphcentcmode', label=NULL,
                                         choices=c('directed', 'undirected')))
                 ),
                 fluidRow(
                   column(4, p('Eigenvector Centrality:', class='stitle')),
                   column(3, p(textOutput('gevcent'), class='snum')),
                   column(4, br())
                 ),
                 fluidRow(
                   column(4, p('Information Centrality:', class='stitle')),
                   column(3, p(textOutput('ginfocent'), class='snum')),
                   column(4, selectInput('ginfocentcmode',label=NULL,
                                         choices=c('weak', 'strong', 'upper',
                                                   'lower')))
                 )


               )),

               h5('Vertex-level descriptive indices',
                  icon('angle-double-left'), id="nodeleveltitle"),
               wellPanel(id="nodelevelbox",
                 fluidRow(
                     column(2,span("Vertex index:")),
                     column(5, numericInput('nodeind', label=NULL, value=1,
                                            min=1))
                     ),
                 tags$hr(),
                 fluidRow(
                   column(2, offset=3, tags$u('Current vertex')),
                   column(3, tags$u('Centrality mode')),
                   column(2, tags$u('Min')),
                   column(2, tags$u('Max'))),
                   fluidRow(
                     column(3, p('Degree:', class='stitle')),
                     column(2, p(textOutput('ndeg'), class='snum')),
                     column(3, selectInput('ndegcmode', label=NULL,
                                           choices=c('indegree', 'outdegree', 'total')),
                            class = "smallselect"),
                     column(2, p(textOutput('ndegmin'), class='snum', align='center')),
                     column(2, p(textOutput('ndegmax'), class='snum', align='center'))
                     ),
                   fluidRow(
                     column(3, p('Betweenness:', class='stitle')),
                     column(2, p(textOutput('nbetw'), class='snum')),
                     column(3, selectInput('nbetwcmode', label=NULL,
                                           choices=c('directed','undirected',
                                                     'endpoints','proximalsrc',
                                                     'proximaltar','proximalsum',
                                                     'lengthscaled', 'linearscaled')),
                            class = "smallselect"),
                     column(2, p(textOutput('nbetwmin'), class='snum')),
                     column(2, p(textOutput('nbetwmax'), class='snum'))
                     ),
                   fluidRow(
                     column(3, p('Closeness:', class='stitle')),
                     column(2, p(textOutput('nclose'), class='snum')),
                     column(3, selectInput('nclosecmode', label=NULL,
                                           choices=c('directed','undirected',
                                                     'suminvdir','suminvundir')),
                            class = "smallselect"),
                     column(2, p(textOutput('nclosemin'))),
                     column(2, p(textOutput('nclosemax')))
                     ),
                   fluidRow(
                     column(3, p('Stress Centrality:', class='stitle')),
                     column(2, p(textOutput('nstress'), class='snum')),
                     column(3, selectInput('nstresscmode', label=NULL,
                                           choices=c('directed','undirected')),
                            class = "smallselect"),
                     column(2, p(textOutput('nstressmin'))),
                     column(2, p(textOutput('nstressmax')))
                     ),
                   fluidRow(
                     column(3, p('(Harary) Graph Centrality:', class='stitle')),
                     column(2, p(textOutput('ngraphcent'), class='snum')),
                     column(3, selectInput('ngraphcentcmode', label=NULL,
                                           choices=c('directed', 'undirected')),
                            class = "smallselect"),
                     column(2, p(textOutput('ngraphcentmin'))),
                     column(2, p(textOutput('ngraphcentmax')))
                     ),
                   fluidRow(
                     column(3, p('Eigenvector Centrality:', class='stitle')),
                     column(2, p(textOutput('nevcent'), class='snum')),
                     column(3, br()),
                     column(2, p(textOutput('nevcentmin'))),
                     column(2, p(textOutput('nevcentmax')))
                     ),
                   fluidRow(
                     column(3, p('Information Centrality:', class='stitle')),
                     column(2, p(textOutput('ninfocent'), class='snum')),
                     column(3, selectInput('ninfocentcmode',label=NULL,
                                           choices=c('weak', 'strong', 'upper',
                                                     'lower')),
                            class = "smallselect"),
                     column(2, p(textOutput('ninfocentmin'))),
                     column(2, p(textOutput('ninfocentmax')))
                     )
                 ),
              conditionalPanel(condition = "output.errstate == '1'",
                               div(class = "error", uiOutput("errbox")))

)



      ),br(),br()
),
 column(4,
     tabsetPanel(id='displaytabs',
       tabPanel(title='Display Options', br(),
          wellPanel(
                conditionalPanel(condition='input.plottabs == "Network Plot"',
                   checkboxInput('iso',
                                 label = 'Display isolates',
                                 value = TRUE),
                   checkboxInput('vnames',
                                 label = 'Display vertex names',
                                 value = FALSE),
                   br(),
                   sliderInput('transp',
                               label = 'Vertex opacity',
                               min = 0, max = 1, value = 1),
                   br(),
                   uiOutput("dynamiccolor"),
                   conditionalPanel(condition="Number(output.attrlevels) > 9",
                     column(10,
                            p(id = "closewarning1", icon(name = "remove"), class = "warning"),
                            div(class = "warning", id = "colorwarning1",
                                span(tags$u("Note:"),
                                     "Color palette becomes a gradient for attributes with more than nine levels.")
                            )
                     )),
                   #span(bsAlert(inputId = 'colorwarning'), style='font-size: 0.82em;'),
                   uiOutput('dynamicsize'),
                   br(),
                   actionButton("refreshplot", icon = icon("refresh"),
                                label = "Refresh Plot", class = "btn-sm"),
                   downloadButton('nwplotdownload',
                                  label = "Download Plot", class = "btn-sm")),
                conditionalPanel(condition='input.plottabs == "Attributes"',
                                 selectInput("attrview", label = "View attributes in:",
                                             choices = c("Large table",
                                                         "Small tables",
                                                         "Plot summaries")),
                                 br(),
                                 uiOutput("attrcheck")
                ),
                conditionalPanel(condition='input.plottabs == "Degree Distribution"',
                   uiOutput("dynamiccmode_dd"),
                   uiOutput("dynamiccolor_dd"),
                   tags$label("Y-axis units:"), br(),
                   actionButton("countButton_dd", label="Count of vertices", class="btn-sm active"),
                   actionButton("percButton_dd", label="Percent of vertices", class="btn-sm"),
                   br(), br(),
                   tags$label('Expected values of null models:'), br(),
                   fluidRow(
                     column(10,
                            checkboxInput('uniformoverlay_dd',
                                   label='Conditional uniform graphs (CUG)',
                                   value=FALSE)
                            ),

                      span(icon('question-circle'), id="cughelper_dd", class="helper",
                           div(id="cughelperbox_dd", class="mischelperbox",
                               "Draws from the distribution of simple random graphs with the same",
                               "fixed density as the observed network. The mean and 95% confidence",
                               "intervals for each degree are plotted."))),
                   fluidRow(
                     column(10,
                            checkboxInput('bernoullioverlay_dd',
                                   label='Bernoulli random graphs (BRG)',
                                   value=FALSE)
                            ),
                     span(icon('question-circle'), id="brghelper_dd", class="helper",
                          div(id="brghelperbox_dd", class="mischelperbox",
                              "Draws from the distribution of simple random graphs with the same",
                              "stochastic tie probability as the observed network.",
                              "The mean and 95% confidence intervals for each degree are plotted."))),
                   br(),
                   downloadButton('degreedistdownload', label = "Download Plot", class="btn-sm")
                  ),
                conditionalPanel(condition='input.plottabs == "Geodesic Distribution"',
                                 tags$label("Y-axis units:"), br(),
                                 actionButton("countButton_gd", "Count of vertex pairs", class="btn-sm active"),
                                 actionButton("percButton_gd", "Percent of vertex pairs", class="btn-sm"),
                                 br(), br(),
                                 tags$label('Expected values of null models:'), br(),
                                 fluidRow(
                                   column(10,
                                     checkboxInput('uniformoverlay_gd',
                                                 label='Conditional uniform graphs (CUG)',
                                                 value=FALSE)
                                          ),
                                   span(icon('question-circle'), id="cughelper_gd", class="helper",
                                        div(id="cughelperbox_gd", class="mischelperbox",
                                            "Draws from the distribution of simple random graphs with the same",
                                            "fixed density as the observed network. The mean and 95% confidence",
                                            "intervals for each degree are plotted."))),
                                 fluidRow(
                                   column(10,
                                     checkboxInput('bernoullioverlay_gd',
                                               label='Bernoulli random graphs (BRG)',
                                               value=FALSE)
                                          ),
                                   span(icon('question-circle'), id="brghelper_gd", class="helper",
                                        div(id="brghelperbox_gd", class="mischelperbox",
                                            "Draws from the distribution of simple random graphs with the same",
                                            "stochastic tie probability as the observed network.",
                                            "The mean and 95% confidence intervals for each degree are plotted."))),
                                 br(),
                                 verbatimTextOutput('infsummary'),
                                 fluidRow(
                                   column(10,
                                          checkboxInput('excludeInfs',
                                                        label=span('Exclude "inf"s from plot'),
                                                        value=FALSE)),
                                   span(icon('question-circle'), id="infhelper_gd", class="helper",
                                        div(id="infhelperbox_gd", class="mischelperbox",
                                            "A pair of nodes without any path connecting",
                                            'it has a geodesic distance of "inf".'))),
                                 br(),
                                 downloadButton('geodistdownload', label= 'Download Plot', class="btn-sm")
                  ),
                conditionalPanel(condition='input.plottabs == "More"',
                                 p("No display options at this time,",
                                   "stay tuned for updates!")
                                 )
                )),
       tabPanel(title='Network Summary', br(),
        verbatimTextOutput('attr2'))
        )
     )
 ),
div(id='plottabhelp', class='helper-btn', icon('question-circle', 'fa-2x')),
 div(class="helper-box", style="display:none",
     p('Use the network plots to gain insight to the observed network.',
       'Edit the display options in the panel on the right and download a PDF of any of the plots.')),
actionLink('plotleft', icon=icon('arrow-left', class='fa-2x'), label=NULL),
actionLink('plotright', icon=icon('arrow-right', class='fa-2x'), label=NULL)
),



# Fit Model ---------------------------------------------------------------

# The output objects for the current dataset and current formula have default values
# specified in server.R to prevent errors from NULL values and so that there are
# helpful messages for the user before they begin entering data.

      tabPanel(title='Fit Model',value='tab4',

          #include progress bar when this tab is loading
           div(class = "busy",
               p("Calculation in progress..."),
               img(src="ajax-loader.gif")
           ),

          fluidRow(
            column(2,
               p('Network:', class="nwlabel"),
               verbatimTextOutput('currentdataset1')
              ),

            column(4,
                   p("ERGM terms:"),
                   div(textInput(inputId="terms", label=NULL, value="edges"),
                       title=paste("Type in term(s) and their arguments.",
                                   "For multiple terms, separate with '+'. ")
                   ),
                   actionButton('addtermButton', 'Add Term(s)', class="btn-primary btn-sm"),
                   actionButton('resetformulaButton', 'Reset Formula', class="btn-sm")


            ),
            column(5,
               tabsetPanel(
                 tabPanel("Term Documentation",
                  br(),
                  div(class="placeholder",
                      fluidRow(
                        column(12,
                               a("Commonly used ergm terms",
                                 href = "http://statnet.github.io/nme/ergmterms.html",
                                 target = "_blank"), br(),
                               a("Term cross-reference tables",
                                 href = "http://cran.r-project.org/web/packages/ergm/vignettes/ergm-term-crossRef.html",
                                 target = "_blank"), br(), br()
                               ),
                        column(6,
                               actionButton("matchingButton", "Compatible terms",
                                            class="btn-sm active"),
                               actionButton("allButton", "All terms",
                                            class="btn-sm")
                               ),
                        column(4, uiOutput("listofterms"))
                      ),
                      fluidRow(
                        column(12,
                               div(id="termdocbox",
                                    uiOutput("termdoc")
                                ),
                                div(id = "termexpand",
                                    icon(name = "angle-double-up"))
                               )

                      )

                  )
                 ),
                 tabPanel("Control Options",
                    div(class = "placeholder",
                    fluidRow(class = "shiftright",
                      column(3, style = "padding-left: 0;",
                        inlineSelectInput('controltype',label = NULL,
                                          choices = c("MCMC"),
                                          style="margin:10px 0px;")),
                      column(5,
                        checkboxInput('controldefault', 'Use default options', value = TRUE))
                    ),
                        conditionalPanel(condition = "input.controltype == 'MCMC'",
                                         class = "shiftright gray", id = "mcmcopt1",
                          fluidRow(
                            column(4,
                                   span("Interval:"),
                                   customNumericInput('MCMCinterval', label = NULL, value = 1024,
                                                      class = "mcmcopt input-mini round"),
                                   title = paste("Number of proposals between sampled statistics.")
                                   ),

                            column(4,
                                   span("Burn-in:"),
                                   customNumericInput('MCMCburnin', label = NULL, value = 16384,
                                                      class = "mcmcopt input-mini round"),
                                   title = paste("Number of proposals before any MCMC sampling is done.",
                                                 "Defaults to 16 times the MCMC interval, unless burn-in is specified after the interval.")
                                   ),

                            column(4,
                                   span("Sample size:"),
                                   customNumericInput('MCMCsamplesize', label = NULL, value = 1024,
                                                      class = "mcmcopt input-mini round"),
                                   title = paste("Number of network statistics, randomly drawn from a given distribution",
                                                 "on the set of all networks, returned by the Metropolis-Hastings algorithm.")
                                   )
                          ),

                          fluidRow(
                              div(span("Other controls:", class = "shiftright"),
                                  customTextInput("customMCMCcontrol", label = NULL, value = "",
                                                  class = "input-small round"),
                                  title = paste("Other arguments to be passed to",
                                       "control.ergm, e.g. MCMC.burnin.retries = 1")
                                  )
                            )),
                        conditionalPanel(condition = "input.controltype == 'MCMLE'",
                                         p("Coming soon"))
                        )))
                     )
            ),
          br(),
         fluidRow(
           column(2,
                  p('Current ergm formula:')),
           column(10,
                  verbatimTextOutput('checkterms_fit'))),
        fluidRow(
           column(2,
                  p('Summary statistics:')),
           column(10,
                  verbatimTextOutput('prefitsum'))),
         fluidRow(column(12,
                         actionButton("fitButton", "Fit Model", class="btn-primary btn-sm"),
                         uiOutput("savemodel"),
                         actionButton("clearmodelButton", label="Clear All Models", class="btn-sm")
         )),
         br(),
         tabsetPanel(id = 'fittingTabs',
           tabPanel('Current Model Summary', br(),
                    verbatimTextOutput('modelfitsum'),
                    downloadButton("modelfitdownload", "Download Summary (.txt)", class="btn-sm")),
           tabPanel('Current Model Fit Report', br(),
                    verbatimTextOutput('modelfit')),
           tabPanel('Model Comparison', br(),
                    verbatimTextOutput('modelcomparison'),
                    downloadButton("modelcompdownload", "Download Comparison (.txt)", class="btn-sm"))
          ), br(),br(),
  div(id='fittabhelp', class='helper-btn', icon('question-circle', 'fa-2x')),
    div(class="helper-box", style="display:none",
      p('Create an ergm formula by typing terms into the text box.',
        'Notice the summary statistics populate for each term added to the formula. ',
        'After fitting the model, the "Fitting" tab will show MCMC iterations (if any) and MLE coefficients,',
        'while the "Summary" tab shows a comprehensive summary of the model fit.',br(),
        'Find more help in the', a('ergm tutorial.',
                        href='http://statnet.csde.washington.edu/workshops/SUNBELT/EUSN/ergm/ergm_tutorial.html',
                        target="_blank"))),
actionLink('fitleft', icon=icon('arrow-left', class='fa-2x'), label=NULL),
actionLink('fitright', icon=icon('arrow-right', class='fa-2x'), label=NULL)
          ),

# MCMC Diagnostics --------------------------------------------------------


tabPanel(title='MCMC Diagnostics', value='tab5',
         #include progress bar when this tab is loading
         div(class = "busy",
             p("Calculation in progress..."),
             img(src="ajax-loader.gif")
         ),

         fluidRow(
           column(2,
                  p('Network:', class="nwlabel"),
                  verbatimTextOutput('currentdataset_mcmc')),
           column(10,
                  div(p('ergm formula:',style="display:inline;"),
                    uiOutput('uichoosemodel_mcmc'), class="nwlabel"),
                  verbatimTextOutput('checkterms_mcmc'))
         ),
         br(),
         tags$hr(),
         tabsetPanel(id='mcmctabs',
           tabPanel('Plot', br(),
                    wellPanel(
                      class = "mcmcwarning",
                      p("Recent changes in the ergm estimation algorithm mean that these plots",
                        "can no longer be used to ensure that the mean statistics from the model match the",
                        "observed network statistics. For that functionality, please use the GOF page.")
                      #p(icon("close"), class = "topright")
                       ),
                    #intercept error and give friendly message when MCMC doesn't run
                    conditionalPanel(condition="output.diagnostics == 'MCMC was not run or MCMC sample was not stored.'",
                                column(1,span(class='helper', id='mcmchelper', icon('question-circle')),
                                       style='width:20px; margin-left:0px'),
                                column(11,pre('MCMC was not run or MCMC sample was not stored.'),
                                       style='margin-left:0px;'),
                                column(3,
                                       div(class='mischelperbox', id='mcmchelpbox',
                                        "MCMC is only run when at least one of the terms in the model represents",
                                        "dyad dependence (e.g., degree terms, or triad related terms).  For",
                                        "models with only dyadic independent terms, estimation relies on",
                                        "traditional maximum likelihood algorithms used for generalized linear",
                                        "models."))
                                     ),
                    uiOutput('diagnosticsplotspace'),
                    downloadButton('mcmcplotdownload',label = 'Download Plots', class="btn-sm")),
           tabPanel('Summary', br(),
                    wellPanel(
                      class = "mcmcwarning",
                      p("Recent changes in the ergm estimation algorithm mean that these plots",
                        "can no longer be used to ensure that the mean statistics from the model match the",
                        "observed network statistics. For that functionality, please use the GOF page.")
                      #p(icon("close"), class = "topright")
                    ),
                    verbatimTextOutput('diagnostics'))
         ),
         div(id='mcmctabhelp', class='helper-btn', icon('question-circle', 'fa-2x')),
         div(class="helper-box", style="display:none",
             p('Check for model degeneracy. When a model converges properly',
               'the MCMC sample statistics should vary randomly around the',
               'observed values at each step, and the difference between the',
               'observed and simulated values of the sample statistics should',
               'have a roughly bell shaped distribution, centered at 0.')),
         actionLink('mcmcleft', icon=icon('arrow-left', class='fa-2x'), label=NULL),
         actionLink('mcmcright', icon=icon('arrow-right', class='fa-2x'), label=NULL)

),

# Goodness of Fit ---------------------------------------------------------


tabPanel(title='Goodness of Fit',value='tab6',

         #include progress bar when this tab is loading
         div(class = "busy",
             p("Calculation in progress..."),
             img(src="ajax-loader.gif")
         ),

         fluidRow(
           column(2,
                  p('Network:', class="nwlabel"),
                  verbatimTextOutput('currentdataset_gof')),
           column(10,
                  div(p('ergm formula:',style="display:inline;"),
                  uiOutput('uichoosemodel_gof'), class="nwlabel"),
                  verbatimTextOutput('checkterms_gof'))
          ),
         p('If you do not specify a term the default formula for undirected
           networks is ', code('~ degree + espartners + distance'), 'and for
           directed networks is ', code('~ idegree + odegree + espartners +
                                        distance'), '.'),
         fluidRow(
            column(2,
                   p("Goodness of fit term:"),
                   selectInput('gofterm', label = NULL,
                               c('Default', 'degree','idegree','odegree',
                                 'distance', 'espartners','dspartners', 'triadcensus',
                                 'model')
                                 )),
            column(1, actionButton('gofButton', 'Run', class="shiftdown"))
            ),
         br(),
     tabsetPanel(
       tabPanel("Current Model", br(),
                fluidRow(
                  column(5,
                         verbatimTextOutput('gofsummary')),
                  column(7,
                         uiOutput('gofplotspace'),
                         downloadButton('gofplotdownload', label = 'Download Plots', class="btn-sm")))
                ),
       tabPanel("Compare Saved Models",align="center", br(),
                uiOutput('gofplotcompspace'),
                fluidRow(align="left",
                         downloadButton('gofplotcompdownload',
                                        label='Download Plots', class="btn-sm"),
                         br())
                )
       ),

   div(id='goftabhelp', class='helper-btn', icon('question-circle', 'fa-2x')),
   div(class="helper-box", style="display:none",
       p('Test how well your model fits the original data by choosing a network',
          'statistic that is not in the model, and comparing the value of this',
          'statistic observed in the original network to the distribution of values',
          'you get in simulated networks from your model.')),
   actionLink('gofleft', icon=icon('arrow-left', class='fa-2x'), label=NULL),
   actionLink('gofright', icon=icon('arrow-right', class='fa-2x'), label=NULL)
),

# Simulations -------------------------------------------------------------


tabPanel(title='Simulations', value='tab7',
         fluidRow(
           column(7,
              fluidRow(
                column(4,
                    p('Network:', verbatimTextOutput('currentdataset_sim'))),
                column(4,
                       p("Number of simulations:"),
                       numericInput('nsims', label = NULL,
                                    min = 1, value = 1)
                       ),
                column(1,
                       actionButton('simButton', 'Simulate', class = "shiftdown")
                       )
                ),
              div(p('ergm formula:',style="display:inline;"),
              uiOutput('uichoosemodel_sim'), class="nwlabel"),
              verbatimTextOutput('checkterms_sim')
              ),
           column(5,
                tabsetPanel(
                  tabPanel("Control Options",
                           fluidRow(
                             column(3, class = "shiftright",
                                    inlineSelectInput('simcontroltype',label=NULL,
                                                      choices=c("MCMC","Parallel"),
                                                      style="margin:10px 0px;")),
                             column(7,
                                    checkboxInput('simcontroldefault','Use default options', value=TRUE))
                           ),
                           conditionalPanel(condition="input.simcontroltype == 'MCMC'",
                                            class="shiftright gray", id = "mcmcopt2",
                             fluidRow(
                                    column(5,
                                        span("Interval:"),
                                        customNumericInput('simMCMCinterval',label=NULL, value=1024,
                                                           class="mcmcopt input-mini round"),
                                        title=paste("Number of proposals between sampled statistics.")
                                        ),
                                    column(5,
                                        span("Burn-in:"),
                                        customNumericInput('simMCMCburnin', label=NULL, value=16384,
                                                           class="mcmcopt input-mini round"),
                                        title=paste("Number of proposals before any MCMC sampling is done.",
                                                    "Defaults to 16 times the MCMC interval, unless burn-in is specified after the interval.")
                                        )

                             ),
                             fluidRow(
                                    div(
                                        span("Other controls:", class = "shiftright"),
                                        customTextInput("simcustomMCMCcontrol", label=NULL, value="",
                                                        class = "input-mini round"),
                                        title=paste("Type in other arguments to be passed to control.simulate,",
                                                    "e.g. MCMC.init.maxedges=200")
                                        )
                             )),
                           conditionalPanel(condition="input.simcontroltype == 'Parallel'",
                                           p("Coming soon"))
                           )
                  )
           )
         ),

         tags$hr(),

         fluidRow(
           column(7,
             tabsetPanel(id="simplotpanel",
             tabPanel("Network Plots", br(),
                 column(5,
                        numericInput('thissim',
                                     label = 'Choose a simulation to plot:',
                                     min = 1, value = 1)
                        ),
                 plotOutput('simplot')
                      ),
             tabPanel("Simulation Statistics",
                      conditionalPanel("output.simnum >1",
                          div(plotOutput('simstatsplot'),
                              title=paste("Statistics from each simulation,",
                                          "plotted over horizontal \n",
                                          "lines of the corresponding target statistics."))
                       )
               )
              ), br(), br()
             ),


         column(4,
                  tabsetPanel(
                    tabPanel('Display Options', br(),
                        conditionalPanel("input.simplotpanel == 'Network Plots'",
                             wellPanel(
                               checkboxInput('iso2',
                                             label = 'Display isolates?',
                                             value = TRUE),
                               checkboxInput('vnames2',
                                             label = 'Display vertex names?',
                                             value = FALSE),
                               br(),
                               sliderInput('transp2',
                                           label = 'Vertex opacity',
                                           min = 0, max = 1, value = 1),
                               br(),
                               uiOutput('dynamiccolor2'),
#                                          span(bsAlert(inputId = 'colorwarning2'), style='font-size: 0.82em;'),
                               uiOutput('dynamicsize2'),
                               downloadButton('simplotdownload',
                                              label = 'Download Plot', class="btn-sm"))
                        ),
                        conditionalPanel("input.simplotpanel == 'Simulation Statistics'",
                               conditionalPanel("output.simnum > 1",
                                         plotOutput('simstatslegend'),
                                         downloadButton('simstatsplotdownload',
                                                        label='Download Plot', class="btn-sm")
                                         ))
                      ),
                    tabPanel('Simulation Summary', br(),
                         wellPanel(
                         conditionalPanel(condition="output.simnum != 1",
                                verbatimTextOutput('simsummary'),
                                verbatimTextOutput('simcoef'),
                                verbatimTextOutput('simstatslabel'),
                                conditionalPanel("output.simnum < 10",
                                  verbatimTextOutput('simstats')),
                                conditionalPanel("output.simnum >= 10",
                                  verbatimTextOutput('simstats2'))
                                            ),
                         conditionalPanel(condition="output.simnum == 1",
                           verbatimTextOutput('simsummary2')
                           ),
                         br(),
                         fluidRow(
                           column(7,
                                downloadButton('simstatsdownload',
                                        label = 'Download Statistics', class="btn-sm")),
                           column(4,
                              div(title=paste0(".txt: Summary of simulations",
                                              " plus full list of statistics. \n",
                                              ".csv: Full list of statistics only."),
                                radioButtons('simstatsfiletype', label=NULL,
                                             choices=c('.txt','.csv'))
                                )
                              )
                           )


                           )
                        )
                    )
           )),
         div(id='simtabhelp', class='helper-btn', icon('question-circle', 'fa-2x')),
         div(class="helper-box", style="display:none",
             p('Choose how many simulations to run and click "Simulate".',
               'Plot any individual simulation, or compare',
               'simulation statistics with target statistics.',
               'Download any of the plots or a .csv file of the',
               'simulation statistics.')),
         actionLink('simleft', icon=icon('arrow-left', class='fa-2x'), label=NULL),
         actionLink('simright', icon=icon('arrow-right', class='fa-2x'), label=NULL)
         ),

# Help --------------------------------------------------------------------


tabPanel(title='Help', value='tab8',
         sidebarLayout(position = 'right',
                       sidebarPanel(
                         h5(tags$u('Resources')),
                         div(title = "Wiki page for statnetWeb",
                             a("About statnetWeb",
                               href = "https://statnet.csde.washington.edu/trac/wiki/statnetWeb",
                               target = "_blank")),
                         div(title=paste("Homepage of the statnet project with tutorials,",
                                         "publications and recent news."),
                             a("About statnet",
                               href = "https://statnet.csde.washington.edu/trac", target = "_blank")
                         ),

                         column(11, offset = 1,
                                span(id="linktitle1",'Key background papers',icon('angle-double-left')),br(),
                                div(id="linkbox1",
                                    a("ergm: Journal of Statistical Software",
                                      href = "http://www.jstatsoft.org/v24/i03/", target = "_blank"),
                                    br(),
                                    a("Using ergm: Journal of Statistical Software",
                                      href = "http://www.jstatsoft.org/v24/i04/", target = "_blank")),

                                span(id="linktitle2",'Tutorials and documentation',icon('angle-double-left')),br(),
                                div(id="linkbox2",
                                    a("ergm tutorial from Sunbelt EUSN 2014 Workshop",
                                      href = "http://statnet.csde.washington.edu/workshops/SUNBELT/EUSN/ergm/ergm_tutorial.html",
                                      target= "_blank"),
                                    br(),
                                    a("ergm documentation on CRAN",
                                      href = "http://cran.r-project.org/web/packages/ergm/ergm.pdf",
                                      target = "_blank")),
                                style="margin-bottom:10px;"),
                         br(),
                         div(a("statnetWeb on GitHub", href="https://github.com/statnet/statnetWeb",
                               target="_blank")),
                         div(a("Shiny: a web application framework for R", href="http://shiny.rstudio.com/",
                               target="_blank"))
                       ),
                       mainPanel(
                         h5(tags$u('Help with statnetWeb')),
                         p("This app is maintained on GitHub. To request new features or report a bug,",
                           "please interact with the",
                           a("repository", href='https://github.com/statnet/statnetWeb',
                             target="_blank"),
                           "or email the statnet_help mailing list (below)."),
                         h5(tags$u('Help with statnet software')),
                         p("The best way to contact us with questions, comments or suggestions",
                           "is through the statnet users group mailing list."),
                         p("To post and receive messages from this list, you need to join.",
                           "See the",
                           a("statnet_help info page",
                             href = "https://mailman.u.washington.edu/mailman/listinfo/statnet_help",
                             target = "_blank"),
                           "for instructions on how to join."),
                         p("You can use the mailing list to:"),
                         tags$ul(
                           tags$li("get help from the statnet development team (and other users)"),
                           tags$li("post questions, comments and ideas to other users"),
                           tags$li("be informed about statnet updates"),
                           tags$li("learn about bugs (and bug fixes)")
                         ),
                         p("Once you have joined the list, you can post your questions and comments to",
                           strong("statnet_help@u.washington.edu")),
                         p("A full list of all messages posted to this list is available",
                           a("here.",
                             href = "https://mailman.u.washington.edu/mailman/private/statnet_help",
                             target = "_blank"))
                         ))
         )


  )
)
