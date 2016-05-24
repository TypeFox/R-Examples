library(RLumShiny)
## UI.R

# load example data
data(ExampleData.DeValues)

# pageWithSidebar contains three panels:
# 1 - headerPanel: basically just for a header
# 2 - sidebarPanel: data input
# 3 - mainPanel: data output
pageWithSidebar(  
  # 1 - title = NULL -> Panel will not be shown
  headerPanel(title = NULL),
  # 2- width = 5 -> refers to twitters bootstrap grid system
  # where the the maximum width is 12 that is to be shared between all
  # elements
  sidebarPanel(width = 5,
               # include a tabs in the input panel for easier navigation
               tabsetPanel(id = "tabs", type = "pill", selected = "Data",
                           # Tab 1: Data input
                           tabPanel("Data",
                                    # informational text
                                    div(align = "center", h5("Data upload")),
                                    # file upload button (data set 1)
                                    fileInput(inputId = "file1", 
                                              label = strong("Primary data set"),
                                              accept="text/plain"),
                                    tooltip(refId = "file1", text = tags$img(src='file_structure.png', width='250px')),
                                    # file upload button (data set 2)
                                    fileInput(inputId = "file2", 
                                              label = strong("Secondary data set"), 
                                              accept="text/plain"),
                                    tooltip(refId = "file2", text = tags$img(src='file_structure.png', width='250px')),
                                    # informational text
                                    div(align = "center", h5("Settings")),
                                    fluidRow(
                                      column(width = 6,
                                             # logical: should NA values be excluded?
                                             checkboxInput(inputId = "naExclude", 
                                                           label = "Exclude NA values",
                                                           value = TRUE),
                                             tooltip(refId = "naExclude", text = "Exclude NA values from the data set prior to any further operations.")
                                      ),
                                      column(width = 6,
                                             # logical: file contains headers?
                                             checkboxInput(inputId = "headers", 
                                                           label = "File contains headers", 
                                                           value = FALSE),
                                             tooltip(refId = "headers", text = tags$img(src='file_containsHeader.png', width='250px'))
                                      )),
                                    # char: columns separated by tab, space, comma
                                    radioButtons("sep", "Separator", selected = "\t", inline = TRUE,
                                                 c("Tab" = "\t",
                                                   "Space" = " ",
                                                   "Comma" = ",",
                                                   "Semicolon" = ";")),
                                    tooltip(refId = "sep", text = tags$img(src='file_sep.png', width='400px'), placement = "auto left"),
                                    hr(),
                                    fluidRow(
                                      column(width = 6,
                                             actionButton(inputId = "refresh", label = "Refresh", icon = icon("refresh")),
                                             tooltip(refId = "refresh", text = "Redraw the plot")
                                      ),
                                      column(width = 6,
                                             actionButton(inputId = "exit", label = "Exit", class = "btn btn-danger")
                                      )
                                    )
                           ),##EndOf::Tab_1
                           
                           # Tab 2: Statistical information
                           tabPanel("Statistics",                             
                                    div(align = "center", h5("Summary")),
                                    fluidRow(
                                      column(width = 6,
                                             checkboxInput(inputId = "summary",
                                                           label = "Show summary",
                                                           value = FALSE),
                                             tooltip(refId = "summary", text = "Adds numerical output to the plot")
                                      ),
                                      column(width = 6,
                                             selectInput(inputId = "sumpos",
                                                         label = "Summary position",
                                                         selected = "topleft",
                                                         choices = list("Subtitle" = "sub",
                                                                        "Center" = "center",
                                                                        Top=c("Top" = "top",
                                                                              "Top left" = "topleft",
                                                                              "Top right"= "topright"),
                                                                        Bottom=c("Bottom" = "bottom",
                                                                                 "Bottom left" = "bottomleft",
                                                                                 "Bottom right" = "bottomright")
                                                         )),
                                             tooltip(refId = "sumpos", attr = "for", text = "Position of the statistical summary. The keyword \"Subtitle\" will only work if no plot subtitle is used.")
                                      )
                                    ),
                                    checkboxGroupInput(inputId = "stats",
                                                       label = "Parameters", 
                                                       selected = c("n","mean"),
                                                       choices = c("n" = "n",
                                                                   "Mean" = "mean",
                                                                   "weighted Mean" = "mean.weighted",
                                                                   "Median" = "median",
                                                                   "weighted Median" = "median.weighted",
                                                                   "rel. Standard deviation" = "sdrel",
                                                                   "abs. Standard deviation" = "sdabs",
                                                                   "rel. Standard error" = "serel",
                                                                   "abs. Standard error" = "seabs",
                                                                   #"25 % Quartile" = "q25", #not implemented yet
                                                                   #"75 % Quartile" = "q75", #not implemented yet
                                                                   "KDEmax"  = "kdemax",
                                                                   "Skewness" = "skewness",
                                                                   "Kurtosis" = "kurtosis",
                                                                   "Confidence interval" = "in.ci")),
                                    tooltip(refId = "stats", text = "Statistical parameters to be shown in the summary"),
                                    br(),
                                    div(align = "center", h5("Datapoint labels")),
                                    div(align = "center", checkboxGroupInput(inputId = "statlabels", inline = TRUE,
                                                                             label = NULL, 
                                                                             choices = c("Min" = "min",
                                                                                         "Max" = "max",
                                                                                         "Median" = "median"))),
                                    tooltip(refId = "statlabels", text = "Additional labels of statistically important values in the plot.")
                           ),##EndOf::Tab_2
                           
                           # Tab 3: input that refer to the plot rather than the data
                           tabPanel("Plot", 
                                    div(align = "center", h5("Title")),
                                    fluidRow(
                                      column(width = 6,
                                             textInput(inputId = "main", 
                                                       label = "Title", 
                                                       value = "Radial Plot")
                                      ),
                                      column(width = 6,
                                             textInput(inputId = "mtext", 
                                                       label = "Subtitle", 
                                                       value = "")
                                      )
                                    ),
                                    div(align = "center", h5("Scaling")),
                                    fluidRow(
                                      column(width = 6,
                                             # inject sliderInput from Server.R
                                             uiOutput(outputId = "centValue"),
                                             tooltip(refId = "centValue", text = "User-defined central value, primarily used for horizontal centering of the z-axis")
                                      ),
                                      column(width = 6,
                                             sliderInput(inputId = "cex", 
                                                         label = "Scaling factor",
                                                         min = 0.5, max = 2, 
                                                         value = 1.0, step = 0.1)
                                      )
                                    ),
                                    selectInput(inputId = "centrality", 
                                                label = "Centrality",
                                                list("Mean" = "mean",
                                                     "Median" = "median", 
                                                     "Weighted mean" = "mean.weighted", 
                                                     "Weighted median" = "median.weighted")),
                                    tooltip(refId = "centrality", attr = "for", text = "Measure of centrality, used for the standardisation, centering the plot and drawing the central line.")
                           ),##EndOf::Tab_3
                           
                           # Tab 4: modify axis parameters
                           tabPanel("Axis",
                                    div(align = "center", h5("X-axis")),
                                    fluidRow(
                                      column(width = 6,
                                             textInput(inputId = "xlab1", 
                                                       label = "Label x-axis (upper)",
                                                       value = "Relative error [%]")
                                      ),
                                      column(width = 6,
                                             textInput(inputId = "xlab2", 
                                                       label = "Label x-axis (lower)",
                                                       value = "Precision")
                                      )
                                    ),
                                    # inject sliderInput from Server.R
                                    uiOutput(outputId = "xlim"),
                                    div(align = "center", h5("Y-axis")),
                                    checkboxInput(inputId = "yticks",
                                                  label = HTML("Show &plusmn;2&sigma; label"),
                                                  value = TRUE),
                                    tooltip(refId = "yticks", text = "Option to hide y-axis labels."),
                                    textInput(inputId = "ylab", 
                                              label = "Label y-axis",
                                              value = "Standardised estimate"),
                                    div(align = "center", h5("Z-axis")),
                                    checkboxInput(inputId = "logz",
                                                  label = "Logarithmic z-axis",
                                                  value = TRUE),
                                    tooltip(refId = "logz", text = "Option to display the z-axis in logarithmic scale."),
                                    textInput(inputId = "zlab", 
                                              label = "Label z-axis",
                                              value = "Equivalent dose [Gy]"),
                                    # inject sliderInput from Server.R
                                    uiOutput(outputId = "zlim"),
                                    sliderInput('curvature', 'Z-axis curvature', 
                                                min=0, max=3,
                                                value=4.5/5.5, step=0.01, round=FALSE),
                                    tooltip(refId = "curvature", attr = "for", text = "User-defined plot area ratio (i.e. curvature of the z-axis). If omitted, the default value (4.5/5.5) is used and modified automatically to optimise the z-axis curvature. The parameter should be decreased when data points are plotted outside the z-axis or when the z-axis gets too elliptic.")
                           ),##EndOf::Tab_4
                           
                           # Tab 5: modify data point representation
                           tabPanel("Datapoints",              
                                    div(align = "center", h5("Primary data set")),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "pch",
                                                         label = "Style",
                                                         selected = "17",
                                                         choices = c("Square"= "1",
                                                                     "Circle"="2",
                                                                     "Triangle point up"="3",
                                                                     "Plus"="4",
                                                                     "Cross"="5",
                                                                     "Diamond"="6",
                                                                     "Triangle point down"="7",
                                                                     "Square cross"="8",
                                                                     "Star"="9",
                                                                     "Diamond plus"="10",
                                                                     "Circle plus"="11",
                                                                     "Triangles up and down"="12",
                                                                     "Square plus"="13",
                                                                     "Circle cross"="14",
                                                                     "Square and Triangle down"="15",
                                                                     "filled Square"="16",
                                                                     "filled Circle"="17",
                                                                     "filled Triangle point up"="18",
                                                                     "filled Diamond"="19",
                                                                     "solid Circle"="20",
                                                                     "Bullet (smaller Circle)"="21",
                                                                     "Custom"="custom"))
                                      ),
                                      column(width = 6,
                                             # show only if custom symbol is desired
                                             conditionalPanel(condition = "input.pch == 'custom'",
                                                              textInput(inputId = "custompch", 
                                                                        label = "Insert character", 
                                                                        value = "?"))
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "color", label = "Datapoint color",
                                                         choices = list("Black" = "black",
                                                                        "Grey" = "grey50",
                                                                        "Red" = "#b22222", 
                                                                        "Green" = "#6E8B3D", 
                                                                        "Blue" = "#428bca",
                                                                        "Custom" = "custom"))
                                      ),
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.color == 'custom'",
                                                              jscolorInput(inputId = "rgb",
                                                                        label = "Choose a color"))
                                      )
                                    ),
                                    div(align = "center", h5("Secondary data set")),
                                    fluidRow(
                                      column(width = 6,
                                             ## DATA SET 2
                                             selectInput(inputId = "pch2",
                                                         label = "Style",
                                                         selected = "17",
                                                         choices = c("Square"= "1",
                                                                     "Circle"="2",
                                                                     "Triangle point up"="3",
                                                                     "Plus"="4",
                                                                     "Cross"="5",
                                                                     "Diamond"="6",
                                                                     "Triangle point down"="7",
                                                                     "Square cross"="8",
                                                                     "Star"="9",
                                                                     "Diamond plus"="10",
                                                                     "Circle plus"="11",
                                                                     "Triangles up and down"="12",
                                                                     "Square plus"="13",
                                                                     "Circle cross"="14",
                                                                     "Square and Triangle down"="15",
                                                                     "filled Square"="16",
                                                                     "filled Circle"="17",
                                                                     "filled Triangle point up"="18",
                                                                     "filled Diamond"="19",
                                                                     "solid Circle"="20",
                                                                     "Bullet (smaller Circle)"="21",
                                                                     "Custom"="custom"))
                                      ),
                                      column(width = 6,
                                             # show only if custom symbol is desired
                                             conditionalPanel(condition = "input.pch2 == 'custom'",
                                                              textInput(inputId = "custompch2", 
                                                                        label = "Insert character", 
                                                                        value = "?"))
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "color2", label = "Datapoint color", 
                                                         selected = "#b22222",
                                                         choices = list("Black" = "black",
                                                                        "Grey" = "grey50",
                                                                        "Red" = "#b22222", 
                                                                        "Green" = "#6E8B3D", 
                                                                        "Blue" = "#428bca",
                                                                        "Custom" = "custom"))
                                      ),
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.color2 == 'custom'",
                                                              jscolorInput(inputId = "rgb2",
                                                                           label = "Choose a color"))
                                      )
                                    )
                           ),##EndOf::Tab_5
                           
                           # Tab 6: add additional lines to the plot
                           tabPanel("Lines",
                                    helpText("Here you can add additional lines."),
                                    # options for custom lines:
                                    # 1 - z-value, 2 - color, 3 - label
                                    # only the options for the first line are shown
                                    numericInput(inputId = "line1", 
                                                 label = strong("Line #1"), 
                                                 value =  NA, min = 0),
                                    tooltip(refId = "line1", text = "Numeric values of the additional lines to be added."),
                                    fluidRow(
                                      column(width = 6, 
                                             HTML("Choose a color<br>"),
                                             jscolorInput(inputId = "colline1")
                                      ),
                                      column(width = 6,                                    
                                             textInput(inputId = "labline1",
                                                       label = "Label",
                                                       value = "")
                                      )
                                    ),
                                    # conditional chain: if valid input (i.e. the z-value is > 0) is provided
                                    # for the previous line, show options for a new line (currently up to eight)
                                    conditionalPanel(condition = "input.line1 > 0",
                                                     numericInput(inputId = "line2", strong("Line #2"), NA, min = 0),
                                                     fluidRow(
                                                       column(width = 6, HTML("Choose a color<br>"),jscolorInput(inputId = "colline2")),
                                                       column(width = 6, textInput("labline2","Label",value = ""))
                                                     )
                                    ),
                                    conditionalPanel(condition = "input.line2 > 0",
                                                     numericInput(inputId = "line3", strong("Line #3"), NA, min = 0),
                                                     fluidRow(
                                                       column(width = 6, HTML("Choose a color<br>"),jscolorInput(inputId = "colline3")),
                                                       column(width = 6, textInput("labline3","Label",value = ""))
                                                     )
                                    ),
                                    
                                    conditionalPanel(condition = "input.line3 > 0",
                                                     numericInput(inputId = "line4", strong("Line #4"), NA, min = 0),
                                                     fluidRow(
                                                       column(width = 6, HTML("Choose a color<br>"),jscolorInput(inputId = "colline4")),
                                                       column(width = 6, textInput("labline4","Label",value = ""))
                                                     )
                                    ),
                                    
                                    conditionalPanel(condition = "input.line4 > 0",
                                                     numericInput(inputId = "line5", strong("Line #5"), NA, min = 0),
                                                     fluidRow(
                                                       column(width = 6, HTML("Choose a color<br>"),jscolorInput(inputId = "colline5")),
                                                       column(width = 6, textInput("labline5","Label",value = ""))
                                                     )
                                    ),
                                    
                                    conditionalPanel(condition = "input.line5 > 0",
                                                     numericInput(inputId = "line6", strong("Line #6"), NA, min = 0),
                                                     fluidRow(
                                                       column(width = 6, HTML("Choose a color<br>"),jscolorInput(inputId = "colline6")),
                                                       column(width = 6, textInput("labline6","Label",value = ""))
                                                     )
                                    ),
                                    
                                    conditionalPanel(condition = "input.line6 > 0",
                                                     numericInput(inputId = "line7", strong("Line #7"), NA, min = 0),
                                                     fluidRow(
                                                       column(width = 6, HTML("Choose a color<br>"),jscolorInput(inputId = "colline7")),
                                                       column(width = 6, textInput("labline7","Label",value = ""))
                                                     )
                                    ),
                                    
                                    conditionalPanel(condition = "input.line7 > 0",
                                                     numericInput(inputId = "line8", strong("Line #8"), NA, min = 0),
                                                     fluidRow(
                                                       column(width = 6, HTML("Choose a color<br>"),jscolorInput(inputId = "colline8")),
                                                       column(width = 6, textInput("labline8","Label",value = ""))
                                                     )
                                    )
                                    
                           ),##EndOf::Tab_6
                           
                           # Tab 7: modify the 2-sigma bar (radial plot), grid (both) and polygon (KDE)
                           tabPanel("Bars & Grid",
                                    div(align = "center", h5("Central line")),
                                    fluidRow(
                                      column(width = 6,
                                             numericInput(inputId = "lwd", 
                                                          label = "Central line width #1", 
                                                          min = 0, max = 5, 
                                                          value = 1)
                                      ),
                                      column(width = 6,
                                             numericInput(inputId = "lwd2", 
                                                          label = "Central line width #2", 
                                                          min = 0, max = 5, 
                                                          value = 1)
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "lty", 
                                                         label = "Line type",
                                                         selected = 2,
                                                         choices = list("Blank" = 0,
                                                                        "Solid" = 1,
                                                                        "Dashed" = 2,
                                                                        "Dotted" = 3,
                                                                        "Dot dash" = 4,
                                                                        "Long dash" = 5,
                                                                        "Two dash" = 6))
                                      ),
                                      column(width = 6,
                                             selectInput(inputId = "lty2", 
                                                         label = "Line type",
                                                         selected = 2,
                                                         choices = list("Blank" = 0,
                                                                        "Solid" = 1,
                                                                        "Dashed" = 2,
                                                                        "Dotted" = 3,
                                                                        "Dot dash" = 4,
                                                                        "Long dash" = 5,
                                                                        "Two dash" = 6))
                                             
                                      )
                                    ),
                                    div(align = "center", HTML("<h5>2&sigma; bar</h5>")),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "bar", label = HTML("2&sigma; bar color"),
                                                         choices = list("Grey" = "grey50",
                                                                        "Custom" = "custom",
                                                                        "None" = "none"))
                                      ),
                                      column(width = 6,
                                             selectInput(inputId = "bar2", label = HTML("2&sigma; bar color #2"),
                                                         choices = list("Grey" = "grey50",
                                                                        "Custom" = "custom",
                                                                        "None" = "none"))
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.bar == 'custom'",
                                                              jscolorInput(inputId = "rgbBar",
                                                                        label = "Choose a color"))
                                      ),
                                      column(width = 6,
                                             
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.bar2 == 'custom'",
                                                              jscolorInput(inputId = "rgbBar2",
                                                                        label = "Choose a color"))
                                      )
                                    ),
                                    sliderInput(inputId = "alpha.bar", 
                                                label = "Transparency",
                                                min = 0, max = 100, 
                                                step = 1, value = 66),
                                    div(align = "center", h5("Grid")),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput("grid", "Grid color",
                                                         list("Grey" = "grey",
                                                              "Custom" = "custom",
                                                              "None" = "none")),
                                             tooltip(refId = "grid", attr = "for", text = "colour of the grid lines (originating at [0,0] and stretching to the z-scale). To disable grid lines, use \"none\".")
                                      ),
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.grid == 'custom'",
                                                              jscolorInput(inputId = "rgbGrid",
                                                                        label = "Choose a color"))
                                      )
                                    ),
                                    sliderInput(inputId = "alpha.grid",
                                                label = "Transparency",
                                                min = 0, max = 100, 
                                                step = 1, value = 100)
                           ),##EndOf::Tab_7
                           
                           tabPanel("Legend",
                                    div(align = "center", h5("Legend")),
                                    fluidRow(
                                      column(width = 6,
                                             checkboxInput(inputId = "showlegend", 
                                                           label = "Show legend", 
                                                           value = FALSE),
                                             tooltip(refId = "showlegend", text = "Legend content to be added to the plot.")
                                      ),
                                      column(width = 6,
                                             selectInput(inputId = "legend.pos",
                                                         label = "Legend position",
                                                         selected = "bottomleft",
                                                         choices = c("Top" = "top",
                                                                     "Top left" = "topleft",
                                                                     "Top right"= "topright",
                                                                     "Center" = "center",
                                                                     "Bottom" = "bottom",
                                                                     "Bottom left" = "bottomleft",
                                                                     "Bottom right" = "bottomright"))
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 6,
                                             textInput(inputId = "legendname", 
                                                       label = "Primary data label", 
                                                       value = "primary data")
                                      ),
                                      column(width = 6,
                                             textInput(inputId = "legendname2", 
                                                       label = "Secondary data label", 
                                                       value = "secondary data")
                                      )
                                    )
                           ),##EndOf::Tab_8
                           
                           # Tab 9: save plot as pdf, wmf or eps
                           tabPanel("Export",
                                    radioButtons(inputId = "fileformat", 
                                                 label = "Fileformat", 
                                                 selected = "pdf",
                                                 choices = c("PDF   (Portable Document Format)" = "pdf",
                                                             "SVG   (Scalable Vector Graphics)" = "svg",
                                                             "EPS   (Encapsulated Postscript)" = "eps")),
                                    textInput(inputId = "filename", 
                                              label = "Filename", 
                                              value = "Radial Plot"),
                                    fluidRow(
                                      column(width = 6,
                                             numericInput(inputId = "imgheight",
                                                          label =  "Image height", 
                                                          value = 7)
                                      ),
                                      column(width = 6,
                                             numericInput(inputId = "imgwidth",
                                                          label = "Image width", 
                                                          value = 7)
                                      )
                                    ),
                                    selectInput(inputId = "fontfamily", 
                                                label = "Font", 
                                                selected = "Helvetica",
                                                choices = c("Helvetica" = "Helvetica",
                                                            "Helvetica Narrow" = "Helvetica Narrow",
                                                            "Times" = "Times",
                                                            "Courier" = "Courier",
                                                            "Bookman" = "Bookman",
                                                            "Palatino" = "Palatino")),
                                    tags$hr(),
                                    downloadButton(outputId = "exportFile", 
                                                   label = "Download plot"),
                                    tags$hr(),
                                    helpText("Additionally, you can download a corresponding .R file that contains",
                                             "a fully functional script to reproduce the plot in your R environment!"),
                                    downloadButton(outputId = "exportScript", 
                                                   label = "Download R script")
                           ),##EndOf::Tab_8
                           
                           # Tab 10: further information
                           tabPanel("About",
                                    hr(),
                                    div(align = "center",
                                        # HTML code to include a .png file in the tab; the image file must be in
                                        # a subfolder called "wwww"
                                        img(src="RL_Logo.png", height = 100, width = 100, alt = "R.Lum"),
                                        p("Links:"),
                                        a(href = "http://www.r-luminescence.de", "R.Luminescence project page", target="_blank"),
                                        br(),
                                        a(href = "https://forum.r-luminescence.de", "Message board", target="_blank"),
                                        br(),
                                        a(href = "http://zerk.canopus.uberspace.de/R.Lum", "Online application", target="_blank"),
                                        br(),hr(),
                                        img(src='GitHub-Mark-32px.png', width='32px', height='32px'),
                                        br(),
                                        a(href = "https://github.com/tzerk/RLumShiny/tree/master/inst/shiny/radialplot", "See the code at GitHub!", target="_blank")
                                    )#/div
                           )##EndOf::Tab_9
               )##EndOf::tabsetPanel
  ),##EndOf::sidebarPanel
  
  # 3 - output panel
  mainPanel(width = 7,
            # insert css code inside <head></head> of the generated HTML file:
            # allow open dropdown menus to reach over the container
            tags$head(tags$style(type="text/css",".tab-content {overflow: visible;}")),
            tags$head(includeCSS("www/style.css")),
            # divide output in separate tabs via tabsetPanel
            tabsetPanel(
              tabPanel("Plot", plotOutput(outputId = "main_plot", height = "500px")),
              tabPanel("Primary data set", dataTableOutput("dataset")),
              tabPanel("Secondary data set", dataTableOutput("dataset2")),
              tabPanel("Central Age Model", dataTableOutput("CAM")),
              tabPanel("R plot code", verbatimTextOutput("plotCode"))
            )###EndOf::tabsetPanel
  )##EndOf::mainPanel
)##EndOf::shinyUI(pageWithSidebar)