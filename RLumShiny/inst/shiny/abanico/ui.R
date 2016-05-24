## UI.R
library(RLumShiny)

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
                                      )
                                    ),
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
                           
                           tabPanel("Login",
                                    div(align = "center", 
                                        fluidRow(
                                          column(width = 2, offset = 5,
                                                 h5("Login")
                                          ),
                                          column(width = 1, offset = 2,
                                                 actionButton("infoButton", "", icon = icon("question")))
                                        )
                                    ),
                                    uiOutput("infobox"),
                                    div(align = "center", htmlOutput("status.msg")),
                                    # Login for own data sets
                                    fluidRow(
                                      column(width = 6,
                                             textInput(inputId = "user", "Username", value = "")
                                      ),
                                      column(width = 6,
                                             textInput(inputId = "pw", "Password", value = "")
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 3,
                                             actionButton("enter", "Enter", class = "btn btn-success")
                                      ),
                                      column(width = 6, offset = 3,
                                             actionButton("newAccount","Create new account")
                                      )
                                    ),
                                    div(align = "center", h5("Data sets")),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "serverData", "Plot data set", choices = "")
                                      ),
                                      column(width = 6,
                                             selectInput(inputId = "del.dataset", "Delete selected data sets", choices = "", multiple = TRUE)
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 6, offset =6,
                                             actionButton("delete", "Delete", class = "btn btn-danger")
                                      )),
                                    div(align = "center", h5("Upload new data set")),
                                    uiOutput("file2server"),
                                    textInput(inputId = "filenameUpload", label = "Filename", "choose a name"),
                                    fluidRow(
                                      column(width = 6,
                                             actionButton("upload", "Upload file", icon = icon("upload"))
                                      ),
                                      column(width = 1, offsett = 5,
                                             #helpPopup("Help", "...", placement='top', trigger='click')
                                             actionButton("logout", "Logout")
                                      )
                                    )
                           ),
                           
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
                                                                                 "Bottom right" = "bottomright"))),
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
                                    tooltip(refId = "statlabels", text = "Additional labels of statistically important values in the plot."),
                                    br(),
                                    div(align = "center", h5("Error bars")),
                                    checkboxInput(inputId = "errorbars",
                                                  label = "Show error bars",
                                                  value = FALSE),
                                    tooltip(refId = "errorbars", text = "Option to show D<sub>e</sub>-errors as error bars on D<sub>e</sub>-points. Useful in combination with hidden y-axis and 2&sigma; bar")
                           ),##EndOf::Tab_2
                           
                           # Tab 3: input that refer to the plot rather than the data
                           tabPanel("Plot", 
                                    div(align = "center", h5("Title")),
                                    fluidRow(
                                      column(width = 6,
                                             textInput(inputId = "main", 
                                                       label = "Title", 
                                                       value = "Abanico Plot")
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
                                             div(id="cent", # DIV with id for tooltip
                                                 uiOutput(outputId = "centValue")
                                             ),
                                             tooltip(refId = "cent", text =  "User-defined central value, primarily used for horizontal centering of the z-axis")
                                      ),
                                      column(width = 6,
                                             # inject sliderInput from Server.R
                                             div(id="bwKDE",
                                                 uiOutput(outputId = "bw")
                                             ),
                                             tooltip(refId = "bwKDE", text = "Bin width of the kernel density estimate")
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 6,
                                             div(id="pratiodiv",
                                                 sliderInput(inputId = "p.ratio", 
                                                             label = "Plot ratio", 
                                                             min=0.25, max=0.90,
                                                             value=0.75, step=0.01, round= FALSE)
                                             ),
                                             tooltip(refId = "pratiodiv", text = "Relative space given to the radial versus the cartesian plot part, default is 0.75.")
                                      ),
                                      column(width = 6,
                                             sliderInput(inputId = "cex", 
                                                         label = "Scaling factor",
                                                         min = 0.5, max = 2, 
                                                         value = 1.0, step = 0.1)
                                      )
                                    ),
                                    br(),
                                    div(align = "center", h5("Centrality & Dispersion")),
                                    checkboxInput(inputId = "addBar", 
                                                  label = "Second 2-Sigma Bar",
                                                  value = FALSE),
                                    # centrality can either be a keyword or numerical input
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "centrality", 
                                                         label = "Centrality",
                                                         list("Mean" = "mean",
                                                              "Median" = "median", 
                                                              "Weighted mean" = "mean.weighted", 
                                                              "Weighted median" = "median.weighted",
                                                              "Custom value" = "custom")),
                                             tooltip(refId = "centrality", text = "Measure of centrality, used for the standardisation, centering the plot and drawing the central line. When a second 2&sigma; bar is plotted the dataset is centered by the median.")
                                      ),
                                      column(width = 6,
                                             selectInput(inputId = "dispersion", 
                                                         label = "Measure of dispersion",
                                                         list("1 sigma" = "sd",
                                                              "2 sigma" = "2sd", 
                                                              "95% quartile" = "qr",
                                                              "Custom quartile" = "custom")),
                                             tooltip(refId = "dispersion", text = "Measure of dispersion, used for drawing the polygon that depicts the spread in the dose distribution.")
                                      )
                                    ),
                                    conditionalPanel(condition = "input.dispersion == 'custom'",
                                                     numericInput(inputId = "cinn",
                                                                  label = "x% quartile",
                                                                  value = 25,
                                                                  min = 0, 
                                                                  max = 100, 
                                                                  step = 1)
                                    ),
                                    fluidRow(
                                      column(width = 6,
                                             conditionalPanel(condition = "input.centrality == 'custom'",
                                                              uiOutput("centralityNumeric"))
                                      ),
                                      
                                      column(width = 6,
                                             conditionalPanel(condition = "input.centrality == 'custom'",
                                                              uiOutput("centralityNumeric2"))
                                      )
                                    ),
                                    div(align = "center", h5("Central line")),
                                    
                                    fluidRow(
                                      column(width = 6,
                                             numericInput(inputId = "lwd", 
                                                          label = "Line width #1", 
                                                          min = 0, max = 5, 
                                                          value = 1)
                                      ),
                                      column(width = 6,
                                             numericInput(inputId = "lwd2", 
                                                          label = "Line width #2", 
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
                                    div(align = "center", h5("Further options")),
                                    fluidRow(
                                      column(width = 6,
                                             checkboxInput(inputId = "rug",
                                                           label = "Add rug",
                                                           value = FALSE),
                                             tooltip(refId = "rug", text = "Option to add a rug to the KDE part, to indicate the location of individual values")
                                      ),
                                      column(width = 6,
                                             checkboxInput(inputId = "rotate",
                                                           label = "Rotate plot",
                                                           value = FALSE),
                                             tooltip(refId = "rotate", text = "Option to rotate the plot by 90&deg;.")
                                      )
                                    ),
                                    checkboxInput(inputId = "kde", label = "KDE", value = TRUE),
                                    tooltip(refId = "kde", text = "Option to add a KDE plot to the dispersion part."),
                                    checkboxInput(inputId = "histogram", label = "Histogram", value = TRUE),
                                    tooltip(refId = "histogram", text = "Option to add a histogram to the dispersion part. Only meaningful when not more than one data set is plotted."),
                                    checkboxInput(inputId = "dots", label = "Dots", value = TRUE),
                                    tooltip(refId = "dots", text = "Option to add a dot plot to the dispersion part. If number of dots exceeds space in the dispersion part, a square indicates this.")
                           ),##EndOf::Tab_3
                           
                           # Tab 4: modify axis parameters
                           tabPanel("Axis",
                                    div(align = "center", h5("X-axis")),
                                    fluidRow(
                                      column(width = 4,
                                             textInput(inputId = "xlab1", 
                                                       label = "Label x-axis (upper)",
                                                       value = "Relative error [%]")
                                      ),
                                      column(width = 4,
                                             textInput(inputId = "xlab2", 
                                                       label = "Label x-axis (lower)",
                                                       value = "Precision")
                                      ),
                                      column(width = 4,
                                             textInput(inputId = "xlab3", 
                                                       label = "Label x-axis (KDE)",
                                                       value = "Density")
                                      )
                                    ),
                                    # inject sliderInput from Server.R
                                    uiOutput(outputId = "xlim"),
                                    br(),
                                    div(align = "center", h5("Y-axis")),
                                    checkboxInput(inputId = "yaxis",
                                                  label = "Show y-axis",
                                                  value = TRUE),
                                    tooltip(refId = "yaxis", text = "Option to hide y-axis labels. Useful for data with small scatter."),
                                    textInput(inputId = "ylab", 
                                              label = "Label y-axis",
                                              value = "Standardised estimate"),
                                    uiOutput("ylim"),
                                    br(),
                                    div(align = "center", h5("Z-axis")),
                                    checkboxInput(inputId = "logz",
                                                  label = "Logarithmic z-axis",
                                                  value = TRUE),
                                    tooltip(refId = "logz", text = "Option to display the z-axis in logarithmic scale."),
                                    textInput(inputId = "zlab", 
                                              label = "Label z-axis",
                                              value = "Equivalent dose [Gy]"),
                                    # inject sliderInput from Server.R
                                    uiOutput(outputId = "zlim")
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
                                                              HTML("Choose a color<br>"),
                                                              jscolorInput(inputId = "jscol1"))
                                      )
                                    ),
                                    br(),
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
                                                              HTML("Choose a color<br>"),
                                                              jscolorInput(inputId = "jscol2"))
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
                                    div(align = "center", h5("Dispersion bar")),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "polygon", 
                                                         label = "Dispersion bar color #1",
                                                         choices = list("Grey" = "grey80",
                                                                        "Custom" = "custom",
                                                                        "None" = "none")),
                                             tooltip(refId = "polygon", attr = "for", text = "Colour of the polygon showing the dose dispersion around the central value.")
                                      ),
                                      column(width = 6,
                                             selectInput(inputId = "polygon2", 
                                                         label = "Dispersion bar color #2",
                                                         choices = list("Grey" = "grey80",
                                                                        "Custom" = "custom",
                                                                        "None" = "none"))
                                      )
                                    ),
                                    fluidRow(
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.polygon == 'custom'",
                                                              jscolorInput(inputId = "rgbPolygon",
                                                                           label = "Choose a color"))
                                      ),
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.polygon2 == 'custom'",
                                                              jscolorInput(inputId = "rgbPolygon2",
                                                                           label = "Choose a color"))
                                      )
                                    ),
                                    sliderInput(inputId = "alpha.polygon", 
                                                label = "Transparency", 
                                                min = 0, max = 100, 
                                                step = 1, value = 66),
                                    br(),
                                    div(align = "center", HTML("<h5>2&sigma; bar</h5>")),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "bar", label = HTML("2&sigma; bar color"),
                                                         choices = list("Grey" = "grey50",
                                                                        "Custom" = "custom",
                                                                        "None" = "none")),
                                             tooltip(refId = "bar", attr = "for", text = "Colour of the bar showing the 2-sigma range of the dose error around the central value.")
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
                                    br(),
                                    div(align = "center", h5("Grid")),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "grid", label = "Grid color",
                                                         selected = "none",
                                                         list("Grey" = "grey90",
                                                              "Custom" = "custom",
                                                              "None" = "none"))
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
                                                step = 1, value = 50)
                           ),##EndOf::Tab_7
                           
                           # Tab 8: add and customize legend
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
                           
                           # Tab 9: Filter data 
                           tabPanel("Filter",
                                    div(align = "center", h5("Primary data set")),
                                    selectInput(inputId = "filter.prim", label = "Choose values to exclude",
                                                choices = "", multiple = TRUE, selected = ""),
                                    div(align = "center", h5("Secondary data set")),
                                    selectInput(inputId = "filter.sec", label = "Choose values to exclude",
                                                choices = "", multiple = TRUE, selected = ""),
                                    actionButton(inputId = "exclude", label = "Exclude")
                           ),##EndOf::Tab_9
                           
                           # Tab 10: Layout
                           tabPanel("Layout",
                                    div(align = "center", h5("Layout")),
                                    div(id = "layout", 
                                        selectInput(inputId = "layout", 
                                                    label = "Choose layout", 
                                                    selected = "default",
                                                    choices = c("Default"="default",
                                                                "Journal"="journal"))
                                    ),
                                    tooltip(refId = "layout", text = "The optional parameter layout allows to modify the entire plot more sophisticated. Each element of the plot can be addressed and its properties can be defined. This includes font type, size and decoration, colours and sizes of all plot items. To infer the definition of a specific layout style cf. get_Layout() or type eg. for the layout type \"journal\" get_Layout(\"journal\"). A layout type can be modified by the user by assigning new values to the list object.")
                           ),
                           
                           # Tab 10: save plot as pdf, wmf or eps
                           tabPanel("Export",
                                    radioButtons(inputId = "fileformat", 
                                                 label = "Fileformat", 
                                                 selected = "pdf",
                                                 choices = c("PDF   (Portable Document Format)" = "pdf",
                                                             "SVG   (Scalable Vector Graphics)" = "svg",
                                                             "EPS   (Encapsulated Postscript)" = "eps")),
                                    textInput(inputId = "filename", 
                                              label = "Filename", 
                                              value = "Abanico Plot"),
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
                                        a(href = "https://github.com/tzerk/RLumShiny/tree/master/inst/shiny/abanico", "See the code at GitHub!", target="_blank")
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
              tabPanel("Primary data set", fluidRow(column(width = 12, dataTableOutput("dataset")))),
              tabPanel("Secondary data set", dataTableOutput("dataset2")),
              tabPanel("Central Age Model", dataTableOutput("CAM")),
              tabPanel("R plot code", verbatimTextOutput("plotCode"))
            )###EndOf::tabsetPanel
  )##EndOf::mainPanel
)##EndOf::shinyUI(pageWithSidebar)