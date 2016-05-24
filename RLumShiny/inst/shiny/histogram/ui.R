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
                                    # logical: should NA values be excluded?
                                    checkboxInput(inputId = "naExclude", 
                                                  label = "Exclude NA values",
                                                  value = TRUE),
                                    tooltip(refId = "naExclude", text = "Exclude NA values from the data set prior to any further operations."),
                                    # logical: file contains headers?
                                    checkboxInput(inputId = "headers", 
                                                  label = "File contains headers", 
                                                  value = FALSE),
                                    tooltip(refId = "headers", text = tags$img(src='file_containsHeader.png', width='250px')),
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
                                    checkboxInput(inputId = "summary",
                                                  label = "Show summary",
                                                  value = FALSE),
                                    tooltip(refId = "summary", text = "Adds numerical output to the plot"),
                                    selectInput(inputId = "sumpos",
                                                label = "Summary position",
                                                selected = "topleft",
                                                choices = c("Top" = "top",
                                                            "Top left" = "topleft",
                                                            "Top right"= "topright",
                                                            "Center" = "center",
                                                            "Bottom" = "bottom",
                                                            "Bottom left" = "bottomleft",
                                                            "Bottom right" = "bottomright")),
                                    tooltip(refId = "sumpos", attr = "for", text = "Position of the statistical summary. The keyword \"Subtitle\" will only work if no plot subtitle is used."),
                                    checkboxGroupInput(inputId = "stats", 
                                                       label = "Parameters", 
                                                       selected = c("n","mean"),
                                                       choices = c("n" = "n",
                                                                   "Mean" = "mean",
                                                                   "Median" = "median",
                                                                   "weighted Mean" = "mean.weighted",
                                                                   "KDEmax" = "kdemax",
                                                                   "rel. Standard deviation" = "sdrel",
                                                                   "abs. Standard deviation" = "sdabs", 
                                                                   "rel. Standard error" = "serel",
                                                                   "abs. Standard error" = "seabs",
                                                                   "weighted Median" = "median.weighted",
                                                                   "Skewness" = "skewness",
                                                                   "Kurtosis" = "kurtosis",
                                                                   "Confidence interval" = "in.ci")),
                                    tooltip(refId = "stats", text = "Statistical parameters to be shown in the summary"),
                                    div(align = "center", h5("Error bars")),
                                    checkboxInput(inputId = "errorBars",
                                                  label = "Show standard error points",
                                                  value = TRUE),
                                    tooltip(refId = "errorBars", text = "Plot the standard error points over the histogram.")
                           ),##EndOf::Tab_2
                           
                           # Tab 1: Data input
                           tabPanel("Plot",
                                    div(align = "center", h5("Title")),
                                    fluidRow(
                                      column(width = 6,
                                             textInput(inputId = "main", 
                                                       label = "Title", 
                                                       value = "Histogram")
                                      ),
                                      column(width = 6,
                                             textInput(inputId = "mtext", 
                                                       label = "Subtitle", 
                                                       value = "")
                                      )
                                    ),
                                    div(align = "center", h5("Histogram bars")),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "barsColor", label = "Bar color",
                                                         selected = "grey80",
                                                         choices = list("White" = "white",
                                                                        "Black" = "black",
                                                                        "Grey" = "grey80",
                                                                        "Red" = "#b22222", 
                                                                        "Green" = "#6E8B3D", 
                                                                        "Blue" = "#428bca",
                                                                        "Custom" = "custom"))
                                      ),
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.barsColor == 'custom'",
                                                              jscolorInput(inputId = "barsRgb",
                                                                        label = "Choose a color"))
                                      )
                                    ),
                                    sliderInput(inputId = "alpha.bars", 
                                                label = "Bar transparency", 
                                                min = 0, max = 100, 
                                                step = 1, value = 66),
                                    br(),
                                    div(align = "center", h5("Normal curve")),
                                    checkboxInput(inputId = "norm",
                                                  label = "Add normal curve",
                                                  value = FALSE),
                                    tooltip(refId = "norm", text = "Add a normal curve to the histogram. Mean and standard deviation are calculated from the input data. If the normal curve is added, the y-axis in the histogram will show the probability density"),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "normalColor", label = "Normal curve color",
                                                         selected = "red",
                                                         choices = list("Black" = "black",
                                                                        "Grey" = "grey50",
                                                                        "Red" = "#b22222", 
                                                                        "Green" = "#6E8B3D", 
                                                                        "Blue" = "#428bca",
                                                                        "Custom" = "custom"))
                                      ),
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.normalColor == 'custom'",
                                                              jscolorInput(inputId = "normalRgb",
                                                                        label = "Choose a color"))
                                      )
                                    ),
                                    div(align = "center", h5("Scaling")),
                                    sliderInput(inputId = "cex", 
                                                label = "Scaling factor",
                                                min = 0.5, max = 2, 
                                                value = 1.0, step = 0.1),
                                    br(),
                                    div(align = "center", h5("Rugs")),
                                    checkboxInput(inputId = "rugs",
                                                  label = "Add rugs",
                                                  value = TRUE),
                                    tooltip(refId = "rugs", text = "Option to add a rug to the KDE part, to indicate the location of individual values"),
                                    fluidRow(
                                      column(width = 6,
                                             selectInput(inputId = "rugsColor", label = "Rugs color",
                                                         choices = list("Black" = "black",
                                                                        "Grey" = "grey50",
                                                                        "Red" = "#b22222", 
                                                                        "Green" = "#6E8B3D", 
                                                                        "Blue" = "#428bca",
                                                                        "Custom" = "custom"))
                                      ),
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.rugsColor == 'custom'",
                                                              jscolorInput(inputId = "rugsRgb",
                                                                        label = "Choose a color"))
                                      )
                                    )
                           ),##EndOf::Tab_9
                           
                           # Tab 4: modify axis parameters
                           tabPanel("Axis",
                                    div(align = "center", h5("X-axis")),
                                    textInput(inputId = "xlab", 
                                              label = "Label x-axis",
                                              value = "Equivalent dose [Gy]"),
                                    # inject sliderInput from Server.R
                                    uiOutput(outputId = "xlim"),
                                    div(align = "center", h5("Y-axis")),
                                    fluidRow(
                                      column(width = 6,
                                             textInput(inputId = "ylab1", 
                                                       label = "Label y-axis (left)",
                                                       value = "Counts")
                                      ),
                                      column(width = 6,
                                             textInput(inputId = "ylab2", 
                                                       label = "Label y-axis (right)",
                                                       value = "Error")
                                      )
                                    )
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
                                             selectInput(inputId = "pchColor", label = "Datapoint color",
                                                         choices = list("Black" = "black",
                                                                        "Grey" = "grey50",
                                                                        "Red" = "#b22222", 
                                                                        "Green" = "#6E8B3D", 
                                                                        "Blue" = "#428bca",
                                                                        "Custom" = "custom"))
                                      ), 
                                      column(width = 6,
                                             # show only if custom color is desired
                                             conditionalPanel(condition = "input.pchColor == 'custom'",
                                                              jscolorInput(inputId = "pchRgb",
                                                                        label = "Choose a color"))
                                      )
                                    )
                           ),##EndOf::Tab_5
                           
                           
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
                                              value = "Histogram"),
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
                                        a(href = "https://github.com/tzerk/RLumShiny/tree/master/inst/shiny/histogram", "See the code at GitHub!", target="_blank")
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
              tabPanel("Data set", dataTableOutput("dataset")),
              tabPanel("Central Age Model", dataTableOutput("CAM")),
              tabPanel("R plot code", verbatimTextOutput("plotCode"))
            )###EndOf::tabsetPanel
  )##EndOf::mainPanel
)##EndOf::shinyUI(pageWithSidebar)