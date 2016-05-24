library(RLumShiny)
# Define UI
pageWithSidebar(
  
  # Application title
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
                                    tooltip(refId = "sumpos", attr = "for", text = "Position of the statistical summary. The keyword \"Subtitle\" will only work if no plot subtitle is used."),
                                    checkboxGroupInput(inputId = "stats",
                                                       label = "Parameters", 
                                                       selected = c("n","mean"),
                                                       choices = c("n" = "n",
                                                                   "Mean" = "mean",
                                                                   "weighted Mean" = "mean.weighted",
                                                                   "Median" = "median",
                                                                   #"weighted Median" = "median.weighted",
                                                                   "rel. Standard deviation" = "sdrel",
                                                                   "abs. Standard deviation" = "sdabs",
                                                                   "rel. Standard error" = "serel",
                                                                   "abs. Standard error" = "seabs",
                                                                   #"25 % Quartile" = "q25", #not implemented yet
                                                                   #"75 % Quartile" = "q75", #not implemented yet
                                                                   "Skewness" = "skewness", #not implemented yet
                                                                   "Kurtosis" = "kurtosis", #not implemented yet
                                                                   "Confidence interval" = "in.ci")),
                                    tooltip(refId = "stats", text = "Statistical parameters to be shown in the summary"),
                                    br(),
                                    div(align = "center", h5("Error range")),
                                    numericInput(inputId = "error",
                                                 label = "Symmetric error range (%)",
                                                 value = 10, min = 0, max = 100, step = 1),
                                    tooltip(refId = "error", text = "Symmetric error range in percent will be shown as dashed lines in the plot. Set error.range to 0 to void plotting of error ranges.")
                           ),##EndOf::Tab_2
                           
                           # Tab 3: input that refer to the plot rather than the data
                           tabPanel("DRT Details", 
                                    div(align = "center", h5("Experimental details")),
                                    numericInput(inputId = "dose", label = "Given dose (primary data set)", value = 2800),
                                    tooltip(refId = "dose", text = "Given dose used for the dose recovery test to normalise data. If only one given dose is provided this given dose is valid for all input data sets (i.e., values is a list). Otherwise a given dose for each input data set has to be provided (e.g., given.dose = c(100,200)). If no given.dose values are plotted without normalisation (might be useful for preheat plateau tests). Note: Unit has to be the same as from the input values (e.g., Seconds or Gray)."),
                                    numericInput(inputId = "dose2", label = "Given dose (secondary data set)", value = 3000),
                                    div(align = "center", h5("Preheat temperatures")),
                                    checkboxInput(inputId = "preheat", label = "Group values by preheat temperature", FALSE),
                                    tooltip(refId = "preheat", text = "Optional preheat temperatures to be used for grouping the De values. If specified, the temperatures are assigned to the x-axis."),
                                    conditionalPanel(condition = 'input.preheat == true', 
                                                     numericInput(inputId = "ph1", "PH Temperature #1", 180, min = 0),
                                                     numericInput(inputId = "ph2", "PH Temperature #2", 200, min = 0),
                                                     numericInput(inputId = "ph3", "PH Temperature #3", 220, min = 0),
                                                     numericInput(inputId = "ph4", "PH Temperature #4", 240, min = 0),
                                                     numericInput(inputId = "ph5", "PH Temperature #5", 260, min = 0),
                                                     numericInput(inputId = "ph6", "PH Temperature #6", 280, min = 0),
                                                     numericInput(inputId = "ph7", "PH Temperature #7", 300, min = 0),
                                                     numericInput(inputId = "ph8", "PH Temperature #8", 320, min = 0)
                                    )
                           ),##EndOf::Tab_3
                           
                           # Tab 3: input that refer to the plot rather than the data
                           tabPanel("Plot", 
                                    div(align = "center", h5("Title")),
                                    fluidRow(
                                      column(width = 6,
                                             textInput(inputId = "main", 
                                                       label = "Title", 
                                                       value = "DRT Plot")
                                      ),
                                      column(width = 6,
                                             textInput(inputId = "mtext", 
                                                       label = "Subtitle", 
                                                       value = "")
                                      )
                                    ),
                                    div(align = "center", h5("Boxplot")),
                                    checkboxInput(inputId = "boxplot", label = "Plot as boxplot", value = FALSE),
                                    tooltip(refId = "boxplot", text = "Optionally plot values, that are grouped by preheat temperature as boxplots. Only possible when preheat vector is specified."),
                                    div(align = "center", h5("Scaling")),
                                    sliderInput(inputId = "cex", 
                                                label = "Scaling factor",
                                                min = 0.5, max = 2, 
                                                value = 1.0, step = 0.1)
                           ),##EndOf::Tab_3
                           
                           # Tab 4: modify axis parameters
                           tabPanel("Axis",
                                    div(align = "center", h5("X-axis")),
                                    textInput(inputId = "xlab", 
                                              label = "Label x-axis",
                                              value = "# Aliquot"),
                                    # inject sliderInput from Server.R
                                    uiOutput(outputId = "xlim"),
                                    br(),
                                    div(align = "center", h5("Y-axis")),
                                    textInput(inputId = "ylab", 
                                              label = "Label y-axis",
                                              value = "Normalised De"),
                                    sliderInput(inputId = "ylim", label = "Range y-axis", 
                                                min = 0, max = 3, 
                                                value = c(0.75, 1.25), 
                                                step = 0.01)
                           ),##EndOf::Tab_4
                           
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
                                                              jscolorInput(inputId = "rgb2",
                                                                           label = "Choose a color"))
                                      )
                                    )
                           ),##EndOf::Tab_5
                           
                           # Tab xy: add and customize legend
                           tabPanel("Legend",
                                    div(align = "center", h5("Legend")),
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
                                    ),
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
                           ),##EndOf::Tab_xy
                           
                           # Tab xy: save plot as pdf, wmf or eps
                           tabPanel("Export",
                                    radioButtons(inputId = "fileformat", 
                                                 label = "Fileformat", 
                                                 selected = "pdf",
                                                 choices = c("PDF   (Portable Document Format)" = "pdf",
                                                             "SVG   (Scalable Vector Graphics)" = "svg",
                                                             "EPS   (Encapsulated Postscript)" = "eps")),
                                    textInput(inputId = "filename", 
                                              label = "Filename", 
                                              value = "DRT Plot"),
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
                                                   label = "Download plot")
                           ),##EndOf::Tab_
                           
                           # Tab xy: further information
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
                                        a(href = "https://github.com/tzerk/RLumShiny/tree/master/inst/shiny/doserecovery", "See the code at GitHub!", target="_blank")
                                    )#/div
                           )##EndOf::Tab_xy
               )
  ),
  
  
  # Show a plot of the generated distribution
  mainPanel(width = 7,
            # insert css code inside <head></head> of the generated HTML file:
            # allow open dropdown menus to reach over the container
            tags$head(tags$style(type="text/css",".tab-content {overflow: visible;}")),
            tags$head(includeCSS("www/style.css")),
            # divide output in separate tabs via tabsetPanel
            tabsetPanel(
              tabPanel("Plot", plotOutput(outputId = "main_plot", height = "400px")),
              tabPanel("Primary data set", fluidRow(column(width = 12, dataTableOutput("dataset")))),
              tabPanel("Secondary data set", dataTableOutput("dataset2"))
            )
  )
)