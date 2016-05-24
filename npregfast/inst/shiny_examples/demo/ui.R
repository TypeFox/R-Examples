detach("package:npregfast")
library(shiny)
#library(shinyjs)
library(miniUI)
#library(wesanderson)
library(npregfast)

shinyUI(fluidPage(
  title = "Demo of npregfast",
  tags$head(includeCSS(file.path('www', 'style.css'))),   
  shinyjs::useShinyjs(),
  
  fluidRow(id = "title-row",
           column(12,
                  h1("Demo of",em(a("npregfast", href = "https://github.com/sestelo/npregfast"))),
                  h4("Example with", a("barnacle", href = "https://github.com/sestelo/npregfast/blob/master/man/barnacle.Rd")," data set"),
                  div("Created by", a("Marta Sestelo", href = "http://sestelo.github.io"),
                      "and", a("Nora M. Villanueva",href = "http://noramvillanueva.github.io"), HTML("&bull;"),
                      "Code on", a("GitHub", href = "https://github.com/sestelo/shiny_npregfast/")
                  )
           )
  ),
  
  
  
  div(id = "loading-content", h2("Loading...")),
  
  fluidRow(id = "app-content",
           column(2, wellPanel(
             class = "settings",
             h4(class = "settings-title", "Estimation"),
             
             selectInput(inputId = "type", 
                         label = "Factor-by-curve interaction?",
                         choices = c("Without" = "without", "With" = "with")),
             
             selectInput(inputId = "kernel",
                         label = "Choose a kernel:",
                         choices = c("Epanechnikov" = "epanech", 
                                     "Gaussian" = "gaussian",
                                     "Triangular" = "triang")),
             
             selectInput(inputId = "poly",
                         label = "Polynomial degree:",
                         choices = c(1, 
                                     2,
                                     3),
                         selected = 3),
             
             radioButtons(inputId = "selband",
                          label = "Bandwidth selection:",
                          choices = c("Cross-validation" = "cv", 
                                      "Manual" = "man"),
                          selected = "cv"),
             
             conditionalPanel(
               condition = "input.selband == 'man'",
               sliderInput(inputId = "band",
                           label = "Bandwidth selection:",
                           min = 0,
                           max = 1,
                           value = 0.5,
                           step = 0.1, 
                           ticks = TRUE,
                           animate = TRUE))
             
           )),
           
           
           
           
           
           column(2, wellPanel(
             class = "settings",
             h4(class = "settings-title", "Graphical"),
             
             conditionalPanel(
               condition = "input.poly == 1",
               checkboxGroupInput(inputId = "der1",
                                  label = "Output:",
                                  choices = c("Conditional mean" = '0'),
                                  selected = '0')),
             
             conditionalPanel(
               condition = "input.poly == 2",
               checkboxGroupInput(inputId = "der2",
                                  label = "Output:",
                                  choices = c("Conditional mean" = '0', 
                                              "First derivative" = '1'),
                                  selected = '0')),
             
             conditionalPanel(
               condition = "input.poly == 3",
               checkboxGroupInput(inputId = "der3",
                                  label = "Output:",
                                  choices = c("Conditional mean" = '0', 
                                              "First derivative" = '1',
                                              "Second derivative" = '2'),
                                  selected = '0')),
             
             
             
             
             div(id = "marginal-settings",
                 shinyjs::colourInput("colmu", "Line color", "#D67236", 
                                      showColour = "background",
                                      palette = "limited",
                                      allowedCols = unlist(wesanderson::wes_palettes),
                                      allowTransparent = FALSE),
                 
                 
                 shinyjs::colourInput("colci", "CI color", "#5B1A18", 
                                      showColour = "background",
                                      palette = "limited",
                                      allowedCols = unlist(wesanderson::wes_palettes),
                                      allowTransparent = FALSE)
             ),
             
             conditionalPanel(
               condition ="input.poly == 1 & input.der1[0] == '0'||input.poly == 2 & input.der2[0] == '0'||input.poly == 3 & input.der3[0] == '0'",
               checkboxInput("show_points", "Show data points", TRUE),
               conditionalPanel(
                 condition ="input.show_points == true",
                 shinyjs::colourInput("pcol", "Points color", "#899DA4", 
                                      showColour = "background",
                                      palette = "limited",
                                      allowedCols = unlist(wesanderson::wes_palettes),
                                      allowTransparent = FALSE)
               )
             )
           )),
           
           
           
           column(8,
                  plotOutput("distPlot", 
                             height = "500px", 
                             width = "100%",
                             click = "plot1_click",
                             brush = brushOpts(id = "plot1_brush"),
                  ),
                  
                  miniButtonBlock(
                   
                    actionButton("exclude_toggle", "Toggle points",
                                 icon = icon("fa fa-codiepie", class = "fa-1x")),
                    
                    actionButton("exclude_reset", "Reset", 
                                 icon = icon("fa fa-refresh", class = "fa-1x")),
                    actionButton(inputId ="info_btn", 
                                 label = "Info", 
                                 icon = icon("fa fa-info-circle", class = "fa-1x"))
                    
                  )
                  
                  
                 # includeMarkdown("plot_shiny.md")
                  
                  
                  
           )
           
  )
))








