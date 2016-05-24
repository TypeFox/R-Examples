shinyUI(pageWithSidebar(
                        headerPanel("Density plot kernel and bandwidth selection"),
                        sidebarPanel(
                                     h5("Sample selection:"),
                                     sliderInput(inputId = "n",
                                                 label = "Sample size:",
                                                 min=5,
                                                 max=1000,
                                                 value=100,
                                                 step=1
                                                ),
                                     selectInput(inputId = "family",
                                                label = "Distribution:",
                                                choices = c("Normal", "Exponential", "Symmetric, long-tailed", "Skew, long-tailed"),
                                                selected = "Normal"),

                                     h5("Kernel and bandwidth selection:"),
                                    selectInput(inputId = "kernel",
                                                label = "Kernel:",
                                                choices = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"),
                                                selected = "gaussian"),
                                    
                                    selectInput(inputId = "algorithm",
                                                label="Bandwidth algorithm",
                                                choices = c("DIY", "nrd0",  "nrd", "ucv", "bcv", "SJ"),
                                                selected = "nrd0"),
                      
                                    ## Display this only if the density is shown
                                    conditionalPanel(condition = "input.algorithm == 'DIY'",
                                                     sliderInput("bw", "Set the bandwidth:", min=0.01, max=2, value=1, step=0.01)
                                                     )
                                    ),
                        mainPanel(
                                  plotOutput(outputId = "main_plot", height = "300px"),
                                  tableOutput(outputId = "summary")
                                  )

                      
))
