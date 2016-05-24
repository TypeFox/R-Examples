shinyUI(pageWithSidebar(
                        headerPanel("Histogram bin selection"),
                        sidebarPanel(
                                    selectInput(inputId = "n",
                                                label = "Sample size (n):",
                                                choices = c(10, 50, 100, 500, 1000, 5000, 10000),
                                                selected = 100),
                                    
                                    selectInput(inputId = "algorithm",
                                                label="Named algorithm",
                                                choices = c("DIY", "Sturges", "Scott", "FD"),
                                                selected = "Sturges"),
                      
                                    ## Display this only if the density is shown
                                    conditionalPanel(condition = "input.algorithm == 'DIY'",
                                                     textInput("n_bins", "Number of observations to view:", "10")
                                                     )
                                    ),
                        mainPanel(
                                  plotOutput(outputId = "main_plot", height = "300px"),
                                  tableOutput(outputId = "summary")
                                  )

                      
))
