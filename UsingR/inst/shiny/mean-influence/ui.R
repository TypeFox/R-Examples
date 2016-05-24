shinyUI(bootstrapPage(
                      h3("Effect on one point on the mean"),
                      p("Median marked in gray, 10 percent trimmed mean in short dash, mean in long dash. Moving the slider repositions the red point."),
                      selectInput(inputId = "n_points",
                                  label = "Number of points",
                                  choices = c(2,3,5,10,20),
                                  selected = 5),
                      
                      plotOutput(outputId = "main_plot", height = "300px"),
                      
                      ## Display this only if the density is shown
                      sliderInput(inputId = "nth",
                                  label = "Position of red point",
                                  min = 0, max = 50, value = 10, step = 1)
                      )
        
)
