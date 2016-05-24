shinyUI(fluidPage(
  titlePanel("Time series and stability vs instability"),
  
  sidebarLayout(position="right",
    sidebarPanel( numericInput("nsize", 
                               label="N Size (choose a value between 10 and 1000",
                               value=300,
                               min = 10,
                               max = 1000,
                               step=1),
                  selectInput("mean", "Mean Behavior",
                              c("Constant" = 1,
                                "Increasing" = 2,
                                "Decreasing" = 3,
                                "Oscillating"= 4)),
                  selectInput("var", "Variance Behavior",
                              c("Constant" = 1,
                                "Increasing" = 2,
                                "Decreasing" = 3,
                                "Oscillating"= 4))
    ),
                  
    mainPanel(p("A time series plot is considered \"stable\" only if the mean and variance
              are both constant across time."),
              plotOutput("plots"),
              textOutput("textvalue"),
              br(),
              br(),
              br(),
              p("Thanks to Omar Chavez for the implementation ideas for this page.")
              
              
              
              
    )          
  )
))