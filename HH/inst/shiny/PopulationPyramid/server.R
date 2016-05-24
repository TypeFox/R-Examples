library(shiny)
library(HH)

## define data once outside the interaction
data(USAge.table)


## Define server logic for plot examples
shinyServer(function(input, output) {

  ## Reactive expression to compose a data frame containing all of the values
  Year <- reactive({
    ## cat(input$year, "\n")
    as.character(input$year)
  })

  ## Show the plot for specified Year
  output$USagePyramidPlot <- renderPlot({
    print(
    likert( ~ Female + Male, data=USAge.table[75:1, , Year()],
           main=Year(),
           scales=list(x=list(limits=c(-2500000, 2500000),
                         at=c(-2,-1,0,1,2)*1000000,
                         labels=c(2,1,0,1,2))),
           xlab="Population in Millions",
           ylab="Age")
      )
  })

  output$plotOutput <- renderUI(
    plotOutput("USagePyramidPlot", width="90%", height=input$px.height)
  )
})
