library(shiny)
library(quipu)

data(potato.quipu)
dat = potato.quipu

# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output) {
  
  
  # Generate a summary of the dataset
  output$quipuPlot <- renderPlot({
    # print(input$acc_id)
    rquipu(dat, 
           a.subset = input$acc_id, 
           layout = input$layout,
           species.name = input$speciesName,
           set.name = input$setName,
           id.label = input$idLabel, 
           node.size = c(as.numeric(input$nodeG1),
                         as.numeric(input$nodeG2),
                         as.numeric(input$nodeG3),
                         as.numeric(input$nodeG4)),
           col.node = c(input$colorG1,
                        input$colorG2,
                        input$colorG3,
                        input$colorG4
                        )
    )
    #plot(1:10,1:10)
  })
  
})
