

library(shiny)

shinyUI(fluidPage(
  titlePanel("Frequency Profile"),
  sidebarLayout(
    
    sidebarPanel(
      fileInput('file', 'Choose File',
                accept=c('.csv','.bin','.fpw')),
      tags$hr(),
      checkboxGroupInput("selected.behaviors", "Behaviors:",
                         c()),
      tags$hr(),
      checkboxInput('ggplot', 'ggplot', FALSE),
      helpText('Note: packages "ggplot2", "reshape2", and "grid" are required.'),
      checkboxInput('panel.in', 'Show left panel', TRUE),
      checkboxInput('panel.out', 'Show right panel', TRUE),
      checkboxInput('multiplot', 'Multiplots', FALSE),
      tags$hr(),
      radioButtons('which', 'Moving function',
                   c('Sum'='sum',
                     'Proportion'='proportion'),
                   'proportion'),
      numericInput("window", "Window length:", min=1, max=100, value=25, step=1),
      radioButtons('unit_length', 'length unit:',
                   c('%'='percent',
                     'bins'='bins')),
      numericInput("step", "Step:", min=1, max=100, value=1,
                   step=1),
      numericInput("resolution", "Resolution:", min=1, max=10, value=1,
                   step=1),
      textInput("units", "Time units:", value = "sec", width = NULL),
      numericInput("tick.every", "Tick every:", min=1, max=10, value=1,
                   step=1),
      numericInput("label.every", "Label every:", min=1, max=10, value=1,
                   step=1),
      tags$hr(),
      downloadButton('downloadData', 'Download Data'),
      tags$hr(),
      downloadButton('downloadPlotPDF', 'Download Plot as PDF'),
      downloadButton('downloadPlotPNG', 'Download Plot as PNG'),
      numericInput("graphWidth", "Width (inches):", min=1, max=20, value=10, step=1),
      numericInput("graphHeight", "Height (inches):", min=1, max=20, value=8, step=1)
    ),
    
    mainPanel(
      plotOutput("distPlot",height = "500px")
    )
    
  )
))


