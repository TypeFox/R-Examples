library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("MFIaggr"),
  sidebarPanel(
    fileInput("input.file", "Choose CSV File (input should contain cycle data)",
              accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    checkboxInput("header", "Header", TRUE),
    radioButtons("csv.type", "Type of csv file",
                 c("Dec: dot (.), Sep: comma (;)" = "csv1",
                   "Dec: comma (.), Sep: semicolon (;)" = "csv2")),
    numericInput("cyc.col", "Column containing the cycle data:", 1,
                 min = 1, step = 1),
    checkboxInput("RSD", "Relative standard deviation:", FALSE),
    checkboxInput("rob", "Median and MAD:", FALSE),
    numericInput("llul.low", "Region of interest - lower border:", 1,
                 min = 1, step = 1),
    numericInput("llul.up", "Region of interest - lower border:", 10,
                 min = 1, step = 1),
    
    downloadButton("download.result", "Download MFIaggr results (with graphics)"),
    br(),
    downloadButton("download.table", "Download MFIaggr results (.csv)"),
    br(),
    p("Lost? Use button below to see an example:"),
    actionButton("run.example", "Run example")
  ),
  mainPanel(
    uiOutput("dynamic.tabset") 
  )
)
)