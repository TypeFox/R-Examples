library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("Comparision of digital PCR experiments"),
  sidebarPanel(
    numericInput("m_1", "Number of molecules (experiment 1):", 90, min = 1, max = 200,
                 step = 1),
    numericInput("nr_1", "Number of runs (experiment 1):", 2, min = 1, max = 5, step = 1),
    numericInput("m_2", "Number of molecules (experiment 2):", 10, min = 1, max = 200,
                 step = 1),
    numericInput("nr_2", "Number of runs (experiment 2):", 2, min = 1, max = 5, step = 1),
    numericInput("m_3", "Number of molecules (experiment 3):", 150, min = 1, max = 200,
                 step = 1),
    numericInput("nr_3", "Number of runs (experiment 3):", 2, min = 1, max = 5, step = 1),
    includeMarkdown("readme.md")
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Input data", verbatimTextOutput("input_data_summ"),
               tableOutput("input_data")),
      tabPanel("Results of test", verbatimTextOutput("test"),
               plotOutput("test_plot"), verbatimTextOutput("test_summ"),
               plotOutput("test_plot_agg"))
    )
  )
))

