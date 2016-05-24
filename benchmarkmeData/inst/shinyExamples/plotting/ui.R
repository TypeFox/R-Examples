library("shiny")

get_tests = function() {
  ui_env = new.env()
  data(past_results, package="benchmarkmeData", envir=ui_env)
  tests = unique(ui_env$past_results$test_group)
  results = get("results", envir=benchmarkmeData:::.bme_env)
  if(!is.null(results)) {
    tests = tests[tests %in% unique(results$test_group)]
  }
  tests
}

opts = c("None", "Byte Optimised", "R Version", "OS", "BLAS Optimised")
ui = fluidPage(
  titlePanel("Benchmark data"),
  br(),
  plotOutput("plot"),
  hr(),
  fluidRow(
    column(3, 
           wellPanel(selectInput("facet_x", "Facet x", opts)
           )),
    column(3, 
           wellPanel(selectInput("facet_y", "Facet y", opts)
           )),
    column(3, 
           wellPanel(selectInput("test", "Benchmark test", get_tests())
           ))
   )
)