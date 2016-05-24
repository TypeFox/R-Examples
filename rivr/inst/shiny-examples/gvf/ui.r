shinyUI(pageWithSidebar(
  headerPanel("Gradually-varied flow"),
  sidebarPanel(
    h4("Model settings"),
    fluidRow(
      column(6,
        radioButtons(inputId = "Cm",
          label = "Units",
          c("US Customary" = 1.486, 
            "Standard International" = 1)
        )
      ),
      column(6,
        selectInput(inputId = "stepdist",
          label = "Step size",
          choices = c(1, 10, 50, 100, 200, 500),
          selected = 50
        )    
      )
    ),
    sliderInput(inputId = "totaldist", 
      label = "Distance to compute", 
      min = 500, 
      max = 5000, 
      value = 2000, 
      step = 500
    ),

    h4("Control section settings"),
    splitLayout(
      numericInput(inputId = "y0",
        label = "Depth",
        value = 1
      ),
      numericInput(inputId = "z0",
        label = "Elevation",
        value = 0
      ),
      numericInput(inputId = "x0",
        label = "Location",
        value = 0
      )
    ),
    
    h4("Channel settings"),
    fluidRow(
      column(8,
        sliderInput(inputId = "Q",
          label = "Flow rate",
          value = 10, 
          min   = 0, 
          max   = 1000,
          step = 5
        ),
        sliderInput(inputId = "n", 
          label = "Manning's roughness", 
          min = 0.001, 
          max = 0.2, 
          value = 0.04, 
          step = 0.0005,
          sep = ""
        ),
        sliderInput(inputId = "B", 
          label = "Bottom width", 
          min = 0, 
          max = 100, 
          value = 10, 
          step = 1
        ),
        sliderInput(inputId = "SS", 
          label = "Sideslope", 
          min = 0, 
          max = 10, 
          value = 0, 
          step= 0.5
        )
      ),
      column(4,
        radioButtons(inputId = "So",
          label = "Slope",
          c("0.00001", "0.00005",
            "0.0001", "0.0005",
            "0.001", "0.005",
            "0.01", "0.05",
            "0.1", "0.5"),
          selected = "0.001"
        )
      )
    )
  ),

  mainPanel(
    h3("This shiny app lets you explore gradually-varied flow profiles with the",
    code("rivr"), "package."),
    br(),
    h4("Don't know where to start? Try one of the profiles below."),
    tableOutput(outputId = "example_table"),
    tabsetPanel(
      tabPanel("Plot", br(), plotOutput(outputId = "main_plot", height = "500px")),
      tabPanel("Data", br(), fluidRow(
        column(3, tableOutput(outputId = "info_table")),
        column(9, tableOutput(outputId = "main_table"))
        )
      ),
      tabPanel("Credits",
      br(),
      h4("Created by", tags$b("Michael Koohafkan")),
      h4("Find me on"),
      h4(tags$ul(
        tags$li(tags$a(href = "https://github.com/mkoohafkan", "Github", 
        target = "_blank")),
        tags$li(tags$a(href = "http://stackoverflow.com/users/2375551/mikeck?tab=profile", 
          "StackOverflow", target = "_blank")),
        tags$li(tags$a(href = "https://www.linkedin.com/in/mkoohafkan", 
        "LinkedIn", target = "_blank"))
      )),
      h4("or check out my personal blog"),
      h4(tags$a(href = "http://www.hydroecology.net", "www.hydroecology.net", 
        target = "_blank"))
      )
    )
  )
))

