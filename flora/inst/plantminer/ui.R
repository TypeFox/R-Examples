shinyUI(fluidPage(
  tags$title("Plantminer - Brazilian Flora Checklist"),
  h1("Plantminer"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 checkboxInput("synonyms", label = "Replace synonyms", value = TRUE),
                 checkboxInput("suggest", label = "Correct misspelled names", value = TRUE),
                 checkboxInput("life.form", label = "Life form", value = FALSE),
                 checkboxInput("habitat", label = "Habitat", value = FALSE),
                 checkboxInput("vernacular", label = "Vernacular names", value = FALSE),
                 checkboxInput("states", label = "Occurrence", value = FALSE),
                 checkboxInput("establishment", label = "Establishment", value = FALSE),
                 sliderInput("distance", label = "Suggestion conservativeness (lower values are less conservative)",
                             min = 0, max = 1, value = 0.9),
                 tags$form(
                   tags$textarea(id="taxa", rows=16, cols=25, "Miconia albicans\nMyrcia lingua\nCofea arabica"),
                   tags$br(),
                   tags$input(type = "Submit"),
                   tags$i("(This may take a while)"))
    ),
    mainPanel(width = 9,
              h5("Data"),
              p("This application is an alternative front end for the",
              tags$a(href = "http://cran.r-project.org/package=flora", "flora"),
              "package for R. All data used here were kindly made available by the ",
              tags$a(href = "http://floradobrasil.jbrj.gov.br", "Brazilian Flora Checklist"),
              " project. Please cite them accordingly. This version of the application uses a database snapshot downloaded on 15 March 2015. Send your suggestions and report bugs to Gustavo Carvalho at gustavo.bio@gmail.com."
              ),
              h5("Usage"),
              p("Usage is simple: paste your taxa without authors in the textbox and hit submit. There is a download button below to export data as a quoted csv file."),
              dataTableOutput(outputId="contents"),
              downloadButton('downloadData', 'Download results in csv format')
    ))
))