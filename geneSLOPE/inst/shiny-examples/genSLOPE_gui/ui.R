library(shiny)

shinyUI(fluidPage(

  tags$head(
    tags$style(HTML("
      .shiny-output-error-validation {
        color: red;
      }
    "))
  ),

  titlePanel(h1("GWAS with SLOPE", align = "center", style="color:#317eac"), windowTitle = "GWAS with SLOPE"),

  column(4, wellPanel(
    fileInput("fileY",
              label = 'Choose file with phenotype',
              accept=c('.phe')),
    checkboxInput('header', 'Header', FALSE),
    radioButtons('sep', label = 'Separator',
                 c(Space=' ', Semicolon=';', Comma=',', Tab='\t'), ' '),
    tags$hr()
    ),
    conditionalPanel(
      condition = "output.phenotypeOk",
      fileInput("file", label = 'Choose .raw file with snps',
                accept=c('*.raw','text/csv'), multiple = TRUE),
      helpText("Max size is 300Mb"),
      fileInput("map.file", label = 'Choose .map file with snp info',
                accept=c('*.map'), multiple = TRUE),
      h5("After setting p-value", align="center"),
      h4("CLICK 'Run'", align = "center", style = "color:blue; font-si16pt"),
      actionButton("go", "Run", icon = icon("youtube-play"), width = "100%"),
      sliderInput("pValCutoff",
                  label = "P value cutoff",
                  min = 0, max = 1,
                  value = 0.05),
      helpText("P-value for initial snp screening"),
      sliderInput("fdr",
                  label = "FDR",
                  min = 0, max = 1,
                  value = 0.1),
      helpText("False Discovery Rate for SLOPE"),
      sliderInput("rho",
                  label = "Correlation",
                  min = 0, max = 1,
                  value = 0.3),
      helpText("Correlation threshold for snps clustering"),
      downloadButton('downloadData', 'Download important snps')
    ),
    tags$hr(),
    textOutput("phenotypeOk")
  ),

  column(5,
         conditionalPanel("output.phenotypeOk",
                          h4("Head of phenotype data"),
                          tableOutput("summary"),
                          h4("Clumping procedure summary"),
                          verbatimTextOutput("clumpSummary"),
                          h4("GWAS results. Important snps"),
                          plotOutput("slopePlot")
         )
  )
))

