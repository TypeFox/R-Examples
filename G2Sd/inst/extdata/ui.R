
library(shiny)

stat.choose <- c("arithmetic","geometric", "folk.ward", "all")
# Define UI for miles per gallon application
shinyUI(fluidPage(theme = "bootstrap3.css",
  
  # Application title
  headerPanel("G2Sd : Grain-size Statistics and Description of Sediment"),  
  sidebarPanel(
    
    span(p(img(src = "R_logo.png", height = 72, width = 72),strong("  G2Sd"), " and ", strong("Shiny"), " are available on CRAN"),style ="text-align:center"),
    
    wellPanel(
      h3(p(strong("Data"))),
      p(fileInput('file1', 'Choose XLS or XLSx File',accept=c('application/vnd.ms-excel', 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet', '.xlsx','.xls')),
        em(helpText("For example of databases,",a("click here",href = "http://regisgallon.wordpress.com/r-software/",target="_blank"))),
        numericInput('sheetindex', 'Sheet N°',1,min=1,max=99),    
        checkboxInput('header', 'Header', TRUE))
      
    ),
    #     ,
    wellPanel(
      h3(p(strong("Statistics "))),
      h5(em("Methods")),
      h4(strong(em(selectInput('method', '', stat.choose)))),
      h5(em("Results displayed")),
      h4(strong(em(selectInput('gran', '', c('Statistics','Index','Texture','All'))))),
      downloadButton('downloadData', 'Download')
    ),
    
    wellPanel(
      h3(p(strong("Graphics"))),
      em(helpText("To save graphics, right-click on the graphic",
                  "and select Save Picture As")),
      h5(em(strong("Histogram plot"))),
      numericInput('from','N° Station from',1),
      numericInput('to','To',1),
      checkboxInput('hist', 'Histogram', TRUE),
      checkboxInput('cum', 'Cumulative curve', TRUE),
      h5(em(strong("Distribution plot"))),
      selectInput('distritype', '', c('fine','large'))
    ),
    
    span(p("For more information, please contact ",  a(span("reg.gallon@gmail.com", style = "color:blue"), href="mailto:reg.gallon@gmail.com?subject=[G2Sd]%20Information"),style ="text-align:center")),  
    span(p(img(src = "mnhn.gif", height = 120, width = 120),img(src = "cnrs.gif", height = 72, width = 72)),style ="text-align:center")  
    
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Table", tableOutput('contents')),
      tabPanel("Statistics", tableOutput('stat')),
      tabPanel("Histogram", plotOutput("plot1",width="auto",height=1900)),
      tabPanel("Distribution", plotOutput("plot2",width="auto",height=800))
      
    ))
))

