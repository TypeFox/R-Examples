#1) global info
library(shiny)
library(sp)
library(rgdal)

colnames = ''
colnames2= ''

# 2) Define server 
#---------------------------------------------------------------------
server <- function(input,output,session) {
  
  dataIn <- reactive({ 
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    
    data <- read.csv(inFile$datapath, header=input$header)
  })
  
  datasetInput <- reactive({
    if(is.null(dataIn()))
      return()
    data <- dataIn()
  })
  
  observe({
    updateSelectInput(session,"Time",
                      choices = c('',colnames(datasetInput())))
    
    updateSelectInput(session,'Lat',
                      choices= c('',colnames(datasetInput())))
    
    updateSelectInput(session,'Long',
                      choices= c('',colnames(datasetInput())))
    
    updateSelectInput(session,'Cov1',
                      choices= c('',colnames(datasetInput())))
    
    updateSelectInput(session,'Cov2',
                      choices= c('',colnames(datasetInput())))
    
    updateSelectInput(session,'Cov3',
                      choices= c('',colnames(datasetInput())))
  })
  
  #Data View 
  output$view <- renderTable({
    head(datasetInput(),n=6L)
  })
  
  #Time Information
  output$time <- renderPrint({
    dataset <- datasetInput()
    if (input$Time=="") return()
    anyNA(dataset[,input$Time])
  })
  
  output$timeClass <- renderPrint({
    dataset <- datasetInput()
    if(input$Time=="") return()
    class(dataset[,input$Time])
    #if(class(dataset[,input$Time])=='numeric') print('Time is numeric.')
    #else print('Time must be numeric')
  })
  
  #Lat/Long Information
  output$lat <- renderPrint({
    dataset <- datasetInput()
    if(input$Lat=="") return()
    class(dataset[,input$Lat])
  })
  
  output$long <- renderPrint({
    dataset <- datasetInput()
    if(input$Long=="") return()
    class(dataset[,input$Long])
  })
  
  #projection Information
  
  newData = reactive(function(){
    df <- datasetInput()
    if(input$Long=="") return()
    if(input$Lat=="") return()
    x <- df[,input$Long]
    y <- df[,input$Lat]
    data.frame(df,x,y)
  })
  
  output$new.table <- renderTable({
    if (is.null(newData())) {return()}
    head(newData(),n=6L)
  })
  
  output$proj <- renderText({
    if(is.null(newData())) return()
    newdata <- newData()
    coordinates(newdata) <- ~x+y
    proj4string(newdata) <- CRS("+proj=longlat")
    proj4string(newdata)
    
  })
  
  output$epsg <- renderText({
    epsg <- paste0('+init=epsg:',input$epsg)
    epsg
  })
  
  projData <- reactive(function(){
    if(is.null(newData())) {return()}
    if(input$projection==TRUE) stop ('')
    projdata <- newData()
    #if(is.na(projdata[,input$Lat])) stop ('')
    coordinates(projdata) <- ~x+y
    proj4string(projdata) <- CRS("+proj=longlat")
    
    if(is.null(input$epsg)) {return()}
    epsgProj <- paste0('+init=epsg:',input$epsg)
    spTransform(projdata,CRS(epsgProj))
  
    })
  
  output$projection <- renderText({
    if(is.null(projData())) {return()}
    class(projData())
  })
  
  #Covariates
  output$cov1 <- renderPrint({
    dataset <- datasetInput()
    if(input$Cov1=="") return()
    anyNA(dataset[,input$Cov1])
  })
  
  output$cov2 <- renderPrint({
    dataset <- datasetInput()
    if(input$Cov2=="") return()
    anyNA(dataset[,input$Cov2])
  })
  
  output$cov3 <- renderPrint({
    dataset <- datasetInput()
    if(input$Cov3=="") return()
    anyNA(dataset[,input$Cov3])
  })
  
  #Export new table
  output$downloadTable <- downloadHandler(        
    filename = function() { paste(input$table_name, '.csv', sep='')},
    content = function(file) {
      if(input$projection==TRUE) {write.csv(datasetInput(),file)}
      else 
        write.csv(projData(), file, row.names=FALSE)
    })
}



# 3) Define UI 
#--------------------------------------------------------------
ui <- fluidPage(
  titlePanel("CRAWL Data Check"),
  
  sidebarLayout(
    
    sidebarPanel(
      fileInput('file1', h4('1. Choose CSV File'),
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain','.csv')),
      #tags$hr(),
      checkboxInput('header', 'Header', TRUE),
      #radioButtons('sep', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),',',
      #             inline=TRUE),
      
      h4('2.Projection Information'),
      checkboxInput('projection', 'Are the data already projected?', FALSE),
      textInput('epsg',label='Define Projection: enter epsg value'),
      
      h4('3. Select Data Columns'),
      selectInput('Time','Select Time Column from Data',
                  choices=c(Choose='',colnames),selectize=FALSE,selected=NULL),
      
      selectInput('Lat', 'Select Latitude Column From Data',
                  choices=c(Choose='',colnames2),selectize=FALSE),
      
      selectInput('Long', 'Select Longitude Column From Data',
                  choices=c(Choose='',colnames2),selectize=FALSE),
      
      selectInput('Cov1','Select Covariates (if any) From Data',
                  choices=c(Choose='',colnames),selectize=FALSE),
      
      selectInput('Cov2','',
                  choices=c(Choose='',colnames),selectize=FALSE),
      
      selectInput('Cov3','',
                  choices=c(Choose='',colnames),selectize=FALSE),
      
      #submitButton('Update Columns'),
      
      h4('4.Export New Dataframe'),
      textInput('table_name','Name your new dataframe'),
      downloadButton('downloadTable','Save data table to .csv')
      
    ),
    
    mainPanel(
      h4("Data Preview"),
      tableOutput('view'),
      
      h4("Time"),
      h5('Does Time Contain NAs?'), verbatimTextOutput('time'),
      h5("Time Class?"), verbatimTextOutput('timeClass'),
      helpText('Note: Time must be either numeric or POSIX. If your data are not in either of these classes, please convert it prior to implementing a crawl model.'),
      
      h4('Lat/Long'),
      h5('Location Class'),verbatimTextOutput('lat'),verbatimTextOutput('long'),
      helpText('Note: This version of the data check app does not project data that has missing Lat/Long values. If you have missing coordinate values, please see the Harbor Seal vignette for details on how to proceed.'),
      
      h4("Covariates"),
      h5('Do Covariates Contain NAs?'),
      verbatimTextOutput('cov1'),
      verbatimTextOutput('cov2'),
      verbatimTextOutput('cov3'),
      
      h4('Projection Information'),
      h5('Current Projection:'),
      verbatimTextOutput('proj'),
      h5('epsg Projection:'),
      verbatimTextOutput('epsg'),
      verbatimTextOutput('projection'),
      
      h4('Updated Dataframe Preview'),
      tableOutput('new.table')
      
    )
  )
)

#-----------------------------------------------------
shinyApp(ui=ui,server = server)

