# load dataset soils
library(exsic)

df = system.file("samples/exsic.csv", package="exsic")
data = read.exsic(df)[1:500,]

html.SBMG = file.path('SBMG.html')
html.NYBG = file.path('NYBG.html')
html.ASPT = file.path('ASPT.html')
html.PK   = file.path('PK.html')

exsic(data, html = html.SBMG, format = format.SBMG)
exsic(data, html = html.NYBG, format = format.NYBG)
exsic(data, html = html.ASPT, format = format.ASPT)
exsic(data, html = html.PK  , format = format.PK)


#############################################
shinyServer(function(input, output) {
  
  # Return the requested dataset
  formatInput <- reactive({
    switch(input$format,
           "SBMG" = html.SBMG,
           "NYBG" = html.NYBG,
           "ASPT" = html.ASPT,
           "PK"   = html.PK
    )
  })
  
  output$exsic <- renderText( {
    readLines(formatInput())
  })
  
  
  output$data <- renderTable({
    data
  })
  
})
            
            