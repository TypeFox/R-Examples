Defunciones<-function(){
  shinyAppDir(system.file("INEGI/Defun", package="Sofi"))
}

#.FdeDatos<-function(){tabPanel("Datos",
#                              h4("Tabla de c\u00f3digos"),
#                              dataTableOutput('tabla1'),
#                              tags$hr(),
#                              h4("Tabla de Datos"),
#                              tableOutput('tabla2')
#)}
