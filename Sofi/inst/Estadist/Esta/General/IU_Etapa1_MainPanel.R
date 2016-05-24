tabsetPanel(
  id = 'Etap01',
tabPanel("Datos",
         h4("Tabla de c\u00f3digos"),
         dataTableOutput('tabla1'),
         tags$hr(),
         h4("Tabla de Datos"),
         tableOutput('tabla2')),

tabPanel("Resumen",
         conditionalPanel(condition = 'input.En === "arc"',
                          fluidRow(
                            column(4,h4("Tabla de Resumen Archivo"),
                                   tableOutput('tabla5')),
                            column(7,offset = 1,h5("Cantidad de registros para la muestra:"),
                                   verbatimTextOutput("num5"))
                          )
         ),
         conditionalPanel(condition = 'input.En === "ecu"',
                          fluidRow(
                            column(5,h4("Tabla de Resumen Ecuaciones"),
                                   tableOutput('tabla51')),
                            column(5,offset = 1,h5("Cantidad de registros para la muestra:"),
                                   verbatimTextOutput("num51"))
                          )
         ),
         conditionalPanel(condition = 'input.En === "man"',
                          fluidRow(
                            column(5,h4("Tabla de Resumen Manual"),
                                   tableOutput('tabla52')),
                            column(5,offset = 1,h5("Cantidad de registros para la muestra:"),
                                   verbatimTextOutput("num52"))
                          )
         )),

tabPanel("Muestra",
         h4("Tabla con la muestra para la Etapa 1"),
         dataTableOutput('tabla4')
)   )