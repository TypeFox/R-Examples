sidebarLayout(
  sidebarPanel(
    conditionalPanel(
      'input.Etap02 === "Datos"',
      fileInput('Etapa2file1', 'Archivo de códigos (en dbf)',
                accept=c('.dbf')),
      helpText('Elegir las variables a utilizar'),
      uiOutput("Etap2CausaA"),
      uiOutput("Etap2Causa1"),
      uiOutput("Etap2Causa2"),
      tags$hr()
      
    ),
    ####
    conditionalPanel(
      'input.Etap02 === "Revisión"',
      ###
      helpText('Presionar evaluar para iniciar el proceso de obtención de casos a revisión'),
      actionButton("E2updat1", "Evaluar"),
      tags$hr(),
      h6("Registros para revisión:"),
      verbatimTextOutput("E2RegRev"),
      tags$hr(),
      checkboxInput("E2RevGua","Solo los Registros para revisión",value = T),
      downloadButton('E2DescarRev', 'Guardar'),
      tags$hr()
    ),
    ####
    conditionalPanel(
      'input.Etap02 === "Tablas"',
      ####
      helpText('Elegir la tabla que decea ver'),
      tags$hr(),
      radioButtons("Tab", "Tipo de tabla:",
                   c("Totales por caso" = "Caso",
                     "Errores a 3 Dígitos" = "Er3D",
                     "Errores a 4 Dígitos" = "Er4D",
                     "Caso 4 Mismo Capítulo" = "C4MC",
                     "Caso 4 Diferente Capítulo" = "C4DF"
                   )),
      conditionalPanel(condition = 'input.Tab === "Caso"',
                       helpText('Si desea guardar la tabla de totales por caso:'),
                       downloadButton('E2DescarCaso', 'Guardar'),
                       tags$hr()
      ),
      conditionalPanel(condition = 'input.Tab === "Er3D"',
                       helpText('Si desea guardar la tabla de Errores a 3 Dígitos:'),
                       downloadButton('E2DescarEr3D', 'Guardar'),
                       tags$hr()             
      ),
      conditionalPanel(condition = 'input.Tab === "Er4D"',
                       helpText('Si desea guardar la tabla de Errores a 4 Dígitos:'),
                       downloadButton('E2DescarEr4D', 'Guardar'),
                       tags$hr()               
      ),
      conditionalPanel(condition = 'input.Tab === "C4MC"',
                       numericInput("E2nMis", "Mínimo de Errores:", 1),
                       helpText('Si desea guardar frecuencias para mismo capítulo:'),
                       downloadButton('E2DescarC4MC', 'Guardar'),
                       tags$hr()               
      ),
      conditionalPanel(condition = 'input.Tab === "C4DF"',
                       numericInput("E2nDif", "Mínimo de Errores:", 1),
                       helpText('Si desea guardar frecuencias para diferente capítulo:'),
                       downloadButton('E2DescarC4DF', 'Guardar'),
                       tags$hr()
      ),
      tags$hr()
    ),
    ####
    conditionalPanel(
      'input.Etap02 === "Población"',
      ####
      helpText('Se requiere el tamaño total de la población, usar archivo de la Etapa 1 sección Resumen'),
      fileInput('Etapa2file2', 'Archivo de datos (texto o csv)',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      
      checkboxInput('header', 'Encabezado', TRUE),
      radioButtons('sep', 'Separado por:',
                   c(Coma=',',
                     Puntoycoma=';',
                     Tabulador='\t'),
                   ','),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"'),
      uiOutput("Etap2Pobla"),
      tags$hr()
    ),
    #####
    conditionalPanel(
      'input.Etap02 === "Limites"',
      #####
      helpText('Intervalos de confianza'),
      actionButton("E2Inter", "Calcular"),
      helpText('Elegir la tabla que decea ver'),
      tags$hr(),
      radioButtons("E2Lim", "Límites o gráficos que desea mostrar:",
                   c("Intervalos de confianza a 3 dígitos" = "IC3D",
                     "Intervalos de confianza a 4 dígitos" = "IC4D",
                     "Gráfico a 3 Dígitos" = "Gr3D",
                     "Gráfico a 4 Dígitos" = "Gr4D"
                   )),
      conditionalPanel(condition = 'input.E2Lim === "IC3D"',
                       helpText('Guardar intervalos de confianza a 3 dígitos:'),
                       sliderInput(inputId = "E2ErrorI3",
                                   label = "Valor para Error (alfa):",
                                   min = .01, max = .2, value = .05, step = 0.01),
                       downloadButton('DescarE2Inter3', 'Guardar'),
                       tags$hr()
      ),
      conditionalPanel(condition = 'input.E2Lim === "IC4D"',
                       helpText('Guardar intervalos de confianza a 4 dígitos:'),
                       sliderInput(inputId = "E2ErrorI4",
                                   label = "Valor para Error (alfa):",
                                   min = .01, max = .2, value = .05, step = 0.01),
                       downloadButton('DescarE2Inter4', 'Guardar'),
                       tags$hr()             
      ),
      conditionalPanel(condition = 'input.E2Lim === "Gr3D"',
                       helpText('Si desea guardar la Gráfico a 3 Dígitos:'),
                       uiOutput("Etap2Int3"),
                       tags$hr()               
      ),
      conditionalPanel(condition = 'input.E2Lim === "Gr4D"',
                       helpText('Si desea guardar la Gráfico a 4 Dígitos:'),
                       uiOutput("Etap2Int4"),
                       tags$hr()               
      ),
      tags$hr(),
      #downloadButton('downloadReportE2', 'Reporte')
      
      radioButtons('format_1', 'Formato del documento', c('HTML', 'Word'),
                   inline = TRUE),
      downloadButton('DescarE2Repot', 'Reporte 1')
      
    )
    
  ),
  ###
  mainPanel(
    tabsetPanel(
      id = 'Etap02',
      ###
      tabPanel("Datos",    
               h4("Tabla de Datos"),
               dataTableOutput('Etapa2Tabla1')
      ),
      
      tabPanel("Revisión",
               h4("Tabla de Revisión"),
               dataTableOutput('Etapa2Tabla2')
      ),
      
      tabPanel("Tablas",
               conditionalPanel(condition = 'input.Tab === "Caso"',
                                h4("Tabla de Totales por Caso"),
                                tableOutput('Etapa2Tabla31')
                                #dataTableOutput('Etapa2Tabla34')
               ),
               conditionalPanel(condition = 'input.Tab === "Er3D"',
                                h4("Tabla de Errores a 3 Dígitos"),
                                tableOutput('Etapa2Tabla32')
               ),
               conditionalPanel(condition = 'input.Tab === "Er4D"',
                                h4("Tabla de Errores a 4 Dígitos"),
                                tableOutput('Etapa2Tabla33')
               ),
               conditionalPanel(condition = 'input.Tab === "C4MC"',
                                h4("Frecuencias para mismo capítulo"),
                                tableOutput('Etapa2TablaMis')
               ),
               conditionalPanel(condition = 'input.Tab === "C4DF"',
                                h4("Frecuencias para diferente capítulo"),
                                tableOutput('Etapa2TablaDif')
               )
      ),
      
      tabPanel("Población",
               h4("Tabla para Análisis"),
               tableOutput('Etapa2TablaTot')
      ),
      
      tabPanel("Limites",
               conditionalPanel(condition = 'input.E2Lim === "IC3D"',
                                h4("Intervalos de confianza a 3 dígitos"),
                                tableOutput('Etapa2Inter3')
               ),
               conditionalPanel(condition = 'input.E2Lim === "IC4D"',
                                h4("Intervalos de confianza a 4 dígitos"),
                                tableOutput('Etapa2Inter4')
               ),
               conditionalPanel(condition = 'input.E2Lim === "Gr3D"',
                                h4("Gráfico a 3 Dígitos"),
                                plotOutput('E2GI3')
               ),
               conditionalPanel(condition = 'input.E2Lim === "Gr4D"',
                                h4("Gráfico a 4 Dígitos"),
                                plotOutput('E2GI4')
               )
               
      )#tabPanel("Limites"
      
    )#tabsetPanel(id = 'Etap02'
  )#mainPanel(
)#sidebarLayout(