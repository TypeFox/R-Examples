library(shiny)

shinyUI(fluidPage(theme = "bootstrap04.css",
  titlePanel(img(src="inegi_horizontal.png", align = "right")),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        'input.Etap01 === "Datos"',
        h2("Valores Extremos", align = "center"),
        tags$hr(),
        fileInput('file', 'Archivo de datos (texto o csv)',
                  accept=c('text/csv',
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        #numericInput("obs", "Primeros casos del archive:", 20),
        tags$hr(),
        
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
        tags$hr(),
        p("Registros en el archivo"),
        verbatimTextOutput("registros")
      ),
      
      conditionalPanel(
        'input.Etap01 === "CerosNA"',
        h3("Datos a cambiar"),
        tags$hr(),
        uiOutput("VNA"),
        tags$hr(),
        actionButton("updat1", "Poner NA's en la información")
      ),
      
      conditionalPanel(
        'input.Etap01 === "Extremos"',
        h3("Resultados"),
        tags$hr(),
        numericInput("ns", "Nivel de Significancia:", 0.975,
                     min = 0.6, max = 1,step = 0.005),
        br(),
        actionButton("updat2", "Obtener Extremos"),
        tags$hr(),
        radioButtons("En", "Datos que quiere ver:",
                     c("Variable Transformada" = "VarTra",
                       "ID de registros" = "reg",
                       "Valores Ext" = "valor",
                       "Manual" = "norm",
                       "Variables Ext" = "vars")),
        tags$hr(),
        p("Valores extremos encontrados"),
        verbatimTextOutput("val"),
        tags$hr(),
        p("Para obtener un archivo .CSV con los registros de los 
          valores extremos presionar en guardar."),
        downloadButton('DescarResum', 'Guardar')
      ),
      
      conditionalPanel(
        'input.Etap01 === "Graficas"',
        h4("Gráfica de densidad de una distribución normal univariada.", align = "center"),
        p("La visualización gráfica de un punto extremo multivariado con más de dos variables 
          no es posible realizarla, sin embargo, en una gráfica bidimensional de la curva de 
          densidad normal estándar, se puede dar una buena idea de la localización de cada 
          variable dentro del conjunto de registros que se están analizando."),
        p("En esta gráfica se presenta con líneas en guiones verticales de colores la 
          ubicación de cada variable de un registro identificado como punto extremo (PE) 
          multivariado, respecto de la curva normal estándar univariada representada en 
          rojo, con el mismo color de las líneas punteadas se proporciona el valor y 
          el nombre de las variables del PE. Se ha marcado con líneas discontinuas 
          verticales negras los umbrales de aceptación para facilitar la identificación de 
          las variables que de forma univariada son valores extremos."),
        tags$hr(),
        p("Para obtener un pequeño reporte"),
        downloadButton('downloadReport', 'Reporte')
      )
      ),

    
    
###########################################################

###########################################################
    mainPanel(#img(src="inegi_horizontal.png"),
      tabsetPanel(
        id = 'Etap01',
        tabPanel("Datos",
                  h4("Tabla de indices"),
                 dataTableOutput('tabla1')
                ),
        
        tabPanel("CerosNA",
                 h4("Tabla de indices"),
                 dataTableOutput('tabla2')
        ),
        
        tabPanel("Extremos",
                 conditionalPanel(condition = 'input.En === "VarTra"',
                                  h4("Tabla de Variables Transformadas"),
                                  tableOutput('tabla30')
                 ),
                 conditionalPanel(condition = 'input.En === "reg"',
                                  h4("Tabla de registros"),
                                  tableOutput('tabla31')
                 ),
                 conditionalPanel(condition = 'input.En === "valor"',
                                  h4("Tabla de Valores"),
                                  tableOutput('tabla32')
                 ),
                 conditionalPanel(condition = 'input.En === "norm"',
                                  h4("Tabla de Normal"),
                                  tableOutput('tabla33')
                 ),
                 conditionalPanel(condition = 'input.En === "vars"',
                                  h4("Tabla de Variables"),
                                  tableOutput('tabla34')
                 )
                 ),
        
        tabPanel("Graficas",
                 h4("Gráfica de extemos"),
                 numericInput("graf", "Gráfica que desea notar:", 1,min = 1),
                 tags$hr(),
                 plotOutput('plot1')
                )
              )
          )
  )
))
