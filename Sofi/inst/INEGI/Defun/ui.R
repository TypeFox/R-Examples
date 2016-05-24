#library(shiny)
library(shinythemes)
options(shiny.deprecation.messages=FALSE)
shinyUI(navbarPage("Defunciones",theme = shinytheme("cerulean"),
####_____________________________________________________________
#___________________________Etapa 1______________________________
####________________________________________________________________                   
tabPanel("Muestra",icon = icon("random"),
    #source(system.file("INEGI/Defun/General/IU_Etapa1.R", package="Sofi"),local=T,encoding="UTF-8")$value
    #source(system.file("INEGI/Defun/General/IU_Etapa1.R", package="Sofi"),local=T)$value
    source('./General/IU_Etapa1.R',local=T,encoding="UTF-8")$value
),
####_____________________________________________________________
#___________________________Etapa 2______________________________
####________________________________________________________________
tabPanel("Revisión",
    source(system.file("INEGI/Defun/General/IU_Etapa2.R", package="Sofi"), local=T, encoding="UTF-8")$value
),#tabPanel("Etapa 2 y 3",
####_____________________________________________________________
#___________________________Etapa 3______________________________
####________________________________________________________________
tabPanel("Experto",
#         source(system.file("INEGI/Defun/General/IU_Etapa2.R", package="Sofi"),local=T)$value
sidebarLayout(
  
  sidebarPanel(
    h5("Códigos sugeridos por el codificador experto."),
  
    uiOutput("Et3_Num_reg"),
    uiOutput("Et3_Causa"),
    
    checkboxGroupInput("Cod_Cor", "Radio buttons:",
                 c("label 1" = "option1",
                   #"label 2" = "option2",
                   "label 3" = "option3",
                   "label 4" = "option4")#,
                 #selected = ""
                 ),
    textInput("Et3_inText",  "Código:", value = ""),
    
    actionButton("Bot_RAM", "Guardar"),
    actionButton("Bot_Ant", "Anterior"),
    actionButton("Bot_Sig", "Siguiente"),
    actionButton("Bot_Guar", "Grabar"),
    
    h5("COD_SEL:"),
    verbatimTextOutput("E3COD"),
    verbatimTextOutput("E3CapTex"),
    #wellPanel(h5("COD_SEL:"),h2(verbatimTextOutput("E3COD"))),
    checkboxGroupInput("E3_Cod_Int", "Ver código:",
                       c("Internet" = "intr",
                         "Tabla" = "tabl",
                         "Tabla Códigos" = "tabl_Codi"), inline = TRUE
                         #"label 3" = "option3",
                         #"label 4" = "option4"
                         #,
                       #selected = ""
    )
    
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    fluidPage(
      fluidRow(
        column(4,
               actionButton("E3_sal", "Salir"),align = "center"),
        column(8,verbatimTextOutput("E3Nul"),align = "center")
      ),
               fluidRow(
                 column(width = 4,verbatimTextOutput("E3DatIde")),
                 column(width = 4,verbatimTextOutput("E3Folea1")),
                 column(width = 4,verbatimTextOutput("E3Folea2"))
                        ),
               fluidRow(
                 column(width = 8,h4("Causa de la defunción",align = "center"),
                        verbatimTextOutput("E3Desc1")),
                 column(width = 2,h4("Duración",align = "center"),
                        verbatimTextOutput("E3DURATION1")),
                 column(width = 2,h4("CIE-10",align = "center"),
                        verbatimTextOutput("E3t_CoA"))
               ),
               fluidRow(
                 column(width = 8,
                        verbatimTextOutput("E3Desc2")),
                 column(width = 2,
                        verbatimTextOutput("E3DURATION2")),
                 column(width = 2,
                        verbatimTextOutput("E3t_CoB"))
               ),
               fluidRow(
                 column(width = 8,
                        verbatimTextOutput("E3Desc3")),
                 column(width = 2,
                        verbatimTextOutput("E3DURATION3")),
                 column(width = 2,
                        verbatimTextOutput("E3t_CoC"))
               ),
               fluidRow(
                 column(width = 8,
                        verbatimTextOutput("E3Desc4")),
                 column(width = 2,
                        verbatimTextOutput("E3DURATION4")),
                 column(width = 2,
                        verbatimTextOutput("E3t_CoD"))
               ),
               fluidRow(
                 column(width = 8,
                        verbatimTextOutput("E3Desc5")),
                 column(width = 2,
                        verbatimTextOutput("E3DURATION5")),
                 column(width = 2,
                        verbatimTextOutput("E3t_CoI")
                        )
               )
      #DT::dataTableOutput("Etapa3Tabla1", width = 900)
  )#fluidPage
  
  ),#mainPanel
  position = "right"
),#sidebarLayout
conditionalPanel("input.E3_Cod_Int == 'intr'",
                 tags$iframe(src="https://eciemaps.mspsi.es/ecieMaps/browser/index_10_2008.html", width = 1200, height=700)#, width = 1200, height=600
                 
                 #tags$iframe(src="https://eciemaps.mspsi.es/ecieMaps/browser/index_10_2008.html", width = 1200, height=600)
                 #tags$iframe(src=paste0("https://eciemaps.mspsi.es/ecieMaps/browser/index_10_2008.html#search=",verbatimTextOutput("E3Link")), width = 1100, height=600) 
                 #cat("Text",verbatimTextOutput("E3Link"))
),
conditionalPanel("input.E3_Cod_Int == 'tabl'",
                 dataTableOutput("Etapa3Tabla1")#, width = 900
),
conditionalPanel("input.E3_Cod_Int == 'tabl_Codi'",
                 dataTableOutput("Etapa3Tabla2")#, width = 900
)
),#tabPanel("Etapa 2 y 3",


####_____________________________________________________________
#___________________________Etapa 4 y 5______________________________
####________________________________________________________________
tabPanel("Final",
         #fluidPage(
    source(system.file("INEGI/Defun/General/IU_Etapa4.R", package="Sofi"),local=T,encoding="UTF-8")$value
)#tabPanel("Etapa 4 y 5",

))
