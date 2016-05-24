#library(shiny)
#library(shinythemes)
#tabPanelAbout <- source("about.r")$value

options(shiny.deprecation.messages=FALSE)
shinyUI(navbarPage("Sofi",#theme = shinytheme("readable"),
####_____________________________________________________________
#___________________________Etapa 1______________________________
####________________________________________________________________                   
#tabPanel("Calculadora de Distribución",icon = icon("random"),
#         source(system.file("Estadist/Estas/General/IU_Distribution_Calculator.R", package="Sofi"),local=T,encoding="UTF-8")$value
#),
tabPanel("Teorema del límite central",icon = icon("random"),
         source(system.file("Estadist/Estas/General/IU_Teorema_limite_central.R", package="Sofi"),local=T,encoding="UTF-8")$value
),
navbarMenu("Estadística 1",
           tabPanel("Diagnosticos",
                    source(system.file("Estadist/Estas/General/IU_Diagnosticos_regresion_lineal.R", package="Sofi"),local=T,encoding="UTF-8")$value),
           tabPanel("Juego de Correlación",
                    source(system.file("Estadist/Estas/General/IU_Juego_de_Correlacion.R", package="Sofi"),local=T,encoding="UTF-8")$value))
#tabPanel("Diagnosticos",icon = icon("random"),
#         source(system.file("Estadist/Estas/General/IU_Diagnosticos_regresion_lineal.R", package="Sofi"),local=T,encoding="UTF-8")$value
#),
#tabPanel("Juego de Correlación",icon = icon("random"),
#         source(system.file("Estadist/Estas/General/IU_Juego_de_Correlacion.R", package="Sofi"),local=T,encoding="UTF-8")$value
#),
#tabPanel("Estadística 1",icon = icon("random"),
#         source(system.file("Estadist/Esta/General/IU_Distribuciones_Var_Aleat.R", package="Sofi"),local=T,encoding="UTF-8")$value
#),

####_____________________________________________________________
#___________________________Etapa 2 y 3______________________________
####________________________________________________________________

####_____________________________________________________________
#___________________________Etapa 4 y 5______________________________
####________________________________________________________________

))
