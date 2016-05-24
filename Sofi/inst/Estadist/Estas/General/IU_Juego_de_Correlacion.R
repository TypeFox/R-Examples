fluidPage(
  
  tags$head(
    tags$style(HTML("
                    @import url('//fonts.googleapis.com/css?family=Lobster|Cabin:400,700');
                    "))
    ),
  
  titlePanel("Juego de Correlación"),
  
  fluidRow(
    column(width=4,
           selectInput("difficulty",label="Dificultad",choices=list("Fácil","Medio","Difícil")),
           h5(textOutput("score")),
           br(),
           checkboxGroupInput(inputId="options", label="Mostrar:",choices=list("Promedios", "línea de la desviación estándar", "Elipse")),
           br(),
           sliderInput("slider",label="La correlación entre X e Y es ...",min=-1,max=1,value=0,step=0.01),
           br(),
           p(textOutput("status1"),style="font-weight=500; color: #000000;"),
           h5(textOutput("status2"),style="font-weight=500; color: #00CC00;"),
           h5(textOutput("status3"),style="font-weight=500; color: #FF0000;"),
           br(),
           actionButton("submit","Enviar"),
           actionButton("newplot","Nuevo gráfico")
    ),
    
    column(width=8,
           plotOutput("plot1",width="450px",height="450px"))
    
  )
  
  
    )