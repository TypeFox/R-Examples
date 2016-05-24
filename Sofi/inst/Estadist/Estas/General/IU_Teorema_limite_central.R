pageWithSidebar(
  
  headerPanel("Teorema del límite central"),
  
  sidebarPanel(
    radioButtons("dist", "Distribución (población):",
                 list("Normal" = "rnorm",
                      "Uniforme" = "runif",
                      "Sesgada a la derecha" = "rlnorm",
                      "Sesgada a la izquierda" = "rbeta")),
    br(),
    
    uiOutput("mu"),
    uiOutput("sd"),
    uiOutput("min"),
    uiOutput("max"),
    uiOutput("skew"),
    
    sliderInput("n", 
                "Tamaño de las muestras:", 
                value = 30,
                min = 2, 
                max = 100),
    br(),
    
    sliderInput("k", 
                "Número de muestras:", 
                value = 20,
                min = 10, 
                max = 100),
    br(),
    helpText("Para saber más sobre el teorema del límite central leer ", 
             a(href="http://es.wikipedia.org/wiki/Teorema_del_l%C3%ADmite_central", target="_blank", "esto,"),
             " Y recuerda visitar la página de", a(href="http://www.inegi.info/sofi", target="_blank", "Sofi.")),
    helpText(a(href="https://github.com/ShinyEd/ShinyEd/tree/master/CLT_mean", target="_blank", "Ver código")),
    helpText(a(href="http://stat.duke.edu/~mc301/shiny/applets.html", target="_blank", "Echa un vistazo a otras aplicaciones")),
    helpText(a(href="https://www.coursera.org/course/statistics", target="_blank", "¿Quieres saber más de forma gratuita?")),
    br(),
    actionButton("sal", "Salir")
  ),
  
  
  
  mainPanel(
    plotOutput("pop.dist"),
    br(),
    plotOutput("sample.dist"),
    div(h3(textOutput("num.samples")), align = "center"),
    br(),
    plotOutput("sampling.dist"),
    br(),
    div(textOutput("sampling.descr"), align = "center"),
    br(),
    div(h5(textOutput("CLT.descr"), align = "center"))
  )
)