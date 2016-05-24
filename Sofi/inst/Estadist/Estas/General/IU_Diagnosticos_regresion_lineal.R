pageWithSidebar(
  
  #  Application title
  headerPanel("Diagnósticos para la regresión lineal simple"),
  
  sidebarPanel(
    
    radioButtons("type", "Seleccione una tendencia:",
                 list("Línea ascendente" = "linear.up",
                      "Línea descendente" = "linear.down",
                      "Curvado hacia arriba" = "curved.up",
                      "Curvado hacia abajo" = "curved.down",
                      "En forma de abanico" = "fan.shaped")),
    br(),
    
    checkboxInput("show.resid", "Mostrar residuales", FALSE),
    
    br(),
    helpText("En esta aplicación que utilizan mínimos cuadrados ordinarios (MCO) para adaptarse a los datos seleccionados con una línea de regresión. 
             Está diseñada para ayudarle a practicar, evaluando si el modelo lineal es un ajuste apropiado a los datos. 
             Las tres parcelas de diagnóstico en la mitad inferior de la página son producidas para ayudarle a identificar patrones indeseables 
             en los residuos que puedan derivarse de las tendencias no lineales en los datos."),
    br(),
    helpText("Para saber más sobre regresión lineal leer ", 
             a(href="http://es.wikipedia.org/wiki/Regresi%C3%B3n_lineal", target="_blank", "esto,"),
             " Y recuerda visitar la página de", a(href="http://www.inegi.info/sofi", target="_blank", "Sofi.")),
    br(),
    helpText(a(href="https://github.com/ShinyEd/ShinyEd/tree/master/CLT_mean", target="_blank", "Ver código")),
    helpText(a(href="http://stat.duke.edu/~mc301/shiny/applets.html", target="_blank", "Echa un vistazo a otras aplicaciones")),
    helpText(a(href="https://www.coursera.org/course/statistics", target="_blank", "¿Quieres saber más de forma gratuita?")),
    br(),
    actionButton("sal2", "Salir")
    ),
  
  
  
  # Show the main display
  mainPanel(
    plotOutput("scatter"),
    br(),
    br(),
    plotOutput("residuals")
  )
  )