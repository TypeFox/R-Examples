library(shiny)

shinyUI(pageWithSidebar(
  
  headerPanel("Ejercicios de distribución"),
  
  sidebarPanel(
    #radio button or dropdown?

    selectInput(inputId = "dist_CalDis",
                label = "Distribución:",
                choices = c("Normal"      = "rnorm",
                            "Binomial"    = "rbinom",
                            "t"           = "rt",
                            "F"           = "rf",
                            "Chi-Cuadrado" = "rchisq"),
                selected = "rnorm"),
    

    br(),

    uiOutput("mean_CalDis"),
    uiOutput("sd_CalDis"),
    uiOutput("df1_CalDis"),
    uiOutput("df2_CalDis"),
    uiOutput("n_CalDis"),
    uiOutput("p_CalDis"),

    br(),
    br(),

    helpText("Modelo:"),
    h4(textOutput("model_CalDis"), align = "center"),
    br(),

    uiOutput("tail_CalDis"),
    uiOutput("lower_bound_CalDis"),
    uiOutput("upper_bound_CalDis"),
    

    uiOutput("a_CalDis"),
    uiOutput("b_CalDis"),
    
    withMathJax(),  # include the MathJax library
    
    conditionalPanel(condition = 'input.dist_CalDis === "rnorm"',
                     p("Para tipificar una variable \\(X\\) que sigue una distribución \\(\\ N(\\mu,\\sigma)\\)  se usa:
      $$Z=\\frac { X-\\mu  }{ \\sigma  } $$"),
                     p("Donde:"),
                     p("\\(\\mu\\) es la media"),
                     p("\\(\\sigma\\) desviación típica"),
                     tags$hr(),
                     helpText("Para saber más sobre la distribución normal leer ", 
                              a(href="http://es.wikipedia.org/wiki/Distribuci%C3%B3n_normal", target="_blank", "esto,"),
                              " Y recuerda visitar la página de", a(href="http://www.inegi.info/sofi", target="_blank", "Sofi."))
                    ),
    conditionalPanel(condition = 'input.dist_CalDis === "rbinom"',
                     p("Recuerda que la función de probabilidad de la distribución binomial es:
      $$p(X=k)=\\left( \\begin{matrix} n \\\\ k \\end{matrix} \\right) { p }^{ k }\\bullet { q }^{ n-k }$$"),
                     p("Donde:"),
                     p("\\(\\ n\\) es el número de pruebas."),
                     p("\\(\\ k\\) es el número de éxitos."),
                     p("\\(\\ p\\) es la probabilidad de éxito."),
                     p("\\(\\ q\\) es la probabilidad de fracaso."),
                     p("Y el coeficiente binomial: $$\\left( \\begin{matrix} n \\\\ k \\end{matrix} \\right) =\\frac { n! }{ k!(n-k)! }$$"),
                     tags$hr(),
                     helpText("Para saber más sobre la distribución binomial leer ", 
                              a(href="http://es.wikipedia.org/wiki/Distribuci%C3%B3n_binomial", target="_blank", "esto,"),
                              " Y recuerda visitar la página de", a(href="http://www.inegi.info/sofi", target="_blank", "Sofi."))
                     ),
    tags$hr(),
        p('Esta aplicación esta en desarrollo gracias a el útil código de ',a(href="https://github.com/mine-cetinkaya-rundel", target="_blank", "Mine Cetinkaya-Rundel"),' asi 
                       como de los paquetes ',a(href="http://cran.rstudio.com/web/packages/shiny/index.html", target="_blank", "shiny"),' y 
      ',a(href="http://cran.rstudio.com/web/packages/openintro/index.html", target="_blank", "openintro")
    ),
    helpText(HTML('<br/>Para dudas y sugerencias escribir a él <a href="mailto:jose.loera@inegi.org.mx?subject=Sofi" title="Enviar correo a Daniel">desarrollador</a>')
    ),
    
    actionButton("sal", "Salir")
    #tags$hr()
    
    ),
  
  mainPanel(
    #uiOutput('report'),
    #div(h4(textOutput("Peso_Est"), align = "center")),
    conditionalPanel(condition = 'input.Ejem_Dis === "Peso_Est"',
                     h3("Peso de estudiantes",align = "center"),
                     h4(textOutput("Doc_Peso"))
    ),
    conditionalPanel(condition = 'input.Ejem_Dis === "Tiro_Arc"',
                     h3("Tiro de arco",align = "center"),
                     h4(textOutput("Doc_Tiro"))
    ),
    conditionalPanel(condition = 'input.Ejem_Dis === "Temp_Est"',
                     h3("Temperatura",align = "center"),
                     h4(textOutput("Doc_Temp"))
    ),
    checkboxInput('Ayuda_visible', 'Mostrar ayuda (alfa)', FALSE),
    uiOutput('Ayuda'),
    #p(textOutput("status1"),style="font-weight=500; color: #000000;"),
    #h5(textOutput("status2"),style="font-weight=500; color: #00CC00;"),
    #h5(textOutput("status3"),style="font-weight=500; color: #FF0000;"),
    br(),
    #actionButton("submit","Enviar"),
    actionButton("newdat","Nuevos datos"),
    actionButton("newEje","Nuevo Ejemplo"),
    tags$hr(),
    
    column(3,
           numericInput("Res_Cuest",
                        "Y tu respuesta es:",
                        value = "0",
                        step=0.01), offset = 1
    ),
    #br(),
    column(8,
           actionButton("Resp_Ejem","Enviar")
    ),
    h5(textOutput("score")),
    p(textOutput("status1"),style="font-weight=500; color: #ddd;"),
    h5(textOutput("status2"),style="font-weight=500; color: #20f;"),
    h5(textOutput("status3"),style="font-weight=500; color: #f00000;"),
    tags$hr(),
    #conditionalPanel(condition = 'input.Ejem_Dis === "Peso_Est"',
    #                 #h3("Peso de estudiantes",align = "center"),
    #                 h4(textOutput("R_Doc_Peso"))
    #),
    #conditionalPanel(condition = 'input.Ejem_Dis === "Tiro_Arc"',
    #                 #h3("Tiro de arco",align = "center"),
    #                 h4(textOutput("R_Doc_Tiro"))
    #),
    #conditionalPanel(condition = 'input.Ejem_Dis === "Temp_Est"',
    #                 #h3("Temperatura",align = "center"),
    #                 h4(textOutput("R_Doc_Temp"))
    #),
    
    #br(),
    
    #h4(textOutput("Error_CalDis")),
    h4(textOutput("area_CalDis"), align = "center"),
    plotOutput("plot_CalDis"),
    #textOutput("area_CalDis")
    #helpText("Some math here, $$Y = \\beta_0 + \\beta_1 x + \\epsilon$$"),
    tags$hr(),
    radioButtons("Ejem_Dis", "Ejemplos: ",
                 c("Peso de estudiante" = "Peso_Est",
                   "Tiro de arco" = "Tiro_Arc",
                   "Temperatura" = "Temp_Est"), inline = TRUE)
    
    #uiOutput('report')
  )
))
