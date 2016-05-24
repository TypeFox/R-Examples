pageWithSidebar(
  
  headerPanel("Calculadora de Distribución"),
  
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
    div(textOutput("model_CalDis"),style="text-indent:20px;font-size:125%;"),
    br(),
    
    uiOutput("tail_CalDis"),
    uiOutput("lower_bound_CalDis"),
    uiOutput("upper_bound_CalDis"),
    
    
    uiOutput("a_CalDis"),
    uiOutput("b_CalDis"),
    
    br(),
    
    helpText(a(href="https://duke.qualtrics.com/SE/?SID=SV_3L8WjmwQo32cVk9", target="_blank", "Rate this app!")),
    helpText(a(href="https://github.com/ShinyEd/ShinyEd/tree/master/dist_calc", target="_blank", "View code")),
    helpText(a(href="http://stat.duke.edu/~mc301/shiny/applets.html", target="_blank", "Check out other apps")),
    helpText(a(href="https://www.coursera.org/course/statistics", target="_blank", "Want to learn more for free?"))),
  
  
  
  mainPanel(
    plotOutput("plot_CalDis"),
    div(textOutput("area_CalDis"), align = "center", style="font-size:150%;")
  )
)