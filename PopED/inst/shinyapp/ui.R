library(shiny)

# Define UI for miles per gallon application
shinyUI(fluidPage(
  
  # Application title
  titlePanel("PopED - Population Experimental Design"),
  
  sidebarLayout(
    sidebarPanel(
      h2("Model Definition"),
      selectInput("struct_model", "Structural Model:",
                  list(
                    "PK: 1-cpt, 1st order abs., single dose, CL param." = "ff.PK.1.comp.oral.sd.CL",
                    "PK: 1-cpt, 1st order abs., single dose, KE param." = "ff.PK.1.comp.oral.sd.KE",
                    "PK: 1-cpt, 1st order abs., multi. dose, CL param." = "ff.PK.1.comp.oral.md.CL",
                    "PK: 1-cpt, 1st order abs., multi. dose, KE param" = "ff.PK.1.comp.oral.md.KE"
                    ,"PKPD: 1-cpt, 1st order abs., multi. dose, CL param., direct effect IMAX" = "ff.PKPD.1.comp.oral.md.CL.imax"
                    ,"PKPD: 1-cpt, single dose, CL param., direct effect EMAX" = "ff.PKPD.1.comp.sd.CL.emax"
                    )),#,width='100%'),
                  
      #br(),
      
      selectInput("bsv_model", "Between Subject Variability Model:",
                  list(
                    "Exponential" = "exp",
                    "Proportional" = "prop",
                    "Additive" = "add"
                  )),
      checkboxInput("per_param", label = "Choose per parameter", value = FALSE),  
      #textInput("param"),
      #br(),
      
      selectInput("ruv_model", "Residual Unexplained Variability Model:",
                  list(
                    "Additive + Proportional" = "feps.add.prop",
                    "Proportional" = "feps.prop",
                    "Additive" = "feps.add"
                    )),
      
      h2("Design Definition"),
      conditionalPanel(
        condition = "input$struct_model == ff.PK.1.comp.oral.sd.CL",
        textInput("xt", "Sample times:", "0.5,1,2,6,24,36,72,120" )
      ),
      
      ##  create plot of model 
      ## plot_model_prediction(poped.db,IPRED=T,DV=T)
      
      #     checkboxInput("IPRED", "Show IPRED", FALSE),
      #     
      #     checkboxInput("DV", "Show DV", FALSE),
      #     checkboxInput("separate.groups", "Separate Groups", FALSE)
      
      ## identifiers
      br(),
      img(src = "poped_splash.png", height = 72, width = 72), 
      a(paste("PopED for R (", packageVersion("PopED"),")",sep=""), 
        href = "http://poped.sf.net"),
      h6("(c) 2014, Andrew C. Hooker, Pharmacometrics Research Group, Uppsala University, Sweden")
    ),
    
    mainPanel(
      #h3(textOutput("caption")),
      #     checkboxInput("smooth", "Smooth"),
      #     conditionalPanel(
      #       condition = "input.smooth == true",
      #       selectInput("smoothMethod", "Method",
      #                   list("lm", "glm", "gam", "loess", "rlm"))
      #     ),
      tabsetPanel(
        tabPanel("Plot model/design", 
                 plotOutput("modelPlot"),
                 #submitButton("Update View"),
                 checkboxInput("IPRED", "Show IPRED", FALSE),
                 checkboxInput("DV", "Show DV", FALSE),
                 checkboxInput("separate.groups", "Separate Groups", FALSE)), 
        tabPanel("Evaluate design", verbatimTextOutput("summary")),
        tabPanel("Optimize design", tableOutput("table")))
      #plotOutput("modelPlot")
    )
  )
))
