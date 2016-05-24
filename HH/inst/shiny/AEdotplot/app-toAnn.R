library(shiny)
library(HH)
data(AEdata)  ## this gets the app started.  Any dataset can be entered once it is running.


shinyApp(
  ui=shinyUI(
    pageWithSidebar(
      headerPanel("AEdotplot"),
      sidebarPanel(
        ## textInput("AEdata", "data.frame", "AEdata", "150px"),
        uiOutput("AEdata"),

        ## textInput("AE",        "Adverse Event",        "AE",   "150px"),
        uiOutput("AE"),

        ## textInput("nAE",       "number AE observed",       "nAE",  "150px"),
        uiOutput("nAE"),

        ## textInput("nTRT",      "number on Treatment",      "nTRT", "150px"),
        uiOutput("nTRT"),

        ## textInput("TRT",       "Treatment",       "TRT",  "150px"),
        uiOutput("TRT"),

        ## textInput("Condition", "Condition",  "",    "150px"),
        uiOutput("Condition"),

        selectInput("sortby", "Sort by",
                    c("Adverse Event",
                      "Percent A",
                      "Percent B",
                      "N A",
                      "N B",
                      "AE A",
                      "AE B",
                      "Relative Risk (RR)",
                      "Asymptotic Standard Error(log(RR))",
                      "CI lower",
                      "CI upper"),
                    selected="Relative Risk (RR)"),
        checkboxInput("AEtable", "Show AE Table", TRUE),
        textInput("px.height", "Pixel Height", "600", "150px")
      ),
      mainPanel(
        uiOutput("plotOutput")
      )
    )
  ),

  server=function(input, output) {
    localdata <- reactive(eval(parse(text=input$AEdata)))


output$AEdata <- renderUI({
  selectInput("AEdata", "data.frame in current directory",
              as.list(c(ls.df(.GlobalEnv))))

})

    varClass <- reactive(
      lapply(get(input$AEdata), class)
    )

    output$AE        <- renderUI({ AEnames <- names(varClass()[sapply(varClass(), function(x) ("factor" %in% x) || ("character" %in% x))])
      selectInput("AE",        "Adverse Event",       as.list(AEnames), "AE"  )  })

    output$nAE       <- renderUI({ nAEnames <- names(varClass()[sapply(varClass(), function(x) ("numeric" %in% x) || ("integer" %in% x))])
      selectInput("nAE",       "number AE observed",  as.list(nAEnames), "nAE" )  })

    output$nTRT      <- renderUI({ nTRTnames <- names(varClass()[sapply(varClass(), function(x) ("numeric" %in% x) || ("integer" %in% x))])
      selectInput("nTRT",      "number on Treatment", as.list(nTRTnames), "nTRT")  })

    output$TRT       <- renderUI({ TRTnames <- names(varClass()[sapply(varClass(), function(x) ("factor" %in% x) || ("character" %in% x))])
      selectInput("TRT",       "Treatment groups",           as.list(TRTnames), "TRT" )  })

    output$Condition <- renderUI({ ConditionNames <- names(varClass()[sapply(varClass(), function(x) ("factor" %in% x) || ("character" %in% x))])
      selectInput("Condition", "Condition",           as.list(c(None=" ", ConditionNames)), " "    )  })

    sortby <- reactive({
      whichcol <- c("Adverse Event"="AE",
                    "Percent A"="PCT A",
                    "Percent B"="PCT B",
                    "N A"="SN A",
                    "N B"="SN B",
                    "AE A"="SAE A",
                    "AE B"="SAE B",
                    "Relative Risk (RR)"="relrisk",
                    "Asymptotic Standard Error(log(RR))"="ase.logrelrisk",
                    "CI lower"="relriskCI.lower",
                    "CI upper"="relriskCI.upper")[input$sortby]
      ab <- strsplit(whichcol, " ")[[1]]
      if (length(ab) == 2) {
        if (ab[2] == "B") {
          whichrows <- 2:1
          whichcol <- ab[1]
        }
        else {
          whichrows <- 1:2
          whichcol <- ab[1]
        }}
      else
        whichrows <- 1

      list(sortbyVar=whichcol, sortbyVarBegin=whichrows[1])
    })

    output$AEdotplot <- renderPlot({
      tmp <- paste(c("AEdotplot(",
                     input$AE,
                     "~",
                     input$nAE,
                     "/",
                     input$nTRT,
                     if (nchar(input$Condition)!=0 && input$Condition!=" ")
                       "|",
                     if (nchar(input$Condition)!=0 && input$Condition!=" ")
                       input$Condition,
                     ", groups=",
                     input$TRT,
                     ", data=",
                     "localdata()",
                     ", sortbyVar=sortby()$sortbyVar",
                     ", sortbyVarBegin=sortby()$sortbyVarBegin",
                     ","), collapse=" ")

      sub <- paste(c(", sub='",
                     input$AE,
                     "~",
                     input$nAE,
                     "/",
                     input$nTRT,
                     if (nchar(input$Condition)!=0 && input$Condition!=" ")
                       "|",
                     if (nchar(input$Condition)!=0 && input$Condition!=" ")
                       input$Condition,
                     ", groups=",
                     input$TRT,
                     ", data=",
                     input$AEdata,
                     "')"), collapse=" ")

      fulltext <- paste(tmp, sub)

      latticePlot <- eval(parse(text=fulltext))
      if (input$AEtable)
        latticePlot
      else
        print(latticePlot, AEtable=FALSE)
    })

    output$plotOutput <- renderUI(
      plotOutput("AEdotplot", height=input$px.height)
    )

  }
)
