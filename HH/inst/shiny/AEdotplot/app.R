library(shiny)
library(HH)
data(AEdata)  ## this gets the app started.  Any data.frame can be entered once it is running.


shinyApp(
  ui=shinyUI(
    pageWithSidebar(
      headerPanel("AEdotplot"),
      sidebarPanel(
        ## textInput("data.frame", "data.frame", "AEdata", "150px"),
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
    localdata <- reactive(eval(parse(text=input$data.frame)))


    output$AEdata <- renderUI({
      dfnames <- unlist(objip(class="data.frame"))
      names(dfnames) <- NULL
      selectInput("data.frame", "data.frame on search list",
                  as.list(dfnames), "AEdata")
    })

    varClass <- reactive(
      lapply(get(input$data.frame), class)
    )

    varNumeric <- reactive(
      as.list(c("", names(varClass()[sapply(varClass(), function(x) ("numeric" %in% x) || ("integer" %in% x))])))
    )

    varFactor <- reactive(
      as.list(c("", names(varClass()[sapply(varClass(), function(x) ("factor" %in% x))])))
    )

    varFactorChar <- reactive(
      as.list(c("", names(varClass()[sapply(varClass(), function(x) ("factor" %in% x) || ("character" %in% x))])))
    )

    output$AE        <- renderUI({
      selectInput("AE",        "Adverse Event",       varFactor(), "AE"  )  })

    output$nAE       <- renderUI({
      selectInput("nAE",       "number AE observed",  varNumeric(), "nAE" )  })

    output$nTRT      <- renderUI({
      selectInput("nTRT",      "number on Treatment", varNumeric(), "nTRT")  })

    output$TRT       <- renderUI({
      selectInput("TRT",       "Treatment groups",           varFactor(), "TRT" )  })

    output$Condition <- renderUI({
      selectInput("Condition", "Condition",           as.list(c(None=" ", varFactorChar())), " "    )  })

    sortby <- reactive({
      whichcol <- c("Adverse Event"="PREF",
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
      validate(if (any(c(input$AE, input$nAE, input$nTRT, input$TRT)=="")) FALSE else NULL)

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
                     ", sortbyRelativeRisk=FALSE",
                     ","), collapse=" ")

      sub <- paste(c(" sub='",
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
                     input$data.frame,
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
