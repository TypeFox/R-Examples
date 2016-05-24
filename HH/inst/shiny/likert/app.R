library(shiny)
library(HH)
data(ProfChal)  ## this gets the app started.  Any data.frame can be entered once it is running.

shinyApp(
  ui=shinyUI(
    pageWithSidebar(
      headerPanel("likert"),
      sidebarPanel(
        ### textInput("data.frame", "data.frame", "ProfChal", "150px"),
        uiOutput("likertData"),

        ##textInput("Question",        "Question",        "Question",   "150px"),
        uiOutput("Question"),

        ## textInput("levels",       "formula for levels",       ".",  "150px"),
        uiOutput("levels"),

        textInput("main",      "Main Title",      "", "150px"),
        textInput("sub",       "Sub Title",       "",  "150px"),

        ## textInput("Condition", "Condition",  "",    "150px"),
        uiOutput("Condition"),

      ##  ## textInput("value",     "Response Count",  "",    "150px"),  ## long data
      ##  uiOutput("value"),


        checkboxInput("as.percent", "Show Percent", TRUE),
        checkboxInput("positive.order", "Positive Order", TRUE),
        checkboxInput("y.relation", "y.relation='free'", TRUE),
        textInput("layoutColumns", "layoutColumns", "", "150px"),
        textInput("layoutRows", "layoutRows", "", "150px"),
        textInput("ylab", "ylab: ' ' empty, '' name", " ", "150px"),
        textInput("px.height", "Pixel Height", "600", "150px")
      ),
      mainPanel(
        uiOutput("plotOutput")
      )
    )
  ),

  server=function(input, output) {
    localdata <- reactive(eval(parse(text=input$data.frame)))


     output$likertData <- renderUI({
      dfnames <- unlist(objip(class="data.frame"))
      names(dfnames) <- NULL
      selectInput("data.frame", "data.frame on search list",
                  as.list(dfnames), "ProfChal")
    })

    varClass <- reactive(
      lapply(get(input$data.frame), class)
    )

    varNumeric <- reactive(
      as.list(c("", names(varClass()[sapply(varClass(), function(x) ("numeric" %in% x) || ("integer" %in% x))])))
    )

    varFactorChar <- reactive(
      as.list(c("", names(varClass()[sapply(varClass(), function(x) ("factor" %in% x) || ("character" %in% x))])))
    )

    varFactorCharLevels <- reactive(
      as.list(c(".", names(varClass()[sapply(varClass(), function(x) ("factor" %in% x) || ("character" %in% x))])))
    )

    output$Question        <- renderUI({
      selectInput("Question",        "Question",       varFactorChar(), "Question"  )  })

    output$levels        <- renderUI({
      selectInput("levels",        "formula for levels",       varFactorCharLevels(), "."  )  })

    output$Condition <- renderUI({
      selectInput("Condition", "Condition",           as.list(c(None=" ", varFactorChar())), " "    )  })

    ## output$value <- renderUI({
    ##   selectInput("value", "value",       as.list(c(None=" ", varNumeric())), " "    )  })


    output$likert <- renderPlot({
      validate(if (any(c(input$Question, input$levels)=="")) FALSE else NULL)

      tmp <- paste(c("likert(",
                     input$Question,
                     "~",
                     input$levels,
                     if (nchar(input$Condition)!=0 && input$Condition!=" ")
                       "|",
                     if (nchar(input$Condition)!=0 && input$Condition!=" ")
                       input$Condition,
                     " , data=",
                     "localdata()",
                     if (input$as.percent) ", as.percent=TRUE",
                     if (input$positive.order) ", positive.order=TRUE",
                     if (input$y.relation) ", scales=list(y=list(relation='free'))",

                     if (nchar(input$Condition)!=0 && input$Condition!=" ") {
                       if (nchar(input$layoutColumns)!=0 && nchar(input$layoutRows)!=0)
                         paste0(", layout=c(", input$layoutColumns, ",", input$layoutRows, ")")
                       else
                         if (nchar(input$layoutColumns)==0 && nchar(input$layoutRows)==0)
                           paste0(", layout=c(", 1, ",", length(levels(localdata()[[input$Condition]])), ")")
                     },

             ##        if (nchar(input$value)!=0 && input$Condition!=" ") paste0(", value='", input$value, "'"),

                     if (nchar(input$ylab)!=0) paste(", ylab='", input$ylab, "'"),
                     if (nchar(input$main)!=0) paste(", main='", input$main, "'"),
                     ","), collapse=" ")



      sub <- if (nchar(input$sub)!=0)
               paste(" sub='",  input$sub, "')")
             else
               paste(c(" sub='",
                       input$Question,
                       "~",
                       input$levels,
                       if (nchar(input$Condition)!=0 && input$Condition!=" ")
                         "|",
                       if (nchar(input$Condition)!=0 && input$Condition!=" ")
                         input$Condition,
                       ", data=",
                       input$data.frame,
                     ##  if (nchar(input$value)!=0 && input$Condition!=" ") paste0(", value=\"", input$value, "\""),
                       "')"), collapse=" ")

      fulltext <- paste(tmp, sub)

      eval(parse(text=fulltext))
    })

    output$plotOutput <- renderUI(
      plotOutput("likert", height=input$px.height)
    )

  }
)
