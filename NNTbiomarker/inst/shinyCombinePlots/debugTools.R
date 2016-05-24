##  For source'ing in server.R, and using in ui.R via uiOutput("debugTools")
## We begin with some convenient assignments and function.
if(interactive())
{
  thisSession <<- session
  reactiveValuesForDebugTools = reactiveValues()
  wasClicked =  function(button) {
    if(exists("input"))
      if(!is.null(button) ) {
        if(button > 0) {
          return(TRUE)
        }
      }
    return(FALSE)
  }
  #  assign("%&%",  function (a, b) paste(a, b, sep = ""))
  catn = function(...) cat(..., "\n")
  #  `%>%` = magrittr::`%>%`

  # Here begins the good stuff.
  output$evaluatedOutputR = renderUI({
    if(wasClicked(input$evalButtonR)) {
      cat("evaluatedOutputR\n")
      evalString = gsub('"', "'", isolate(input$evalStringR)) # replace all DQ with SQ.
      cat("R code evalString is ", evalString, "\n")
      result = capture.output(eval(parse(text=evalString)))
      print(str(result))
      result = paste(unlist(result), collapse = "\\n")
      result = gsub("<", "\\<", result)
      result = gsub(">", "\\>", result)
      print(str(result))
      cat("R code result is ", result, "\n")
      div(tags$script(
        paste0(  'alert("', result, '");')
      ))
      #       div(
      #         paste(collapse="<br>",
      #               gsub(" ", "&nbsp;", capture.output(eval(parse(text=evalString)))))
      #         ))
      ## You have to isolate input$evalStringR; otherwise each character typed calls this callback.
      ## The following might be useful later for up-arrowing through past expressions.
      #   if(is.null(reactiveValuesForDebugTools$evalStringHistory))
      #     reactiveValuesForDebugTools$evalStringHistory = character(0)
      #  reactiveValuesForDebugTools$evalStringHistory = c(reactiveValuesForDebugTools$evalStringHistory, evalString)
    }
  })

  outputPreambleJS <<- 'window.Shiny.shinyapp.$bindings.'
  # EXAMPLE:  window.Shiny.shinyapp.$bindings.selTxt.firstChild.nodeValue
  inputPreambleJS <<- 'window.Shiny.shinyapp.$inputValues.'
  wrapperToGetKeys <<- function(x) "Object.keys(" %&% x %&% ")"  # Currently not used.
  observerPreambleToggles = observe({
    input$prependInputPreambleToggle
    input$prependOutputPreambleToggle
    try({
      evalString = isolate(input$evalStringJS)
      if(wasClicked(input$prependInputPreambleToggle)) {
        if(substr(evalString, 1, nchar(inputPreambleJS)) != inputPreambleJS)
          evalString = paste0(inputPreambleJS, evalString)
      }
      else ## Remove inputPreambleJS
        evalString = gsub(inputPreambleJS, "", evalString, fixed=TRUE)

      if(wasClicked(input$prependOutputPreambleToggle)) {
        if(substr(evalString, 1, nchar(outputPreambleJS)) != outputPreambleJS)
          evalString = paste0(outputPreambleJS, evalString)
      }
      else ## Remove outputPreambleJS
        evalString = gsub(outputPreambleJS, "", evalString, fixed=TRUE)
      isolate( { reactiveValuesForDebugTools$evalStringJS = evalString } )
      catn("Responding to preamble toggles, evalString=", evalString)
      updateTextInput(thisSession, "evalStringJS", label="", value=reactiveValuesForDebugTools$evalStringJS)
      # You need to specify the label arg too. The default, NULL, doesn't cut it.
    })
  })

  #output$evaluatedOutputJS = renderText({
  #shinyalert("JS output is in a popup alert window, if there was no error. Otherwise nothing happens")
  # }
  # )

  output$JSevaluation = renderUI({
    if(wasClicked(input$evalButtonJS) ) {
      evalString = gsub('"', "'", isolate(input$evalStringJS)) # replace all DQ with SQ.
      div(list(tags$script(
        # 'alert(', '"HERE IS JS"', ')'     # THIS WORKS!
        # 'alert(eval(', '"1+2"', '))'       # THIS WORKS!
        paste0(
          'alert(eval("', evalString, '"))'       # THIS WORKS!
        )
      )))
    }
    # TRY THIS SOME TIME, to avoid creating an alert window for the JS output:
    #document.getElementById("demo").innerHTML = ... ;
  })

  output$shiny.trace.text = renderText({
    eval(options(shiny.trace=input$traceCheckbox), envir = .GlobalEnv);
    cat("shiny.trace: ", options("shiny.trace")[[1]], "\n")
    if( options("shiny.trace")[[1]] != input$traceCheckbox)
      cat('Error: options("shiny.trace")[[1]] should equal input$traceCheckbox', "\n");
    paste("trace=", input$traceCheckbox)
  })   #### OK this works now.

  observe({
    if(! is.null(input$errorOptionCheckbox)) {
      #     errorFunction = ifelse(input$errorOptionCheckbox,
      #                            (function () recover() ),
      #                            (function () .rs.breakOnError(userOnly = FALSE) ))
      #    print(errorFunction)
      if(input$errorOptionCheckbox==TRUE)
        eval(options(shiny.error = (function () recover() )),
             envir=.GlobalEnv)
      else
        eval(options(shiny.error = (function () .rs.breakOnError(userOnly = FALSE) )),
             envir=.GlobalEnv)
    }
  })

  output$debugTools = renderUI({
    div(style="background:darkGrey",
        singleton(tags$script(paste(
          "outputPreambleJS = '", outputPreambleJS, "';")))
        ,
        checkboxInput(inputId='debugToolsCheckbox', value=FALSE,
                      label=em(strong("Debugging aids"))),
        conditionalPanel(
          'input.debugToolsCheckbox',
          fluidRow(
            column(1,
                   tagAppendAttributes(style="color: blue",
                                       checkboxInput(inputId="traceCheckbox",
                                                     value=FALSE,
                                                     label=textOutput("shiny.trace.text")
                                       ))),
            column(1,
                   tagAppendAttributes(style="color: blue",
                                       checkboxInput(inputId="errorOptionCheckbox",
                                                     value=FALSE,
                                                     label="shiny.error: recover?")
                   )),
            column(4,  #class='well container-fluid',
                   # tagAppendAttributes(
                   "R code", HTML("&rarr;") ,
                   span(style="color:red; vertical-align='center'",
                        tags$textarea(id = "evalStringR", value="'getwd()'") ),
                   actionButton("evalButtonR",
                                HTML("<font color='red'> evaluate R</font>"),
                                #style='text-align:"right"; color:white',
                                uiOutput(outputId="evaluatedOutputR")
                   )
            )
            ,
            column(6,
                   column(8,
                          "JS code", HTML("&rarr;") ,
                          span(style="color:red; vertical-align='top'",
                               tags$textarea(id = "evalStringJS",
                                             value="'JS code here'") ),
                          actionButton("evalButtonJS",
                                       HTML("<font color='red'> evaluate JS</font>")
                                       #  style="display: flex; justify-content:flex-end;")
                                       # Cool! Too bad it doesn't work.
                          ) ),
                   column(4,
                          checkboxInput(inputId="prependOutputPreambleToggle",
                                        value=FALSE,
                                        label="prepend\nOutputPreamble")
                          , checkboxInput(inputId="prependInputPreambleToggle",
                                          value=FALSE,
                                          label="prepend\nInputPreamble")
                   ),
                   uiOutput(outputId='JSevaluation')
            )
          )
        )
    )
  })
}
