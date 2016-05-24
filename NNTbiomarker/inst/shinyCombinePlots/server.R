appName = "shinyCombinePlots"
cat("Launching ", appName, "\n")

shinyServerFunction =
  function(input, output, session) {
    thisSession <<- session
    rValues = reactiveValues()

    #source("conveniences.R", local=TRUE)
    source("debugTools.R", local=TRUE)
    source("contraBayesPlot.R", local=TRUE)


    observe({
      updateNumericInput(session, inputId="NtruePositives",
                         value=round(input$Npositives / input$NNTpos))
    })
    observe({
      updateNumericInput(session, inputId="NtrueNegatives",
                         value=round(input$Nnegatives * (1 - 1/input$NNTneg)))
    })

    output$intervalsProspective = renderTable({
      catn("In output$intervalsProspective")
      data = NNTintervalsProspective(
        Npositives = input$Npositives,
        Nnegatives = input$Nnegatives,
        NtruePositives = input$NtruePositives,
        NtrueNegatives = input$NtrueNegatives,
        prev = input$prevalence
      )
      print(data)
      xtable::xtable(digits=3, data)
    })

    observe({
      sensitivity = rValues$parameterTable[1, "sensitivity"]
      updateNumericInput(session, inputId="NposCases",
                         value=round(input$Ncases * sensitivity))
    })
    observe({
      specificity = rValues$parameterTable[1, "specificity"]
      updateNumericInput(session, inputId="NnegControls",
                         value=round(input$Ncontrols * specificity))
    })

    output$intervalsRetrospective = renderTable({
      catn("In output$intervalsRetrospective")
      data = NNTintervalsRetrospective(
        Ncases = input$Ncases,
        Ncontrols = input$Ncontrols,
        NposCases = input$NposCases,
        NposControls = input$Ncontrols - input$NnegControls,
        prev = input$prevalence
      )
      print(data)
      xtable::xtable(digits=3, data)
    })


    observe({
      if(!is.null(input$contraBayesPlot_click)) {
        ppv = input$contraBayesPlot_click$x
        npv = input$contraBayesPlot_click$y
        nnts = pv.to.NNT(ppv = ppv, npv = npv)
        catn("nnts observed: ", nnts[[1]], nnts[[2]])
        if(all(!is.nan(nnts))) {
          updateNumericInput(session, inputId="NNTneg",
                             value=(nnts[[2]]))
          updateNumericInput(session, inputId="NNTpos",
                             value=(nnts[[1]]))
        }
       }
    })

    source("plotDiscomfort.R", local=TRUE)

    output$plotDiscomfort = renderPlot({
      plotDiscomfort(drawPosNeg=FALSE,
                     NNTlower = input$NNTlower,
                     NNTupper = input$NNTupper)
    }
    #, height=280
    )

    PPVderived = reactive({1/input$NNTpos})
    output$PPVderived = renderText({PPVderived()})
    NPVderived = reactive({1 - 1/input$NNTneg})
    output$NPVderived = renderText({NPVderived() })

    output$plotNNTgoals = renderPlot({
      plotDiscomfort(drawPosNeg=TRUE,
                     NNTlower = input$NNTlower,
                     NNTupper = input$NNTupper,
                     NNTpos = input$NNTpos,
                     NNTneg = input$NNTneg)
    }
    #, height=280
    )
  }

#debug(shinyServerFunction)
shinyServer(func=shinyServerFunction)


