####  shinyElicit server
appName = "shinyElicit"
cat("Launching ", appName, "\n")

require("shiny")
require("xtable")
require("NNTbiomarker")
require("knitr")

shinyServerFunction =
  function(input, output, session) {
    thisSession <<- session

    #source("conveniences.R", local=TRUE)
    source("debugTools.R", local=TRUE)
    source("contraBayesPlot.R", local=TRUE)
    source("report.R", local=TRUE)


    observe({
      if(!is.null(input$contraBayesPlot_click)) {
        ppv = input$contraBayesPlot_click$x
        npv = input$contraBayesPlot_click$y
        nnts = pv.to.NNT(ppv = ppv, npv = npv)
        catn("nnts observed: ", nnts[[1]], nnts[[2]])
        if(all(!is.nan(nnts))) {
          updateNumericInput(session, inputId="NNTneg",
                             value=round(nnts[[2]]))
          updateNumericInput(session, inputId="NNTpos",
                             value=round(nnts[[1]]))
        }
       }
    })
    rValues = reactiveValues(
      stepsTable = stepsTableInitial )
#     stepsTable = data.frame(stringsAsFactors = FALSE,
#                             `Done?` = character(0),
#                             `Stepping stone` = character(0), Question = character(0)
#     )
#

    output$steps = renderTable({
      #catn("Calling renderTable on stepsTable");
      rValues$stepsTable
    })

    obs = function(number) {
      assign("obs" %&% number,
             observe({
               newValue <- input[["stepStatus" %&% number]]
               if( ! is.null(newValue)) {
                 catn("Toggling stepStatus" %&% number %&% " = " %&% newValue)
                 isolate({ # necessary, or else crash!
                   rValues$stepsTable[number, "Completed?"] = newValue
                   catn("New value in stepsTable: ",
                        rValues$stepsTable[number, "Completed?"] )
                 })
               }
             }
             )
      )
      #     cat("obs" %&% number %&% " is where? ",
      #         find("obs" %&% number))
      #     print(getAnywhere("obs" %&% number))
    }
    obsList = lapply(1:nrow(stepsTableInitial), obs)
    # find(obs1) # not found!
    #print(identical(obs1, obs2))
    #print(identical(obs1, obsList[[1]]))
    #print(identical(obsList[[1]], obsList[[7]]))
    ###
    ### ### must assign, or else only the last will take.

    observe({
      if(input$who == "" | input$Option_Treat == "" | input$Option_Wait == "")
        try(disableActionButton("stepStatus1", session))
    })
    observe({
      try(
        if(!all(sapply(1:nrow(stepsTableInitial),
                       function(n) "Done"==
                       input[["stepStatus" %&% n]])))
          disableActionButton("reportButton", session)
        else
          enableActionButton("reportButton", session)
      )
    })

    source("plotDiscomfort.R", local=TRUE)

    NNTgap = 1
#     observe( {
#       updateNumericInput(session, inputId="NNTpos",
#                          value=min(input$NNTpos,
#                                    isolate(input$NNTlower-NNTgap)),
#                          max = isolate(input$NNTlower-NNTgap))
#     })
#     observe( {
#       updateNumericInput(session, inputId="NNTlower",
#                          value=max(input$NNTlower,
#                                    isolate(input$NNTpos+NNTgap)),
#                          min = isolate(input$NNTpos+NNTgap))
#     })
#     observe( {
#       updateNumericInput(session, inputId="NNTneg",
#                          value=max(input$NNTneg,
#                                    isolate(input$NNTupper+NNTgap)),
#                          min = isolate(input$NNTupper+NNTgap))
#     })
#     observe( {
#       updateNumericInput(session, inputId="NNTupper",
#                          value=min(input$NNTupper,
#                                    isolate(input$NNTneg-NNTgap)),
#                          max = isolate(input$NNTneg-NNTgap))
#     })
    output$plotDiscomfort = renderPlot({
      plotDiscomfort(drawPosNeg=FALSE,
                     NNTpos = input$NNTpos,
                     NNTneg = input$NNTneg,
                     NNTlower = input$NNTlower,
                     NNTupper = input$NNTupper)
    }
    #, height=280
    )

    PPVderived = reactive({1/input$NNTpos})
    output$PPVderived = renderText({PPVderived()})
    NPVderived = reactive({1 - 1/input$NNTneg})
    output$NPVderived = renderText({NPVderived() })

    ## Prospective design:

    observe({
      cat("PROSPECTIVE INTERVALS (entering)\n")
      if( ! identical(input$NpatientsProspective, numeric(0))) {
        sampleSize = print(input$NpatientsProspective)
        Npositives = print((input$percentPositive/100 * sampleSize))
        Nnegatives = print(sampleSize - Npositives)
        NtruePositives = print( PPVderived() * Npositives)
        NtrueNegatives = print( NPVderived() * Nnegatives)
        prev = input$prevalence
        ProspectiveIntervals = NNTintervalsProspective(
            Npositives = round(Npositives),
            Nnegatives = round(Nnegatives),
            NtruePositives = round(NtruePositives),
            NtrueNegatives = min(Nnegatives - 1, round(NtrueNegatives)),
            prev = prev
          )
        print(str(rValues$ProspectiveIntervals))
        rValues$ProspectiveIntervals = ProspectiveIntervals
      }
      try({
        rValues$PPV_ProspectiveInterval = rValues$ProspectiveIntervals[
          c("lower boundary", "upper boundary"), "PPV"]
        rValues$NPV_ProspectiveInterval = rValues$ProspectiveIntervals[
          c("lower boundary", "upper boundary"), "NPV"]
        rValues$NNTpos_ProspectiveInterval = rValues$ProspectiveIntervals[
          c("lower boundary", "upper boundary"), "NNTpos"]
        rValues$NNTneg_ProspectiveInterval = rValues$ProspectiveIntervals[
          c("lower boundary", "upper boundary"), "NNTneg"]
      })
    })

#    output$PPVconfidenceInterval = renderText({rValues$})

###### Retrospective design ####

    observe({
      try({
        sesp = sesp.from.NNT(
          NNTpos = input$NNTpos,
          NNTneg = input$NNTneg,
          prev = input$prevalence
        )
        rValues$sensitivity = sesp["se"]
        rValues$specificity = sesp["sp"]
        rValues$sensitivityPercent = round(100 * rValues$sensitivity)
        rValues$specificityPercent = round(100 * rValues$specificity)
        rValues$NposCases = round(input$samplesizeCases* rValues$sensitivity)
        rValues$NposControls = round(input$samplesizeControls * (1-rValues$specificity))

        rValues$RetrospectiveIntervals = NNTintervalsRetrospective(
          Ncases = input$samplesizeCases,
          Ncontrols = input$samplesizeControls,
          NposCases = rValues$NposCases,
          NposControls = rValues$NposControls,
          prev = input$prevalence)
        rValues$Se_RetrospectiveInterval = rValues$RetrospectiveIntervals[
          c("lower boundary", "upper boundary"), "sensitivity"]
        rValues$Sp_RetrospectiveInterval = rValues$RetrospectiveIntervals[
          c("lower boundary", "upper boundary"), "specificity"]
        rValues$NNTpos_RetrospectiveInterval = rValues$RetrospectiveIntervals[
          c("lower boundary", "upper boundary"), "NNTpos"]
        rValues$NNTneg_RetrospectiveInterval = rValues$RetrospectiveIntervals[
          c("lower boundary", "upper boundary"), "NNTneg"]
      })
    })

    output$plotNNTgoals = renderPlot({
      plotDiscomfort(drawPosNeg=TRUE,
                     NNTlower = input$NNTlower,
                     NNTupper = input$NNTupper,
                     NNTpos = input$NNTpos,
                     NNTneg = input$NNTneg)
    }
    #, height=280
    )


    wasClicked =  function(button) {
      if(exists("input"))
        if(!is.null(button) ) {
          if(button > 0) {
            return(TRUE)
          }
        }
      return(FALSE)
    }
  }
#debug(shinyServerFunction)
shinyServer(func=shinyServerFunction)


