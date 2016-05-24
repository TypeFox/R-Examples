autoFillObserver = observe({
  cat("==> autoFillObserver\n")
  if(wasClicked(input$autoFill) ) {
    updateTextInput(session, "biomarkerReportTitle",
                    value="DEMO report for Statistical Considerations")
    updateTextInput(session, "ClinicalScenario",
                    value="Prognosis of cutaneous T cell
                    lymphoma (CTCL). In early stages of CTCL, patients (Stages IA-IIA)
                    usually do well and have slowly progressive disease, which does not
                    require aggressive therapy associated with substantial side effects.
                    However, about 15% of these patients have unexpected progressive course
                    and rapid demise.")
    updateTextInput(session, "who",
                    value="CTCL patients in Stages IA-IIA")
    updateTextInput(session, "Option_Treat",
                    value="Immediate aggressive treatment")
    updateTextInput(session, "Option_Wait",
                    value="Standard care")
    updateTextInput(session, "SpecificBenefit",
                    value="a biomarker progression risk model
                    that is able to classify patients into high and low risk groups will
                    enable personalized and more aggressive therapy for the patients at
                    highest risk for progression.")
    updateTextInput(session, "BestToTreatDescription", value = "with rapid progression")
    updateNumericInput(session, "prevalence", value = 0.15)
    updateNumericInput(session, "NNTlower", value = 3)
    updateNumericInput(session, "NNTupper", value = 50)
    updateNumericInput(session, "NNTpos", value = 2)
    updateNumericInput(session, "NNTneg", value = 30)
    updateNumericInput(session, "samplesizeCases", value = 22)
    updateNumericInput(session, "samplesizeControls", value = 40)
    for(stepNum in 1:env_sectionHeader$number)
      updateRadioButtons(session, "stepStatus" %&% stepNum, selected="Done")
  }
})
assembleReportObserver = observe({
  cat("==> assembleReportObserver\n")

  ### Only react when the reportButton is clicked.
  if(wasClicked(input$reportButton)) {
    ClinicalScenario = input$ClinicalScenario
    cat("input$options = ", capture.output(input$options), '\n')
    Option_Treat <<- input$Option_Treat
    Option_Wait <<- input$Option_Wait
    NNTlower <<- input$NNTlower
    NNTupper <<- input$NNTupper
    NNTpos <<- input$NNTpos
    NNTneg <<- input$NNTneg
    Benefit <<- input$benefit

    cat('getOption("markdown.HTML.options")',
        capture.output(getOption("markdown.HTML.options")), '\n')
    cat('getOption("markdown.extensions")',
        capture.output(getOption("markdown.extensions")), '\n')
    knit2html("www/Steps.Rmd", output = "www/Steps.html")
    # browseURL("Steps.html") ### Fails at shinyapp.io.
    #window.open("Steps.html"); # JS works if file is in www folder.
    print(getwd())
    print(dir())
  }
})
