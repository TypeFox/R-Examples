output$contraBayesPlot <- renderPlot({
  achievable.se.sp(the.prev = input$prevalence)
  catn("contraBayesPlot lines ", PPVderived(), NPVderived())
  if(appName=="shinyElicit")
    abline(v=PPVderived(), h=NPVderived(), col="green", lwd=3)
})

output$parameterTable = renderTable({
  data = hoverData()
  if(is.null(data)) {
    sesp = sesp.from.NNT(NNTpos = input$NNTpos,
                         NNTneg=input$NNTneg,
                         prev=input$prevalence)
    data = data.frame(row.names = "",
      sensitivity=sesp["se"],
      specificity=sesp["sp"],
      PPV=1/input$NNTpos,
      NPV=1 - 1/input$NNTneg,
      NNTpos=input$NNTpos,
      NNTneg=input$NNTneg
    )
  }
  rValues$parameterTable = data
  xtable::xtable(digits=3, data)
})

hoverData = reactive({
  #catn("hoverData")
  # TODO: currently assumes  axes = "pv",
  # as opposed to NNTpos, NNTneg.
  ppv = input$contraBayesPlot_hover$x
  npv = input$contraBayesPlot_hover$y
  the.prev = input$prevalence
  nnt = NNT.from.pv(ppv = ppv, npv=npv)
  sesp = sesp.from.pv(ppv = ppv, npv=npv, prev=the.prev)
  result = try(silent = TRUE, (data.frame(
    sensitivity=sesp["se"], specificity=sesp["sp"],
    PPV=ppv,
    NPV=npv,
    NNTpos=nnt[1],
    NNTneg=nnt[2]
  )))
  if(class(result) == "try-error")
    return (NULL)
  rownames(result) = "+"
  #names(result) = "value"
  result
})
#catn("is.reactive(hoverData)= ", is.reactive(hoverData))

observe({
    catn("hoverinfo ", capture.output(input$contraBayesPlot_hover))
    if(!is.null(input$contraBayesPlot_hover))
       rValues$hoverinfo <- hoverData()
})
