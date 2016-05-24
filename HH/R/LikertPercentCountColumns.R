LikertPercentCountColumns <-
  function(x, data,
           px=list( ## defaults designed for long QuestionName values
             LL=c(.00,  .50), ## and 7in x 7in window
             LP=c(.50,  .70),
             ML=c(.50,  .51),  ## arbitrary, visually center the labels and legend
             RP=c(.71,  .87),
             RL=c(.87, 1.00)),
           ...,
           QuestionName="Question",
           as.percent="Capture and then ignore this argument",
           positive.order=FALSE) {

    percentPlot <- likert(x, data,
                          as.percent="noRightAxis", ...,
                          positive.order=positive.order)
    ## percentPlot

    data2 <- data
    if (positive.order) {
      pct.order <- percentPlot$y.limits
      if (is.list(pct.order))
        pct.order <- percentPlot$y.limits[[1]]
      data2[[QuestionName]] <- factor(data2[[QuestionName]],
                                      levels=rev(pct.order))
    }

    countPlot <- likert(x, data2,
                        as.percent=FALSE, ...,
                        positive.order=FALSE,
                        rightAxis=TRUE)

    as.TwoTrellisColumns5(percentPlot, countPlot, px=px)
}
