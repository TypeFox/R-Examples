calc_skewness <-
function(input, FRQ_val = NULL, HWE_val = NULL, cal_val = NULL, imp_val = NULL, ...) {
  if(!is.null(c(FRQ_val, HWE_val, cal_val, imp_val))) { effect <- input$EFFECT[HQ_filter(data = input, FRQ_val = FRQ_val, HWE_val = HWE_val, cal_val = cal_val, imp_val = imp_val, ...)]
  } else { if(is.vector(input)) { effect <- input } else { effect <- input$EFFECT } }
  if(any(is.na(effect))) { effect <- effect[!is.na(effect)] }
  return(if(length(effect) == 0L) { NA } else { sum( (effect - mean(effect))^3) / ((length(effect)-1) * sd(effect)^3 ) } )
}
