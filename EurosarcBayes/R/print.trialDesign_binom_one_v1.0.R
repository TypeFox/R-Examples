

.print.trialDesign_binom_one=function(x, ...){

  cat("Trial design and properties for single arm, single endpoint trial designs\n\n")

  m=data.frame(reviews=x@reviews)
  if(sum(is.finite(x@success))>0){m$success=x@success}
  if(sum(is.finite(x@failure))>0){m$failure=x@failure}

  print(m,row.names = FALSE)

  cat("\np0     :",x@p0)
  cat("\np1     :",x@p1)
  cat("\nAlpha  :",x@alpha)
  cat("\nPower  :",x@power)
  cat("\nExp(p0):",x@exp.p0)
  cat("\nExp(p1):",x@exp.p1)
  cat("\nEta    :",x@eta)
  cat("\nZeta   :",x@zeta,"\n")

}

setMethod(f="print",signature=c(x="trialDesign_binom_one"),definition=.print.trialDesign_binom_one)

.show.trialDesign_binom_one=function(object){

  cat("Trial design and properties for single arm, single endpoint trial designs\n\n")

  m=data.frame(reviews=object@reviews)
  if(sum(is.finite(object@success))>0){m$success=object@success}
  if(sum(is.finite(object@failure))>0){m$failure=object@failure}

  print(m,row.names = FALSE)

  cat("\np0     :",object@p0)
  cat("\np1     :",object@p1)
  cat("\nAlpha  :",object@alpha)
  cat("\nPower  :",object@power)
  cat("\nExp(p0):",object@exp.p0)
  cat("\nExp(p1):",object@exp.p1)
  cat("\nEta    :",object@eta)
  cat("\nZeta   :",object@zeta,"\n")

}

setMethod(f="show",signature="trialDesign_binom_one",definition=.show.trialDesign_binom_one)
