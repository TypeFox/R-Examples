`TDT.pvalue` <-
function(pP,qP, pX,xX,qX, pQ,qQ, exact = FALSE)               #Spielman
  SMN.pvalue(pP+(2*pX+xX)+pQ, qP+(2*qX+xX)+qQ, exact = exact) #  (1993)
