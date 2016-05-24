cQE <-
function(FC,equation="raw",cor=2,quadratic=FALSE,data=NULL) {
  if(equation=="raw") {
    log.cQE<-2.749*log10(FC*sqrt(cor))-1.104
    cQE<-round(10^log.cQE,1)
    ppe.err<-cQE*0.2563
    upper.cQE<-round(cQE+ppe.err,1)
    lower.cQE<-round(cQE-ppe.err,1)
    if(quadratic) {
      cor.FC <- log10(FC*sqrt(cor))
      log.qcQE <- -0.04856 * cor.FC^2 + 2.92291 * cor.FC - 1.24954
      qcQE <- round(10^log.qcQE, 1)
      ppe.err <- qcQE * 0.2537316
      upper.qcQE <- round(qcQE + ppe.err, 1)
      lower.qcQE <- round(qcQE - ppe.err, 1)
    }
  }
  if(equation=="phylocor") {
    log.cQE<-2.754*log10(FC*sqrt(cor))-1.097
    cQE<-round(10^log.cQE,1)
    ppe.err<-cQE*0.2503
    upper.cQE<-round(cQE+ppe.err,1)
    lower.cQE<-round(cQE-ppe.err,1)
    if(quadratic) {
      cor.FC <- log10(FC*sqrt(cor))
      log.qcQE <- -0.0585856 * cor.FC^2 + 2.9629676 * cor.FC - 1.2646675
      qcQE <- round(10^log.qcQE, 1)
      ppe.err <- qcQE * 0.2469842
      upper.qcQE <- round(qcQE + ppe.err, 1)
      lower.qcQE <- round(qcQE - ppe.err, 1)
    }
  }
  if(quadratic) res <- cbind(data,log.cQE,cQE,lower.cQE,upper.cQE,log.qcQE,qcQE,lower.qcQE,upper.qcQE)
  else res <- cbind(data,log.cQE,cQE,lower.cQE,upper.cQE)
  return(res)
}
