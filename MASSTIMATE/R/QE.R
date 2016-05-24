QE <-
function(HFC=NULL,HC,FC,equation="raw",quadratic=FALSE,data=NULL) {
  if(is.null(HFC)) HFC <- HC+FC
  if(equation=="raw") {
    log.QE<-2.749*log10(HFC)-1.104
    QE<-round(10^log.QE,1)
    ppe.err<-QE*0.2563
    upper.QE<-round(QE+ppe.err,1)
    lower.QE<-round(QE-ppe.err,1)
    if(quadratic) {
      log.hfc <- log10(HFC)
      log.qQE <- -0.04856 * log.hfc^2 + 2.92291 * log.hfc - 1.24954
      qQE <- round(10^log.qQE, 1)
      ppe.err <- qQE * 0.2537316
      upper.qQE <- round(qQE + ppe.err, 1)
      lower.qQE <- round(qQE - ppe.err, 1)
    }
  }
  if(equation=="phylocor") {
    log.QE<-2.754*log10(HFC)-1.097
    QE<-round(10^log.QE,1)
    ppe.err<-QE*0.2503
    upper.QE<-round(QE+ppe.err,1)
    lower.QE<-round(QE-ppe.err,1)
    if(quadratic) {
      log.hfc <- log10(HFC)
      log.qQE <- -0.0585856 * log.hfc^2 + 2.9629676 * log.hfc - 1.2646675
      qQE <- round(10^log.qQE, 1)
      ppe.err <- qQE * 0.2469842
      upper.qQE <- round(qQE + ppe.err, 1)
      lower.qQE <- round(qQE - ppe.err, 1)
    }
  }
  if(quadratic) res <- cbind(data,log.QE,QE,lower.QE,upper.QE,log.qQE,qQE,lower.qQE,upper.qQE)
  else res <- cbind(data,log.QE,QE,lower.QE,upper.QE)
  return(res)
}
