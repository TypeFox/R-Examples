bipeds <-
function(FC,cQE.eq="raw",cQE.cor=2,data=NULL) {    
  #cQE estimate
  if(cQE.eq=="raw") {
    #linear equation
    log.masstimate<-2.749*log10(FC*sqrt(cQE.cor))-1.104
    cQE<-round(10^log.masstimate,1)
    ppe.err<-cQE*0.2563
    upper.cQE<-round(cQE+ppe.err,1)
    lower.cQE<-round(cQE-ppe.err,1)
    #quadratic equation
    cor.FC <- log10(FC*sqrt(cQE.cor))
    qcQE.masstimate <- -0.04856 * cor.FC^2 + 2.92291 * cor.FC - 1.24954
    qcQE <- round(10^qcQE.masstimate, 1)
    ppe.err <- qcQE * 0.2537316
    upper.qcQE <- round(qcQE + ppe.err, 1)
    lower.qcQE <- round(qcQE - ppe.err, 1)
  }
  if(cQE.eq=="phylocor") {
    #linear equation
    log.masstimate<-2.754*log10(FC*sqrt(cQE.cor))-1.097
    cQE<-round(10^log.masstimate,1)
    ppe.err<-cQE*0.2503
    upper.cQE<-round(cQE+ppe.err,1)
    lower.cQE<-round(cQE-ppe.err,1)
    #quadratic equation
    cor.FC <- log10(FC*sqrt(cQE.cor))
    qcQE.masstimate <- -0.0585856 * cor.FC^2 + 2.9629676 * cor.FC - 1.2646675
    qcQE <- round(10^qcQE.masstimate, 1)
    ppe.err <- qcQE * 0.2469842
    upper.qcQE <- round(qcQE + ppe.err, 1)
    lower.qcQE <- round(qcQE - ppe.err, 1)
  }
  #Anderseon et al 1985 estimate
  AHR1985<-round(0.16*(FC)^2.73,2)
  #Christiansen and Farina 2004 estimate
  log.estimate<-2.738*log10(FC)-3.607
  CF2004<-round((10^log.estimate)*1000,2)
  #Campbell and Marcus 1992 estimate
  log.estimate<-2.411*log10(FC)-0.065
  CM1992<-round(10^log.estimate,2)
  return(cbind(data,cQE,lower.cQE,upper.cQE,qcQE,lower.qcQE,upper.qcQE,AHR1985,CF2004,CM1992))
}
