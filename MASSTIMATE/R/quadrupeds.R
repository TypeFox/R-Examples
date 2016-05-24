quadrupeds <-
function(HC,FC,QE_MR.eq="raw",data=NULL) {
  HFC <- HC+FC
  #QE and MR estimate
  if(QE_MR.eq=="raw") {
    #linear equation
    QE.masstimate<-2.749*log10(HFC)-1.104
    QE<-round(10^QE.masstimate,1)
    QE.ppe.err<-QE*0.2563
    upper.QE<-round(QE+QE.ppe.err,1)
    lower.QE<-round(QE-QE.ppe.err,1)
    #quadratic equation
    log.hfc <- log10(HFC)
    qQE.masstimate <- -0.04856 * log.hfc^2 + 2.92291 * log.hfc - 1.24954
    qQE <- round(10^qQE.masstimate, 1)
    qQE.ppe.err <- qQE * 0.2537316
    upper.qQE <- round(qQE + qQE.ppe.err, 1)
    lower.qQE <- round(qQE - qQE.ppe.err, 1)
    #multiple regression
    MR.masstimate<-1.78*log10(HC)+0.939*log10(FC)-0.215
    MR<-round(10^MR.masstimate,1)
    MR.ppe.err<-MR*0.24932
    upper.MR<-round(MR+MR.ppe.err,1)
    lower.MR<-round(MR-MR.ppe.err,1)
  }
  if(QE_MR.eq=="phylocor") {
    #linear equations
    QE.masstimate<-2.754*log10(HFC)-1.097
    QE<-round(10^QE.masstimate,1)
    QE.ppe.err<-QE*0.2503
    upper.QE<-round(QE+QE.ppe.err,1)
    lower.QE<-round(QE-QE.ppe.err,1)
    #quadratic equation
    log.hfc <- log10(HFC)
    qQE.masstimate <- -0.0585856 * log.hfc^2 + 2.9629676 * log.hfc - 1.2646675
    qQE <- round(10^qQE.masstimate, 1)
    ppe.err <- qQE * 0.2469842
    upper.qQE <- round(qQE + ppe.err, 1)
    lower.qQE <- round(qQE - ppe.err, 1)
    #multiple regression
    MR.masstimate<-1.54*log10(HC)+1.195*log10(FC)-0.234
    MR<-round(10^MR.masstimate,1)
    MR.ppe.err<-MR*0.24624
    upper.MR<-round(MR+MR.ppe.err,1)
    lower.MR<-round(MR-MR.ppe.err,1)
  }
  #Anderseon et al 1985 estimate
  AHR1985<-round(0.078*(HFC)^2.73,2)
  #Mazzette et al 2004 estimate
  log.estimate<-2.955*log10(FC)-4.166
  MCF2004<-round((10^log.estimate)*1000,2)
  
  return(cbind(data,QE,lower.QE,upper.QE,qQE,lower.qQE,upper.qQE,MR,lower.MR,upper.MR,AHR1985,MCF2004))
}
