laterAFCM <- function (data, scannf=FALSE, nf=2, saveDatadisj=FALSE, fileDatadisj="Datadisj.csv"
                       , saveSumcolDatadisj=FALSE, fileSumcolDatadisj="SumcolDatadisj.csv"
                       , saveDataburt=FALSE, fileDataburt="Databurt.csv"
                       , saveContributions=FALSE, fileContributions="Contributions.csv")
{
  Data<-as.data.frame.matrix(data)
  for (i in 1:ncol(Data)) {
      Data[[i]]<-as.factor(Data[[i]])
  }

  # Disjunctive table
  Datadisj <- acm.disjonctif(Data)
  if (saveDatadisj == "csv") {write.csv(Datadisj, file = fileDatadisj)} else {}
  if (saveDatadisj == "csv2") {write.csv2(Datadisj, file = fileDatadisj)} else {}
  
  SumcolDatadisj<-colSums(Datadisj)
  if (saveSumcolDatadisj == "csv") {write.csv(SumcolDatadisj, file = fileSumcolDatadisj)} else {}
  if (saveSumcolDatadisj == "csv2") {write.csv2(SumcolDatadisj, file = fileSumcolDatadisj)} else {}
  
  # Burt table
  Databurt <- acm.burt(Data,Data)
  if (saveDataburt == "csv") {write.csv(Databurt, file = fileDataburt)} else {}
  if (saveDataburt == "csv2") {write.csv2(Databurt, file = fileDataburt)} else {}
  
  # Inertia
  if (scannf == "TRUE") {dev.new()} else {}
  acmData <- dudi.acm(Data, scannf=scannf, nf=nf)
  contri <- inertia.dudi(acmData,row.inertia=TRUE,col.inertia=TRUE)
  Contributions<-cbind(contri$col.abs, contri$col.rel,contri$col.cum)
  if (saveContributions == "csv") {write.csv(Contributions, file = fileContributions)} else {}
  if (saveContributions == "csv2") {write.csv2(Contributions, file = fileContributions)} else {}
  results<-c()
  results$Datadisj<-Datadisj
  results$SumcolDatadisj<-SumcolDatadisj
  results$Databurt<-Databurt
  results$Contributions<-Contributions
  results
}

