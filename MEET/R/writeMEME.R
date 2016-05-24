writeMEME <-function(matriu, m){
  
  suma<-apply(matriu,2,function(X){sum(X=="-")})
  matriu<-matriu[,suma==0]
  primer<-as.numeric(colnames(matriu))

  PWM<-CalculPWM(matriu)

  final<-rep("-", times=nrow(PWM))
  PWM<-cbind(PWM, final)
  lletres<-c("A", "C", "G", "T", " ")
  names<-c("PO",1:ncol(matriu))
  PWM<-rbind(lletres, PWM)
  PWM.write<-data.frame(PWM, row.names=names)
  motif.head<-"ID "
  motif.head<-paste(motif.head, m)
  write(motif.head, file="motif.dat")
  write.table(PWM.write, "motif.dat", append=TRUE, quote=FALSE, col.names=FALSE, row.names=TRUE)
  write("XX", file="motif.dat", append=TRUE)

  return(primer[1])
 
}

