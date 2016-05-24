cutT<-function(T, step, mult=2, min=mult, max=+Inf){
  # cut T in Nstep regular sub intervals such that T/Nstep ~step
  # mult : Nstep = mult * i

  Nstep<- as.integer(T %/% step)
  Nstep<- as.integer(((Nstep %/% mult) + 1) * mult)
  Nstep<-ifelse(Nstep<min, min, Nstep)
  Nstep<-ifelse(Nstep>max, max, Nstep)
  data.frame(NstepT=Nstep, stepT=T/Nstep)
}


cutTfromto<-function(fromT, toT, step, mult=2, min=mult, max=+Inf){
  # cut T in Nstep regular sub intervals such that T/Nstep ~step
  # mult : Nstep est une multiple de mult

  Nstep<- as.integer((toT - fromT)  %/% step)
  Nstep<- as.integer(((Nstep %/% mult) + 1) * mult)
  Nstep<-ifelse(Nstep<min, min, Nstep)
  Nstep<-ifelse(Nstep>max, max, Nstep)
  data.frame(NstepT=Nstep, stepT=(toT - fromT) /Nstep)
}




