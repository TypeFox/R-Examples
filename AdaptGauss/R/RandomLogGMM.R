RandomLogGMM=function(Means,SDs,Weights,IsLogDistribution,TotalNoPoints=1000){
  # GMM = RandomLogGMM(Means,SDs,Weights,IsLogDistribution,TotalNoPoints)
  # genereierung von ZufalsDaten, die einer Mischung von Vereilungen aus
  # Gauss & Log-Normalen folgt
  # INPUT
  # Means(1:L), SDs(1:L), Weights(1:L) die Paramter von der Verteilungen
  # IsLogDistribution(1:L) gibt an ob die einzelverteilung einer (generalisierten)Lognormaverteilung ist
  #             wenn IsLogDistribution(i)==0 dann Mix(i) = Weights(i) * N(Means(i),SDs(i)
  #             wenn IsLogDistribution(i)==1 dann Mix(i) = Weights(i) * LogNormal(Means(i),SDs(i)
  # Die Gesamtverteilung ergibst sich als Summe der GMM(1:L)
  # 
  # OPTIONAL
  # TotalNoPoints              Mix enthalet so viele Punkte am Schluss.
  #                              Default: TotalNrOfPoints ca 1000
  #
  # NOTA: die Log-Normalverteilung ist generaisiert d.h.: L(Means,SDs) = sign(Means)*lognrnd(abs(Means),SDs);
  #author: ?, RG
  #1.Eitor: MT 06/2015
  L= length(Means)
  AnzPoints = round(Weights*TotalNoPoints*1.1)
  # Verteilung der Mischungen  erzeugen
  Mix = c() # init
  for(d in 1:L){
    if(IsLogDistribution[d]==1){
      # Mixi = symlogrnd(Means[d],SDs[d],AnzPoints[d],1)
      # Mix(i) als LogNormal erzeugen
      #temp <- symlognSigmaMue(Means[d],SDs[d])
      variance<-log(SDs[d]*SDs[d]/(Means[d]*Means[d])+1)
      sig<-sqrt(variance)
      mu<-log(abs(Means[d]))-0.5*variance
      temp=list(variance=variance,sig=sig,mu=mu)
      mu <- temp$mu
      sig <- temp$sig
      n <- AnzPoints[d]*1
      Mixi <- matrix(sign(Means) * rlnorm(n, abs(mu), sig), AnzPoints[d], 1)
    }else{ # Mormalverteilung
      #Mixi = normrnd(Means[d],SDs[d],AnzPoints[d],1) # Mix(i) als Gauss erzeugen
      x <- 1*AnzPoints[d]
      Mixi <- matrix(rnorm(n=x, mean=Means[d], sd=SDs[d]), nrow=AnzPoints[d], ncol=1)     
    }
    Mix = c(Mix,Mixi)
  }  # for d
  # hier enthaelt Mix die der Vereilung entsprechende Punkte  
  # Dureinanderwuerfeln und auf gewuenschte Anzahl bringen
  Ind = sample(c(1:TotalNoPoints),TotalNoPoints)
  GMM=Mix[Ind]
  return(GMM)  
}