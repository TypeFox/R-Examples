#***********************************************************************************************************************************************
#*  
#*  (C) 2010     Andrzej Dudek    Uniwersytet Ekonomiczny we Wrocławiu
#*  
#*  Mapy Kohonena dla danych symbolicznych opisanych zmiennymi interwałowymi
#*  Skrypt do książki:
#*  "Podejście symboliczne w analizie danych ekonomicznych" powstałej w ramach projektu badawczego habilitacyjnego N N111 105 234
#*  
#*  Kod poniższy może być modyfikowany, kopiowany i rozprowadzany na warunkach li-cencji GPL 2 (http://gnu.org.pl/text/licencja-gnu.html), 
#*  a w szczególności pod warunkiem umieszczenia w zmodyfikowanym pliku widocznej informacji o dokonanych zmianach, wraz z datą ich dokonania. 
#*  
#***********************************************************************************************************************************************
kohonen.SDA<-function(data, rlen = 100, alpha = c(0.05, 0.01)){
  individualsNo<-dim(data)[1]
  prot<-data[sample(1:41,16),,]
  dane<-data
  clas<-rep(0,nrow(data))
  for(z in 1:rlen)   {
  print(z)
    for(i in 1:individualsNo){
      d<--1
      for(j in 1:dim(prot)[1]){
          dis<-2^(3-1)*(dist(rbind(dane[i,,1],prot[j,,1]))^2+dist(rbind(dane[i,,2],prot[j,,2]^2)))
          if(d==-1 || dis<d){
            d<-dis
            clas[i]<-j
          }
        
      for(j in 1:dim(prot)[1]){
        prot[j,,1]<-prot[j,,1]+(0.05-(z-1)*0.0005)*exp(-dist(rbind(dane[i,,1],prot[j,,1]))/2)
        prot[j,,2]<-prot[j,,2]+(0.05-(z-1)*0.0005)*exp(-dist(rbind(dane[i,,2],prot[j,,2]))/2)
      }
      resul<- list(clas=clas,prot=prot)
      resul
      }
    }
  }

  plot(NULL,xlim=c(0,4),ylim=c(0,4),axes=FALSE,xlab="",ylab="")
  for(i in 0:3){
  for(j in 0:3){
    segments(i,j,i,j+1)
    segments(i,j,i+1,j)
    segments(i,j+1,i+1,j+1)
    segments(i+1,j,i+1,j+1)
    t<-NULL
    for(k in 1:length(resul$clas)){
      if(resul$clas[k]==4*i+j+1){
        t<-c(t,k)
      }
    }
    if(!is.null(t)){
      text(i+0.5,j+0.5,paste(t,collapse=";"))
    }
  }
  }
  resul
}


