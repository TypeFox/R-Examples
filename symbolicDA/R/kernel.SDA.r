#***********************************************************************************************************************************************
#*  
#*  (C) 2010     Andrzej Dudek    Uniwersytet Ekonomiczny we Wrocławiu
#*  
#*  Jądrowa analiza dyskryminacyna dla danych symbolicznych
#*  Skrypt do książki:
#*  "Podejście symboliczne w analizie danych ekonomicznych" powstałej w ramach projektu badawczego habilitacyjnego N N111 105 234
#*  
#*  Kod poniższy może być modyfikowany, kopiowany i rozprowadzany na warunkach li-cencji GPL 2 (http://gnu.org.pl/text/licencja-gnu.html), 
#*  a w szczególności pod warunkiem umieszczenia w zmodyfikowanym pliku widocznej informacji o dokonanych zmianach, wraz z datą ich dokonania. 
#*  
#***********************************************************************************************************************************************
kernel.SDA<-function(sdt,formula,testSet,h,...)
{
  f<-.parseFormula(formula,sdt)
  clusters<-f$classes
  nrClusters<-max(clusters[testSet])
  p<-nk<-rep(0,nrClusters)
  pred<-rep(0,nrow(sdt$individuals)) 
  dist<-as.matrix(dist.SDA(sdt,variableSelection=f$variableSelection,...))
  for (i in testSet)
  {
    for(cl in 1:nrClusters)
    {
      p[cl]<-0
      nk[cl]<-sum(clusters==cl)
      for (j in (1:nrow(sdt$individuals))[-testSet])
      {
          if (dist[i,j]<=h)
          p[clusters[i]]<-p[clusters[i]]+1;
      }
    }
    #print(p)
    p=p/nk;
    pred[i]<-(1:nrClusters)[p==max(p)][1]
    #print(c(i,pred[i]))
  }
  pred[testSet]
}