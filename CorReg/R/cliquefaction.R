# '  cliquefaction d'un graphe
cliquefaction<-function(Z){
  Z=Matrix(Z)
  p=ncol(Z)
  Zno=Z+t(Z)#le graphe devient non orient?
  Zcliq=Zno
  for(k in 1:p){
    Zno=Zno%*%(Z+t(Z))
    Zno[Zno>1]=1
    Zcliq=Zcliq+Zno
  } 
  Zcliq[Zcliq>1]=1
  return(Zcliq)
}