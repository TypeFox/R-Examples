lefrig2par<-function(et)
{
#from left,right representation to parent representation

d<-length(et$N)
len<-length(et$mean)

parent<-matrix(NA,len,1)
level<-matrix(0,len,1)
volume<-matrix(0,len,1)
center<-matrix(0,d,len)

parent[1]<-0              #root has no children
level[1]<-1
volume[1]<-1

pino<-matrix(0,len,1)
pinoind<-1
pino[1]<-1
while (pinoind>0){
   cur<-pino[pinoind]
   pinoind<-pinoind-1
    
   if (et$left[cur]!=0){
      parent[et$left[cur]]<-cur
      parent[et$right[cur]]<-cur
      level[et$left[cur]]<-level[cur]+1
      level[et$right[cur]]<-level[cur]+1
      volume[et$left[cur]]<-volume[cur]/2
      volume[et$right[cur]]<-volume[cur]/2
      center[,et$left[cur]]<-rep(1,d)          #left is furhtest from origo
      center[,et$right[cur]]<-rep(0,d)
   }
   
   while (et$left[cur]>0){
       #
       # laita oikea pinoon, 
       #
       oikea<-et$right[cur]
       pinoind<-pinoind+1
       pino[pinoind]<-oikea
       #
       # go to left
       #
       cur<-et$left[cur]
       #
       if (et$left[cur]!=0){
          parent[et$left[cur]]<-cur
          parent[et$right[cur]]<-cur
          level[et$left[cur]]<-level[cur]+1
          level[et$right[cur]]<-level[cur]+1
          volume[et$left[cur]]<-volume[cur]/2
          volume[et$right[cur]]<-volume[cur]/2
          center[,et$left[cur]]<-rep(1,d)       #left is furthest from origo
          center[,et$right[cur]]<-rep(0,d)
       }
   }
}

parent2<-matrix(0,len,1)
level2<-matrix(0,len,1)
volume2<-matrix(0,len,1)
center2<-matrix(0,d,len)
codeb<-matrix(0,d,len)
laskuri<-0
for (i in 1:len){
  if (!is.na(parent[i])){
     laskuri<-laskuri+1
     parent2[laskuri]<-parent[i]
     level2[laskuri]<-level[i]
     volume2[laskuri]<-volume[i]
     center2[,laskuri]<-center[,i]
     codeb[laskuri]<-i
  }
}
parent2<-parent2[1:laskuri]
level2<-level2[1:laskuri]
volume2<-volume2[1:laskuri]
center2<-center[,1:laskuri]
codeb<-codeb[1:laskuri] 

i<-2
while (i<=laskuri){
   cod<-parent2[i]
   j<-1
   while ((j<=laskuri) && (cod!=codeb[j])){
      j<-j+1
   }
   parent2[i]<-j
   i<-i+1
}

return(list(parent=parent2,level=level2,center=center2,volume=volume2))
}









