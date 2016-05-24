intersec<-function(taso,endind,cur,uni){
#Makes from a set of rectangles "cur" a larger rectangle.
#For a given rectangle in cur, we make intersection with
#rectangles in uni, starting with rectangle after the point
#indicated in "endind". 
#Result is a set of k over "taso" rectangles, where k is the number
#of rectangles in "uni". 
#uni has the basic sets, we have in cur all the (taso-1)-fold
#intersections of uni, we want to form taso-fold intersections,
#in endind we have index of the last rectangle in (taso-1)-fold
#intersection: we have to form intersections with all the rest
#rectangles in uni. Thus result has the size: how many subsets of
#size taso, we can take from a set of size k. 
#
#taso is integer >=1
#endind is l-vector
#cur is l*(2*d)-matrix
#uni is k*(2*d)-matrix, huom oletetaan etta k>1!!!!!
#
#Return NA if there is no intersectio, otherwise
#list(set=tulos,endind=newendind)
#
k<-length(uni[,1])   #rows of uni
d2<-length(uni[1,])  #col of uni is the 2*d
tulrow<-choose(k,taso)
#tulrow<-gamma(k+1)/(gamma(k-taso+1)*gamma(taso+1))#k yli taso,gamma(k)=(k-1)!
newendind<-matrix(0,tulrow,1)
tulos<-matrix(0,tulrow,d2)
ind<-0        #indeksi to tulos and newendind
if (dim(t(cur))[1]==1) a<-1 else a<-length(cur[,1])  #rows of cur
for (i in 1:a){
  if (endind[i]<k){
    for (j in (endind[i]+1):k){
      if (a==1) apu<-leikkaa(cur,uni[j,])
        else apu<-leikkaa(cur[i,],uni[j,])
         #for (l in 1:d){
         #  tulos[ind,2*l-1]<-max(cur[i,2*l-1],uni[j,2*l-1])
         #  tulos[ind,2*l]<-min(cur[i,2*l],uni[j,2*l])
         #}      
      if (!is.na(apu)){  #if there is intersection, save the result
         ind<-ind+1
         tulos[ind,]<-apu
         newendind[ind]<-j
      }
    }
  }
}
if (ind==0) palauta<-NA
else{
  tulos<-tulos[1:ind,]
  newendind<-newendind[1:ind]
  palauta<-list(set=tulos,endind=newendind)
}
return(palauta) 
}







