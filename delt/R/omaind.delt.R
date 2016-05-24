omaind.delt<-function(v){
#v on vektori, palautetaan indeksi jossa vektorin pienin arvo 
#
lkm<-length(v)
i<-1
while ((i<lkm) && (is.na(v[i]))) i<-i+1
if ((i==lkm) && (is.na(v[lkm]))) y<-1
 else
 if ((i==lkm) && (!is.na(v[lkm]))) y<-lkm
  else{
  apuu<-i
  valapu<-v[apuu]
  while (i<lkm){
    i<-i+1
    if ((!is.na(v[i])) && (v[i] < valapu)){
      apuu<-i
      valapu<-v[i]
    }
  }
y<-apuu
  }
return(y)
}
