omaver<-function(a,b){
#(a<b), NA=infty
tulos<-F
if (is.na(a)) tulos<-F 
else{ if (is.na(b)) tulos<-T
      else{ if (a<b) tulos<-T}
}
return(tulos)
}
