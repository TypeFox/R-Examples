omamindelt<-function(a,b){
#min(a,b), NA=infty, 
tulos<-F
if (is.na(a)) tulos<-b 
else{ if (is.na(b)) tulos<-a
      else{ if (a<=b) tulos<-a else tulos<-b}
}
return(tulos)
}
