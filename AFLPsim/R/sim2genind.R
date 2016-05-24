sim2genind<-  function(x){
      x$S <- NULL
      x$SelMarkers <- NULL
      x<-x[!is.na(x)]
    if(class(x[[1]])=="list"){
      x<-unlist(x,recursive=F)
      x <- lapply(x, function(x)x[,-1 ])
    }
    N<-length(x)
    N2<-N
    if (names(x)[N]=="S")
      x$S<-NULL
    popnames<-names(x)
    n<-numeric(N2)
    for (i in 1:N2){ 
        n[i]<-nrow(x[[i]])
      }
    pop<-rep(popnames[1:N2],n[1:N2])
      
     X<-do.call("rbind",x)
     genindobj<-df2genind(X, pop=pop, ploidy=2, type="PA",  ncode=1)
      genindobj
}
