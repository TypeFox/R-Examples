points.segmented<-function(x, term, interc=TRUE, link=TRUE, 
    rev.sgn=FALSE, transf=I, ...){
#--------------
        f.U<-function(nomiU, term=NULL){
        #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
        #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
            k<-length(nomiU)
            nomiUsenzaU<-strsplit(nomiU, "\\.")
            nomiU.ok<-vector(length=k)
            for(i in 1:k){
                nomi.i<-nomiUsenzaU[[i]][-1]
                if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".")
                nomiU.ok[i]<-nomi.i
                }
          if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
          return(nomiU.ok)
        }
#-------------
      if(missing(term)){
          if(length(x$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-x$nameUV$Z}
               }
      opz<-list(...)
      nameV<- x$nameUV$V[f.U(x$nameUV$V, term)]
      psii<- x$psi[nameV, "Est."]
      d<-data.frame(a=psii)
      names(d)<-term
      opz$y<-broken.line(x,d, se.fit=FALSE, interc=interc, link=link)[[1]]
      if(rev.sgn) psii<- -psii
      opz$x<- psii 
      if(is.null(opz$cex)) opz$cex<-1.5
      if(is.null(opz$lwd)) opz$lwd<-2
      opz$y<-do.call(transf, list(opz$y))
      do.call(points, opz)
      invisible(NULL)
      }

      
      

      