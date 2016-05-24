draw.history<-function(obj,term,...){
#show.history() se c'e' stato boot restart potrebbe produrre un grafico 2x1 di "dev vs it" and "no.of distinct vs it"
#--
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
#--        
      if(missing(term)){
          if(length(obj$nameUV$Z)>1 ) {stop("please, specify `term'")}
               else {term<-obj$nameUV$Z}
               }
      range.ok<-obj$rangeZ[,term]
      #id.ok<-grep(paste("\\.",term,"$",sep=""),  rownames(obj$psi),value=FALSE)
      id.ok<- f.U(rownames(obj$psi), term)
      est.psi<-obj$psi[id.ok,2]
  if(length(obj$psi.history)==5) { #boot (non-autom)
        par(mfrow=c(1,2))
        plot(obj$psi.history$all.selected.ss, type="b", xlab="bootstrap replicates", 
            ylab="RSS  (selected values)", xaxt="n", pch=20)
        axis(1,at=1:length(obj$psi.history$all.selected.ss),cex.axis=.7)        
        #unicita' delle soluzioni
        if(is.vector(obj$psi.history$all.selected.psi)){
            psi.matr<-m<-matrix(obj$psi.history$all.selected.psi, ncol=1)
            } else {
              psi.matr<-m<-obj$psi.history$all.selected.psi[,id.ok,drop=FALSE]
              }
        
        for(i in 1:nrow(m)) m[i,]<-apply(psi.matr[1:i,,drop=FALSE],2,function(xx)length(unique(xx)))
        m<-t(t(m)+.1*(0:(ncol(m)-1)))
        matplot(1:nrow(m),m, pch=1:ncol(m), type="b", col=1:ncol(m), 
            ylab="no. of distinct solutions",xlab="bootstrap replicates", xaxt="n")
        axis(1,at=1:nrow(m),cex.axis=.7)        
        } else {
    if(all(diff(sapply(obj$psi.history, length))==0)){ #non-boot, non-autom
      A<-t(matrix(unlist(obj$psi.history),nrow=nrow(obj$psi),byrow=FALSE))
      colnames(A)<-rownames(obj$psi)
      matplot(1:nrow(A),A[,id.ok],type="o",pch=1:length(est.psi),col=1,
        xlab="iterations", ylab=paste("breakpoint ","(",term,")",sep=""),
        ylim=range.ok, xaxt="n",...)
      axis(1,at=1:nrow(A),cex.axis=.7)
      #if(rug) points(rep(1))
      abline(h=est.psi,lty=3)
      } else { #non-boot, Autom
          id.iter<-rep(1:length(obj$psi.history), times=sapply(obj$psi.history, length))
          psi.history<-unlist(obj$psi.history)
          nomi<-unlist(sapply(obj$psi.history, names))
          d<-data.frame(iter=id.iter, psi=psi.history, nomi=nomi)
          #associa i nomi delle componenti di $psi.history (che sono indici 1,2,..) con i nomi della variabile term
          ii<-unique(names(obj$psi.history[[length(obj$psi.history)]])[id.ok])
          if(length(ii)>1) stop("some error in the names?..")            
          with(d[d$nomi==ii,], plot(iter, psi,
                    xlab="iterations", ylab=paste("breakpoint ","(",term,")",sep=""),
                    xaxt="n",...))
                    axis(1,at=unique(d$iter),cex.axis=.7)
            #se vuoi proprio associare le stime tra le diverse iterazioni 
            #(per poi unire nel grafico i punti con le linee. Ovviamente alcune linee saranno interrotte)
#            for(i in 1:length(obj$psi.history)) {
#                a<-obj$psi.history[[i]]
                
#                for(j in 1:length(est.psi)){  
#                    psij<-est.psi[j]
                    #a<- ..names match
#                    r[i,j]<-a[which.min(abs(a-psij))]
#                    a<-setdiff(a, r[i,j])                    
#                    }
#                }        
        }
        }
     } #end_fn
      
