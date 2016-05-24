`confint.segmented` <-
function(object, parm, level=0.95, rev.sgn=FALSE, var.diff=FALSE, digits=max(3, getOption("digits") - 3), ...){
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
#        if(!"segmented"%in%class(object)) stop("A segmented model is needed")
        if(var.diff && length(object$nameUV$Z)>1) {
            var.diff<-FALSE
            warning("var.diff set to FALSE with multiple segmented variables", call.=FALSE)
            }
        #nomi delle variabili segmented:
        if(missing(parm)) {
          nomeZ<- object$nameUV[[3]]
          if(length(rev.sgn)==1) rev.sgn<-rep(rev.sgn,length(nomeZ))
          }
             else {
                if(! parm %in% object$nameUV[[3]]) {stop("invalid parm")}
                  else {nomeZ<-parm}
                  }
        if(length(rev.sgn)!=length(nomeZ)) rev.sgn<-rep(rev.sgn, length.out=length(nomeZ))
        rr<-list()
        z<-if("lm"%in%class(object)) abs(qt((1-level)/2,df=object$df.residual)) else abs(qnorm((1-level)/2))
        for(i in 1:length(nomeZ)){ #per ogni variabile segmented `parm' (tutte o selezionata)..
            #nomi.U<-grep(paste("\\.",nomeZ[i],"$",sep=""),object$nameUV$U,value=TRUE)
            #nomi.V<-grep(paste("\\.",nomeZ[i],"$",sep=""),object$nameUV$V,value=TRUE)
            nomi.U<- object$nameUV$U[f.U(object$nameUV$U, nomeZ[i])]
            nomi.V<- object$nameUV$V[f.U(object$nameUV$V, nomeZ[i])]
            m<-matrix(,length(nomi.U),3)
            rownames(m)<-nomi.V
            colnames(m)<-c("Est.",paste("CI","(",level*100,"%",")",c(".l",".u"),sep=""))
            for(j in 1:length(nomi.U)){ #per ogni psi della stessa variabile segmented..
                    sel<-c(nomi.V[j],nomi.U[j])
                    V<-vcov(object,var.diff=var.diff)[sel,sel] #questa e' vcov di (psi,U)
                    b<-coef(object)[sel[2]] #diff-Slope
                    th<-c(b,1)
                    orig.coef<-drop(diag(th)%*%coef(object)[sel]) #sono i (gamma,beta) th*coef(ogg)[sel]
                    gammma<-orig.coef[1]
                    est.psi<-object$psi[sel[1],2]
                    V<-diag(th)%*%V%*%diag(th) #2x2 vcov() di gamma e beta
                    se.psi<-sqrt((V[1,1]+V[2,2]*(gammma/b)^2-2*V[1,2]*(gammma/b))/b^2)
                    r<-c(est.psi, est.psi-z*se.psi, est.psi+z*se.psi)
                    if(rev.sgn[i]) r<-c(-r[1],rev(-r[2:3]))
                    m[j,]<-r
                    } #end loop j (ogni psi della stessa variabile segmented)
            #CONTROLLA QUESTO:..sarebbe piu' bello
            if(nrow(m)==1) rownames(m)<-"" else m<-m[order(m[,1]),]
            if(rev.sgn[i]) {
                #m<-m[nrow(m):1,]
                rownames(m)<-rev(rownames(m))
                }
            rr[[length(rr)+1]]<-signif(m,digits)
            } #end loop i (ogni variabile segmented)
        names(rr)<-nomeZ
        return(rr)
          } #end_function


