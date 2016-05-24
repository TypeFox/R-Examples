predict.segmented<-function(object, newdata, ...){
#rev: 30/10/2013: it seems to work correctly, even with the minus variable (null right slope..)
#rev: 14/4/2014 now it works like predict.lm/glm
#BUT problems if type="terms" (in realta' funziona, il problema e' che
#     restituisce una colonna per "x", "U.x", "psi.x".. (Eventualmente si dovrebbero sommare..)
#if(!is.null(object$orig.call$offset)) stop("predict.segmented can not handle argument 'offset'. Include it in formula!")

  dummy.matrix<-function(x.values, x.name, obj.seg, psi.est=TRUE){
    #given the segmented fit 'obj.seg' and a segmented variable x.name with corresponding values x.values,
    #this function simply returns a matrix with columns (x, (x-psi)_+, -b*I(x>psi))
    #or  ((x-psi)_+, -b*I(x>psi)) if obj.seg does not include the coef for the linear "x"
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
        n<-length(x.values)
        #le seguenti righe selezionavano (ERRONEAMENTE) sia "U1.x" sia "U1.neg.x" (se "x" e "neg.x" erano segmented covariates)
        #nameU<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$U, value = TRUE)
        #nameV<- grep(paste("\\.",x.name,"$", sep=""), obj.seg$nameUV$V, value = TRUE)
        nameU<-obj.seg$nameUV$U[f.U(obj.seg$nameUV$U,x.name)]
        nameV<-obj.seg$nameUV$V[f.U(obj.seg$nameUV$V,x.name)]

        diffSlope<-coef(obj.seg)[nameU]
        est.psi<-obj.seg$psi[nameV,2]

        k<-length(est.psi)

        PSI <- matrix(rep(est.psi, rep(n, k)), ncol = k)
        newZ<-matrix(x.values, nrow=n,ncol=k, byrow = FALSE)

        dummy1<-pmax(newZ-PSI,0)
        
        if(psi.est){
          V<-ifelse(newZ>PSI,-1,0)
          dummy2<- if(k==1) V*diffSlope  else V%*%diag(diffSlope) #t(diffSlope*t(-I(newZ>PSI)))
          newd<-cbind(x.values,dummy1,dummy2)
          colnames(newd)<-c(x.name,nameU, nameV)
          } else {
          newd<-cbind(x.values,dummy1)
          colnames(newd)<-c(x.name,nameU)
          }
        if(!x.name%in%names(coef(obj.seg))) newd<-newd[,-1,drop=FALSE]
        return(newd)
    }
#--------------------------------------------------------------
  if(missing(newdata)){
    newd.ok<-model.frame(object)
  } else {
  #devi trasformare la variabili segmented attraverso dummy.matrix()
  nameU<-object$nameUV$U
  nameV<-object$nameUV$V
  nameZ<-object$nameUV$Z
  n<-nrow(newdata)
  r<-NULL
  for(i in 1:length(nameZ)){
      x.values<-newdata[[nameZ[i]]]
      DM<-dummy.matrix(x.values, nameZ[i], object)
      r[[i]]<-DM
      }  
  newd.ok<-data.frame(matrix(unlist(r), nrow=n, byrow = FALSE))
  names(newd.ok)<- unlist(sapply(r, colnames))
  idZ<-match(nameZ, names( newdata))
  newdata<-cbind(newdata[,-idZ, drop=FALSE], newd.ok) #  newdata<-subset(newdata, select=-idZ)
  newdata<-cbind(newdata, newd.ok)
  }
  class(object)<-class(object)[-1]
  f<-predict(object, newdata=newdata, ...)
  #f<-if(inherits(object, what = "glm", which = FALSE)) predict.glm(object, newdata=newd.ok, ...) else predict.lm(object, newdata=newd.ok, ...)
  return(f)
#sommare se "terms"?
  }  
  
  
  
  
