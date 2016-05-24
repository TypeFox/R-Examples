bigRR.formula <-
function(formula=NULL, y=NULL, X=NULL, Z=NULL, data=NULL,shrink=NULL,weight = NULL,
		family = gaussian(link = identity), lambda = NULL,impute = FALSE, 
		tol.err = 1e-6, tol.conv = 1e-8, only.estimates = FALSE, GPU = FALSE,...) {
Call<-match.call()
la<-lambda
o.er<-tol.err
o.es<-only.estimates
trms<-all.vars(formula)
if(!is.null(shrink)){
if(!is.vector(shrink)) stop("shrink must be numeric or a character vector")
  }
if(is.numeric(shrink)) {
  if(length(shrink)>=(length(trms))) {
      if(!("."%in%trms)) stop("Number of shrinkage parameter exceeds total number of parameters")
      }
  if("."%in%trms){
      if(is.null(data)) stop("A formula containing a dot (.) must accopmpany with a data frame")
      if(ncol(data)<=length(trms)) stop("Number of shrinkage parameter exceeds total number of parameters")
      vnames<-names(data)
      skvar<-vnames[shrink]
      } else{
       skvar<-trms[(shrink+1)]
      }

  } else{
  skvar<-shrink
  if(!(all(skvar)%in%trms)) stop("Parameter to shrink is not in the model formula")
  }
  fm1<-paste(c(deparse(formula),paste(skvar,collapse="-")),collapse="-")
  fmf<-as.formula(fm1)
  fmr<-as.formula(paste(c("~",paste(skvar,collapse="+"),"-1"),collapse=""))
  fm<-model.frame(fmf,data=data)
  X<-model.matrix(attr(fm,"terms"),data=fmf)
  row.names(X)<-NULL
  y<-model.response(fm)
  fz<-model.frame(fmr,data)
  Z<-model.matrix(attr(fz,"terms"),data=fmr)
  row.names(Z)<-NULL
res<-bigRR(y=y,X=X,Z=Z,lambda = la,tol.err = o.er, only.estimates=o.es)
res$Call<-Call
return(res)
}

