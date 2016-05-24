CLIP.profile <-
function(obj=NULL, variable, data, which,  firth=TRUE, weightvar, control=logistf.control(), 
  offset=NULL, from=NULL, to=NULL, steps=101, legacy=FALSE, keep=FALSE){
    require(mice)
#    require(session)      ###?? required?
    old<-legacy
    if (!is.null(obj)) {
       if (is.mira(obj)) {
          # assuming a mira object with fits on each imputed data set
#          data.imp<-texteval(paste("",obj$call$data,"",sep=""))
          data.imp<-eval(obj$call$data)    #############GEHT NOCH NICHT!!
          nimp<-data.imp$m          
          data<-lapply(1:nimp, function(x) complete(data.imp, action=x))
          fits<-obj$analyses
          formula<-as.formula(fits[[1]]$call$formula)
          }
       else {
          # assuming as input a list of data sets (data) and a list of fits (obj)
          fits<-obj
          if(missing(data)) if(is.null(fits[[1]]$data)) stop("Please provide data either as list of imputed data sets or by calling logistf on the imputed data sets with dataout=TRUE.\n")
                            else data<-lapply(1:length(fits), function(X) fits[[X]]$data)

          formula<-as.formula(fits[[1]]$call$formula)
          nimp<-length(data)
          } 
       }
     else  {
            nimp<-length(data)
            stop(paste("Please provide a list of fit objects or a mira object with fits on the ",nimp," imputed data sets.\n"))
            }     
                 
   
     
 
     imputations<-nimp
     if (missing(which) & missing(variable))
        stop("You must specify which (a one-sided formula) or variable (the name of the variable).")
    if (missing(control))
        control <- logistf.control()
    res<-matrix(0,length(grid),3)
    
    nperimp<-unlist(lapply(1:imputations,function(x) nrow(data[[x]])))
    
    variable.names<-colnames(model.matrix(formula, data = data[[1]]))
#    mat.data<-matrix(lapply(1:imputations,function(x) matrix(unlist(data[[x]]),nperimp[[x]],ncol(data[[x]]))),sum(nperimp),ncol(data[[1]]))
    

    imputation.indicator<-unlist(sapply(1:imputations, function(x) rep(x,nperimp[x]))[TRUE])       #[TRUE]makes vector out of matrix
    
     mat.data<-matrix(0,sum(nperimp),ncol(data[[1]]))   ### copy list of data set into a matrix
     for(i in 1:imputations) mat.data[imputation.indicator==i,1:ncol(data[[i]])]<-as.matrix(data[[i]])

    
  #  if(missing(weightvar)) {
  #     weightvar<-"weights"
  #     mat.data[,ncol(mat.data)]<-rep(1,nrow(mat.data))
  #   }
    big.data<-data.frame(mat.data)
    colnames(big.data)<-colnames(data[[1]])

     
    k<-length(variable.names)  
    xyw<-matrix(0,sum(nperimp),k+2)
    xyw[,1:k]<-  model.matrix(formula, data = big.data)
    xyw[,k+1]<- model.response(model.frame(as.formula(formula),data=big.data))
    if (!missing(which)) {
     cov.name <- variable.names
     cov.name2 <- colnames(model.matrix(which, data = big.data))
     pos <- match(cov.name2, cov.name)
     variable <- cov.name2
    }
    else {
     cov.name <- variable.names 
     cov.name2 <- variable
     pos <- match(cov.name2, cov.name)
     }
#    fit <- logistf.fit(x, y, weight = weight, offset = offset,
#        firth = firth, control = control)
#    std.pos <- diag(fit$var)[pos]^0.5
#    coefs <- fit$beta
#    covs <- fit$var
#    n <- nrow(x)
#    cov.name <- labels(x)[[2]]

if (missing(weightvar)) {
# data<-lapply(1:imputations, function(z) {
#  data[[z]]$weightvar<-rep(1,nrow(data[[z]]))
#  data[[z]]
#  })
#  weightvar<-"weightvar"
  xyw[,k+2]<-rep(1,nrow(xyw))
  }
else xyw[,k+2]<-big.data[,weightvar]
posweight<-k+2
  
 
lf<-function(index) logistf.fit(y=xyw[index,k+1], x=xyw[index,1:k], weight=xyw[index,k+2])
                       
if (is.null(fits))   fits<-lapply(1:imputations, function(z) lf(imputation.indicator==z))

if(is.null(from)){
 lower.collect<-        (unlist(lapply(1:imputations,function(x) fits[[x]]$ci.lower[pos])))
 from<-min(lower.collect) ###von dem den index nehmen und davon das PL CI ausrechnen
} 

if(is.null(to)){
 upper.collect<-        (unlist(lapply(1:imputations,function(x) fits[[x]]$ci.upper[pos])))
 to<-max(upper.collect)
} 



estimate<-mean(unlist(lapply(1:imputations,function(x) fits[[x]]$coefficients[pos])))

iter<-numeric(0)

loglik<-unlist(lapply(1:imputations, function(x) fits[[x]]$loglik[2]))
beta<-t(matrix(unlist(lapply(1:imputations,function(x) fits[[x]]$coefficients)),k,imputations))


lpdf<-function(zz,z) logistf.pdf(x=xyw[imputation.indicator==zz,1:k], y=xyw[imputation.indicator==zz,k+1], 
   weight=xyw[imputation.indicator==zz,k+2], beta=beta[zz,],loglik=loglik[zz],
   pos=pos, firth=firth, offset=offset, control=control, b=z, old=old)$pdf

 z_seq<-seq(from, to, (to-from)/steps)

 if(keep==FALSE){
   f=function(z)  mean(unlist(lapply(1:imputations, function(zz) lpdf(zz,z))))
   ### lasse z laufen von from nach to
   ### evaluiere f an allen z's und errechne daraus profile
   pdf_mat<-NULL
   profile.mat<-NULL
   pdf_seq<-unlist(lapply(z_seq, function(Z) f(Z)))
   } else {
   f=function(z)  unlist(lapply(1:imputations, function(zz) lpdf(zz,z)))
#   pdf_mat<-matrix(0,steps+1,imputations)
   pdf_mat<-matrix(unlist(lapply(z_seq, function(Z) f(Z))), imputations, steps+1)   
   profile.mat<- -qnorm(pdf_mat)**2
   pdf_seq<-apply(pdf_mat,2,mean)
   }
   

 chisq_seq<- -qnorm(pdf_seq)**2
 res<-list(beta=z_seq, cdf=pdf_seq, profile=chisq_seq, cdf.matrix=pdf_mat, profile.matrix=profile.mat, call=match.call())
 attr(res,"class")<-c("logistf.CLIP.profile","logistf.profile")
 res
}

