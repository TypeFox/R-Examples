
W.inv<- function(W, symmetric=TRUE,inverse=TRUE){
   eW <- eigen(W, symmetric=symmetric)
   d <- eW$values
   if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
       stop("'W' is not positive definite")
   else d[d<=0]<- ifelse(inverse, Inf, 0)
   A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
   A # t(A)%*%A = W^{-1}
}

# adapted from lm.gls in MASS
lmGls<- function (formula, data, A, ...) {
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$A <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    yy <- model.response(m)
       y<- A%*%yy
    xx <- model.matrix(Terms, m, contrasts)
       x<- A%*%xx; colnames(x)[1]<- "(Intercept)"
    dtf<- data.frame(y=y,x)
    fit<- lm(y~.-1, data=dtf, ...)

    fit
}

# generalized least squares test
scanOne.0 <-
   function(y,
            x,
            prdat,
            cov,
            intcovar = NULL,
            test = c("None","F","Chisq"))
{
# prdat$pr: n by ? by ? matrix, allele probabilities
# vc: object from estVC or aicVC
# test: "Chisq", "F" or "Cp"
   diag.cov<- diag(cov)
   if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
      if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
         weights<- NULL
      }else weights<- 1/diag.cov
   }else weights<- NA
   gcv<- W.inv(cov)
   test<- match.arg(test)

   nsnp<- dim(prdat$pr)[3]
   if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
   model.par<- vector("list",nsnp)
      names(model.par)<- prdat$snp
   P<- rep(Inf,nsnp)
      names(P)<- prdat$snp
   V<- P
   if(is.null(intcovar)){
      if(!missing(x)){
         oTmp<- data.frame(y=y,x)
      }else{
         oTmp<- data.frame(y=y)
      }
      if( !is.null(weights[1]) && is.na(weights[1]) ){
         g0<- lmGls(y~.,data=oTmp,A=gcv)
      }else{
         g0<- lm(y~.,data=oTmp,weights=weights)
      }
      if(test=="None"){
         P0<- logLik(g0)
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,prdat$pr[,-1,k])
            }else{
               oTmp<- data.frame(y=y,a=prdat$pr[,-1,k])
            }

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(y~.,data=oTmp,A=gcv)
            }else{
               g<- lm(y~.,data=oTmp,weights=weights)
            }
            model.par[[k]]<- g$coef
            P[k]<- logLik(g)
            V[k]<- sum(g$res^2)
         }
         P<- 2*(P-P0)
      }else{
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,prdat$pr[,-1,k])
            }else{
               oTmp<- data.frame(y=y,prdat$pr[,-1,k])
            }

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(y~.,data=oTmp,A=gcv)
            }else{
               g<- lm(y~.,data=oTmp,weights=weights)
            }
            model.par[[k]]<- g$coef
            P[k]<- anova(g0,g,test=test)$P[2]
            V[k]<- sum(g$res^2)
         }
      }
      V<- sum(g0$res^2) - V
         V<- V/sum(anova(g0)[,"Sum Sq"])
   }else{
      if(!missing(x)){
         oTmp<- data.frame(y=y,x,intcovar)
      }else{
         oTmp<- data.frame(y=y,intcovar)
      }
      if( !is.null(weights[1]) && is.na(weights[1]) ){
         g0<- lmGls(y~.,data=oTmp,A=gcv)
      }else{
         g0<- lm(y~.,data=oTmp,weights=weights)
      }
      if(test=="None"){
         P0<- logLik(g0)
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,prdat$pr[,-1,k])
            }else{
               oTmp<- data.frame(y=y,intcovar,prdat$pr[,-1,k])
            }
            nc<- ncol(oTmp); nq<- ncol(prdat$pr[,-1,k])-1
            str<- paste(paste("(",paste(colnames(oTmp)[nc-nq-(nint:1)],collapse="+"),")",sep=""),
                        paste("(",paste(colnames(oTmp)[(nc-nq):nc],collapse="+"),")",sep=""),
                        sep=":")
            str<- paste("y~.+",str,sep="")

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(formula(str),data=oTmp,A=gcv)
            }else{
               g<- lm(formula(str),data=oTmp,weights=weights)
            }
            model.par[[k]]<- g$coef
            P[k]<- logLik(g)
            V[k]<- sum(g$res^2)
         }
         P<- 2*(P-P0)
      }else{
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,prdat$pr[,-1,k])
            }else{
               oTmp<- data.frame(y=y,intcovar,prdat$pr[,-1,k])
            }
            nc<- ncol(oTmp); nq<- ncol(prdat$pr[,-1,k])-1
            str<- paste(paste("(",paste(colnames(oTmp)[nc-nq-(nint:1)],collapse="+"),")",sep=""),
                        paste("(",paste(colnames(oTmp)[(nc-nq):nc],collapse="+"),")",sep=""),
                        sep=":")
            str<- paste("y~.+",str,sep="")

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(formula(str),data=oTmp,A=gcv)
            }else{
               g<- lm(formula(str),data=oTmp,weights=weights)
            }
            model.par[[k]]<- g$coef
            P[k]<- anova(g0,g,test=test)$P[2]
            V[k]<- sum(g$res^2)
         }
      }
      V<- sum(g0$res^2) - V
         V<- V/sum(anova(g0)[,"Sum Sq"])
   }

   list(snp=prdat$snp,
        chr=prdat$chr,
        dist=prdat$dist,
        p=P,
        v=V*100,
        parameters=model.par)
}

scanOne.1 <-
   function(y,
            x,
            prdat,
            cov,
            intcovar = NULL,
            test = c("None","F","Chisq"))
{
# prdat$pr: n by 3 by ? matrix, conditional probabilities
# vc: object from estVC or aicVC
# test: "Chisq", "F" or "Cp"
   diag.cov<- diag(cov)
   if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
      if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
         weights<- NULL
      }else weights<- 1/diag.cov
   }else weights<- NA
   gcv<- W.inv(cov)
   test<- match.arg(test)

   nsnp<- dim(prdat$pr)[3]
   if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
   model.par<- vector("list",nsnp)
      names(model.par)<- prdat$snp
   P<- rep(Inf,nsnp)
      names(P)<- prdat$snp
   V<- P
   if(is.null(intcovar)){
      if(!missing(x)){
         oTmp<- data.frame(y=y,x)
      }else{
         oTmp<- data.frame(y=y)
      }
      if( !is.null(weights[1]) && is.na(weights[1]) ){
         g0<- lmGls(y~.,data=oTmp,A=gcv)
      }else{
         g0<- lm(y~.,data=oTmp,weights=weights)
      }
      if(test=="None"){
         P0<- logLik(g0)
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }else{
               oTmp<- data.frame(y=y,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(y~.,data=oTmp,A=gcv)
            }else{
               g<- lm(y~.,data=oTmp,weights=weights)
            }
            model.par[[k]]<- g$coef
            P[k]<- logLik(g)
            V[k]<- sum(g$res^2)
         }
         P<- 2*(P-P0)
      }else{
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }else{
               oTmp<- data.frame(y=y,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(y~.,data=oTmp,A=gcv)
            }else{
               g<- lm(y~.,data=oTmp,weights=weights)
            }
            model.par[[k]]<- g$coef
            P[k]<- anova(g0,g,test=test)$P[2]
            V[k]<- sum(g$res^2)
         }
      }
      V<- sum(g0$res^2) - V
         V<- V/sum(anova(g0)[,"Sum Sq"])
   }else{
      if(!missing(x)){
         oTmp<- data.frame(y=y,x,intcovar)
      }else{
         oTmp<- data.frame(y=y,intcovar)
      }
      if( !is.null(weights[1]) && is.na(weights[1]) ){
         g0<- lmGls(y~.,data=oTmp,A=gcv)
      }else{
         g0<- lm(y~.,data=oTmp,weights=weights)
      }
      if(test=="None"){
         P0<- logLik(g0)
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }else{
               oTmp<- data.frame(y=y,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }
            nc<- ncol(oTmp); nq<- ncol(prdat$pr[,-1,k])-1
            str<- paste(paste("(",paste(colnames(oTmp)[nc-nq-(nint:1)],collapse="+"),")",sep=""),
                        paste("(",paste(colnames(oTmp)[(nc-nq):nc],collapse="+"),")",sep=""),
                        sep=":")
            str<- paste("y~.+",str,sep="")

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(formula(str),data=oTmp,A=gcv)
            }else{
               g<- lm(formula(str),data=oTmp,weights=weights)
            }
            model.par[[k]]<- g$coef
            P[k]<- logLik(g)
            V[k]<- sum(g$res^2)
         }
         P<- 2*(P-P0)
      }else{
         for(k in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }else{
               oTmp<- data.frame(y=y,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
            }
            nc<- ncol(oTmp); nq<- ncol(prdat$pr[,-1,k])-1
            str<- paste(paste("(",paste(colnames(oTmp)[nc-nq-(nint:1)],collapse="+"),")",sep=""),
                        paste("(",paste(colnames(oTmp)[(nc-nq):nc],collapse="+"),")",sep=""),
                        sep=":")
            str<- paste("y~.+",str,sep="")

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(formula(str),data=oTmp,A=gcv)
            }else{
               g<- lm(formula(str),data=oTmp,weights=weights)
            }
            model.par[[k]]<- g$coef
            P[k]<- anova(g0,g,test=test)$P[2]
            V[k]<- sum(g$res^2)
         }
      }
      V<- sum(g0$res^2) - V
         V<- V/sum(anova(g0)[,"Sum Sq"])
   }

   list(snp=prdat$snp,
        chr=prdat$chr,
        dist=prdat$dist,
        p=P,
        v=V*100,
        parameters=model.par)
}

scanOne.2 <-
   function(y,
            x,
            gdat,
            cov,
            intcovar = NULL,
            numGeno = FALSE,
            test = c("None","F","Chisq"))
{
# gdat: n by ? matrix, marker data. Markers in columes!!!
# vc: object from estVC or aicVC
# intcover: covariates that interact with QTL
# test: "Chisq", "F" or "Cp"
   diag.cov<- diag(cov)
   if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
      if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
         weights<- NULL
      }else weights<- 1/diag.cov
   }else weights<- NA
   gcv<- W.inv(cov)
   test<- match.arg(test)
   if(numGeno){
      num.geno<- I
   }else num.geno<- as.factor

   nsnp<- dim(gdat)[2]
   if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
   model.par<- vector("list",nsnp)
      names(model.par)<- colnames(gdat)
   P<- rep(Inf,nsnp)
      names(P)<- colnames(gdat)
   V<- P
   if(is.null(intcovar)){
      if(!missing(x)){
         oTmp<- data.frame(y=y,x)
      }else{
         oTmp<- data.frame(y=y)
      }
      if( !is.null(weights[1]) && is.na(weights[1]) ){
         g0<- lmGls(y~.,data=oTmp,A=gcv)
      }else{
         g0<- lm(y~.,data=oTmp,weights=weights)
      }
      if(test=="None"){
         P0<- logLik(g0)
         for(j in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,snp=num.geno(gdat[,j]))
            }else{
               oTmp<- data.frame(y=y,snp=num.geno(gdat[,j]))
            }

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(y~.,data=oTmp,A=gcv)
            }else{
               g<- lm(y~.,data=oTmp,weights=weights)
            }
            model.par[[j]]<- g$coef
            P[j]<- logLik(g)
            V[j]<- sum(g$res^2)
         }
         P<- 2*(P - P0)
      }else{
         for(j in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,snp=num.geno(gdat[,j]))
            }else{
               oTmp<- data.frame(y=y,snp=num.geno(gdat[,j]))
            }

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(y~.,data=oTmp,A=gcv)
            }else{
               g<- lm(y~.,data=oTmp,weights=weights)
            }
            model.par[[j]]<- g$coef
            P[j]<- anova(g0,g,test=test)$P[2]
            V[j]<- sum(g$res^2)
         }
      }
      V<- sum(g0$res^2) - V
         V<- V/sum(anova(g0)[,"Sum Sq"])
   }else{
      if(!missing(x)){
         oTmp<- data.frame(y=y,x,intcovar)
      }else{
         oTmp<- data.frame(y=y,intcovar)
      }
      if( !is.null(weights[1]) && is.na(weights[1]) ){
         g0<- lmGls(y~.,data=oTmp,A=gcv)
      }else{
         g0<- lm(y~.,data=oTmp,weights=weights)
      }
      if(test=="None"){
         P0<- logLik(g0)
         for(j in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,snp=num.geno(gdat[,j]))
            }else{
               oTmp<- data.frame(y=y,intcovar,snp=num.geno(gdat[,j]))
            }

            nc<- ncol(oTmp)
            str<- paste(colnames(oTmp)[nc-(nint:1)],colnames(oTmp)[nc],collapse="+",sep=":")
            str<- paste("y~.+",str,sep="")

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(formula(str),data=oTmp,A=gcv)
            }else{
               g<- lm(formula(str),data=oTmp,weights=weights)
            }
            model.par[[j]]<- g$coef
            P[j]<- logLik(g)
            V[j]<- sum(g$res^2)
         }
         P<- 2*(P-P0)
      }else{
         for(j in 1:nsnp){
            if(!missing(x)){
               oTmp<- data.frame(y=y,x,intcovar,snp=num.geno(gdat[,j]))
            }else{
               oTmp<- data.frame(y=y,intcovar,snp=num.geno(gdat[,j]))
            }

            nc<- ncol(oTmp)
            str<- paste(colnames(oTmp)[nc-(nint:1)],colnames(oTmp)[nc],collapse="+",sep=":")
            str<- paste("y~.+",str,sep="")

            if( !is.null(weights[1]) && is.na(weights[1]) ){
               g<- lmGls(formula(str),data=oTmp,A=gcv)
            }else{
               g<- lm(formula(str),data=oTmp,weights=weights)
            }
            model.par[[j]]<- g$coef
            P[j]<- anova(g0,g,test=test)$P[2]
            V[j]<- sum(g$res^2)
         }
      }
      V<- sum(g0$res^2) - V
         V<- V/sum(anova(g0)[,"Sum Sq"])
   }

   list(p=P,
        v=V*100,
        parameters=model.par)
}

scanOne<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intcovar = NULL,
            numGeno = FALSE,
            test = c("None","F","Chisq"),
            minorGenoFreq = 0,
            rmv = TRUE)

{
   if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.")
   if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
         stop("x: missing or infinite data points not allowed.")
   UseMethod("scanOne")
}

scanOne.default<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            intcovar = NULL,
            numGeno = FALSE,
            test = c("None","F","Chisq"),
            minorGenoFreq = 0,
            rmv = TRUE)
{
   if(!is.null(vc)){
      if(is.element("bgv",attr(vc,"class"))){
         nb<- length(vc$par) - sum(vc$nnl)
         nr<- nrow(vc$y)
         cov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:vc$nv)
            if(vc$nnl[i]) cov<- cov + vc$v[[i]]*vc$par[nb+vc$nn[i]]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("vc should be a matrix.")
         if(!is.numeric(vc)) stop("vc should be a numeric matrix.")
         cov<- vc
      }
   }else cov<- diag(nrow(as.matrix(y)))
   if(!is.null(prdat)){
      if(is.element("addEff",class(prdat))){
         pv<- scanOne.0(y=y,x=x,prdat=prdat,cov=cov,intcovar=intcovar,test=test)
      }else{
         pv<- scanOne.1(y=y,x=x,prdat=prdat,cov=cov,intcovar=intcovar,test=test)
      }
   }else{
      if(any(is.na(gdat)))
         stop("There are missing genotypes...")
      tb<- sort(union(as.matrix(gdat),NULL))
      tbf<- NULL
      for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
         if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
      tbf<- apply(tbf,2,min)
      idx<- (tbf < nrow(gdat)*minorGenoFreq)
      if(sum(idx)>0){
         if(rmv){
            gdat<- as.matrix(gdat)
            tb<- sort(union(as.matrix(gdat),NULL))
            tbf<- NULL
            for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
               if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
            tbf<- apply(tbf,2,min)
            idx<- (tbf < nrow(gdat)*minorGenoFreq)
            gdat<- gdat[,!idx]
            rm(tb,tbf,ii,idx)
         }else{
            cat("minor genotype frequency is too small at one or more SNPs.\n")
            return(NULL)
        }
      }

      gdat<- as.data.frame(gdat)
      pv<- scanOne.2(y=y,x=x,gdat=gdat,cov=cov,intcovar=intcovar,numGeno=numGeno, test=test)
   }

   class(pv)<- c("scanOne",test)
   pv
}

print.scanOne <-
   function(x,...)
{
   tt<- x; class(tt)<- NULL
      tt$parameters<- NULL
      tt<- as.data.frame(tt)
   if(length(tt$p)>5){
      cat("Test statistic:\n")
      print(tt[1:5,])
      cat("... ...\n\n")

      cat("Coefficients:\n")
      print(x$par[1:5])
      cat("... ...\n\n")
   }else{
      cat("Test statistic:\n")
      print(tt)
      cat("\n")

      cat("Coefficients:\n")
      print(x$par)
   }
}

scanTwo.1 <-
   function(y,
            x,
            prdat,
            cov)
{
   diag.cov<- diag(cov)
   if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
      if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
         weights<- NULL
      }else weights<- 1/diag.cov
   }else weights<- NA
   gcv<- W.inv(cov)
   nsnp<- dim(prdat$pr)[3]
   P<- matrix(NA,nrow=nsnp,ncol=nsnp)
      rownames(P)<- colnames(P)<- prdat$snp
   if(nsnp<1) return(NULL)

   if(!missing(x)){
      oTmp.xy<- data.frame(y=y,x)
   }else{
      oTmp.xy<- data.frame(y=y)
   }
   for(i in 1:(nsnp-1)){
      for(k in (i+1):nsnp){
         xTmp<- data.frame(a1=prdat$pr[,1,i] - prdat$pr[,3,i],
                           d1=prdat$pr[,2,i],
                           a2=prdat$pr[,1,k] - prdat$pr[,3,k],
                           d2=prdat$pr[,2,k])
         oTmp<- cbind(oTmp.xy, oTmp)

         if( !is.null(weights[1]) && is.na(weights[1]) ){
            g0<- lmGls(y~.,data=oTmp,A=gcv)
         }else{
            g0<- lm(y~.,data=oTmp,weights=weights)
         }
         if( !is.null(weights[1]) && is.na(weights[1]) ){
            g<- lmGls(y~(a1+d1)*(a2+d2) + .,data=oTmp,A=gcv)
         }else{
            g<- lm(y~(a1+d1)*(a2+d2) + .,data=oTmp,weights=weights)
         }
         P[i,k]<- 2*(logLik(g)-logLik(g0))
      }
   }

   P
}

scanTwo.2 <-
   function(y,
            x,
            gdat,
            cov,
            numGeno)
{
   if(numGeno){
      num.geno<- I
   }else num.geno<- as.factor
   diag.cov<- diag(cov)
   if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
      if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
         weights<- NULL
      }else weights<- 1/diag.cov
   }else weights<- NA
   gcv<- W.inv(cov)

   nsnp<- dim(gdat)[2]
   P<- matrix(NA,nrow=nsnp,ncol=nsnp)
      rownames(P)<- colnames(P)<- colnames(gdat)
   if(nsnp<1) return(NULL)

   if(!missing(x)){
      oTmp.xy<- data.frame(y=y,x)
   }else{
      oTmp.xy<- data.frame(y=y)
   }
   for(i in 1:(nsnp-1)){
      for(j in (i+1):nsnp){
         oTmp<- data.frame(snp1=num.geno(gdat[,i]),
                           snp2=num.geno(gdat[,j]))
         oTmp<- cbind(oTmp.xy, oTmp)

         if( !is.null(weights[1]) && is.na(weights[1]) ){
            g0<- lmGls(y~.,data=oTmp,A=gcv)
         }else{
            g0<- lm(y~.,data=oTmp,weights=weights)
         }
         if( !is.null(weights[1]) && is.na(weights[1]) ){
            g<- lmGls(y~snp1*snp2 + .,data=oTmp,A=gcv)
         }else{
            g<- lm(y~snp1*snp2 + .,data=oTmp,weights=weights)
         }
         P[i,j]<- 2*(logLik(g)-logLik(g0))
      }
   }

   P
}

scanTwo<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            numGeno = FALSE,
            minorGenoFreq = 0,
            rmv = TRUE)

{
   if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.")
   if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
         stop("x: missing or infinite data points not allowed.")
   UseMethod("scanTwo")
}

scanTwo.default<- 
   function(y,
            x,
            gdat,
            prdat = NULL,
            vc = NULL,
            numGeno = FALSE,
            minorGenoFreq = 0,
            rmv = TRUE)
{
   if(!is.null(vc)){
      if(is.element("bgv",attr(vc,"class"))){
         nb<- length(vc$par) - sum(vc$nnl)
         nr<- nrow(vc$y)
         cov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:vc$nv)
            if(vc$nnl[i]) cov<- cov + vc$v[[i]]*vc$par[nb+vc$nn[i]]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("vc should be a matrix.")
         if(!is.numeric(vc)) stop("vc should be a numeric matrix.")
         cov<- vc
      }
   }else cov<- diag(nrow(as.matrix(y)))
   if(!is.null(prdat)){
      pv<- scanTwo.1(y=y,x=x,prdat=prdat,cov=cov)
   }else{
      if(any(is.na(gdat)))
         stop("There are missing genotypes...")
      tb<- sort(union(as.matrix(gdat),NULL))
      tbf<- NULL
      for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
         if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
      tbf<- apply(tbf,2,min)
      idx<- (tbf < nrow(gdat)*minorGenoFreq)
      if(sum(idx)>0){
         if(rmv){
            gdat<- as.matrix(gdat)
            tb<- sort(union(as.matrix(gdat),NULL))
            tbf<- NULL
            for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
               if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
            tbf<- apply(tbf,2,min)
            idx<- (tbf < nrow(gdat)*minorGenoFreq)
            gdat<- gdat[,!idx]
            rm(tb,tbf,ii,idx)
         }else{
            cat("minor genotype frequency is too small at one or more SNPs.\n")
            return(NULL)
        }
      }

      gdat<- as.data.frame(gdat)
      pv<- scanTwo.2(y=y,x=x,gdat=gdat,cov=cov,numGeno=numGeno)
   }

   class(pv)<- "scanTwo"
   pv
}

# generalized least squares estimates
gls<- function(formula,data=NULL,vc=NULL){
   if(is.null(data)){
      xx<- model.matrix(formula)
   }else xx<- model.matrix(formula,data)
   yy<- model.response(model.frame(formula,data))

   nr<- nrow(xx)
   if(!is.null(vc)){
      if(is.element("bgv",attr(vc,"class"))){
         nb<- length(vc$par) - sum(vc$nnl)
         nr<- nrow(vc$y)
         cov<- matrix(0,nrow=nr,ncol=nr)
         for(i in 1:vc$nv)
            if(vc$nnl[i]) cov<- cov + vc$v[[i]]*vc$par[nb+vc$nn[i]]
      }else{
         if(is.data.frame(vc)) vc<- as.matrix(vc)
         if(!is.matrix(vc)) stop("vc should be a matrix.")
         if(!is.numeric(vc)) stop("vc should be a numeric matrix.")
         cov<- vc
      }
   }else cov<- diag(nrow(as.matrix(yy)))
   A<- W.inv(cov)

   x<- A%*%xx; colnames(x)[1]<- "(Intercept)"
   y<- A%*%yy
   dtf<- data.frame(y=y,x)
   mdl<- lm(y~.-1, data=dtf)
   mdl$data<- dtf

#   print(logLik(mdl))
   summary(mdl)$coeff
}

