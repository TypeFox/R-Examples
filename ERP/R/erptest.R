erptest <-
function(dta,design,design0=NULL,method="BH",alpha=0.05,pi0=1) {
   
   method = match.arg(method,choices=c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"))
   if (is.null(design0)) design0 = matrix(1,nrow=nrow(dta),ncol=1)
   erpdta = as.matrix(dta)
   design = as.matrix(design)
   design0 = as.matrix(design0)
   if (typeof(erpdta)!="double") stop("ERPs should be of type double")
   if (nrow(erpdta)!=nrow(design)) stop("dta and design should have the same number of rows")
   if (nrow(erpdta)!=nrow(design0)) stop("dta and design0 should have the same number of rows")
   if (ncol(design)<=ncol(design0)) stop("design0 should have fewer columns than design")
   idsignal = NULL
   for (j in 1:ncol(design)) {
      cj = apply(design0,2,function(x,y) all(x==y),y=design[,j])
      if (all(!cj)) idsignal = c(idsignal,j)
   }
   if (length(idsignal)!=(ncol(design)-ncol(design0))) stop("the null model design0 should be nested into the non-null model design")
   if (typeof(alpha)!="double") stop("alpha should be of type double")
   if ((alpha<=0)|(alpha>=1)) stop("alpha should be in ]0,1[, typically 0.05")
   if (typeof(pi0)!="double") stop("pi0 should be of type double")
   if ((pi0<=0)|(pi0>1)) stop("pi0 should be in ]0,1]")

   n = nrow(erpdta)                  # sets the value for n
   T = ncol(erpdta)                  # sets the value for T 

   pdesign = solve(t(design)%*%design)%*%t(design)
   Proj = design%*%pdesign
   Proj0 = design0%*%solve(t(design0)%*%design0)%*%t(design0)
   rdf1 = nrow(design)-ncol(design)
   rdf0 = nrow(design0)-ncol(design0)

   beta = (pdesign%*%erpdta)[idsignal,]
   if (length(idsignal)==1) beta = matrix(beta,nrow=1)

   res1 = erpdta-Proj%*%erpdta
   scer1 = as.vector(t(rep(1,n))%*%res1^2)
   res0 = erpdta-Proj0%*%erpdta
   scer0 = as.vector(t(rep(1,n))%*%res0^2)
   fstat = ((scer0-scer1)/(rdf0-rdf1))/(scer1/rdf1)
   pval = pf(fstat,df1=rdf0-rdf1,df2=rdf1,lower.tail=FALSE)

   if (is.null(pi0)) pi0 = pval.estimate.eta0(pval,diagnostic.plot=FALSE)

   correctedpval = pi0*p.adjust(pval,method=method)
   significant = which(correctedpval<=alpha)
   return(list(pval=pval,correctedpval=correctedpval,significant=significant,pi0=pi0,test=fstat,df1=rdf1,df0=rdf0,signal=beta,r2=(1-1/(1+fstat*((rdf0-rdf1)/rdf1)))))
}
