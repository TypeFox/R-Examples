kktc <- function(par, fval, ngr, nHes, nbm, maxfn=FALSE, control=list()) {
# Provide a check on Kuhn-Karush-Tucker conditions based on quantities
# already computed. Some of these used only for reporting.
#
# Input:
#  par = a single vector of starting values
#  fval = objective function value
#  ngr = gradient evaluated at parameters par
#  nHes = Hessian matrix evaluated at the parameters par
#  nbm = number of active bounds and masks from gHgenb (gHgen returns 0)
#  maxfn = logical TRUE if we want to maximize the function. Default FALSE.
#  control = list of controls, currently, 
#            kkttol=1e-3, kkt2tol=1e-6, ktrace=FALSE
#
# Output: A list of four elements, namely,
#  gmax = max abs gradient element
#  evratio = ratio of smallest to largest eigenvalue of estimated Hessian
#  kkt1 = logical flag: TRUE if gradient KKT test is satisfied to tolerances
#         otherwise FALSE. Note that decisions are sensitive to scaling and tolerances
#  kkt2 = logical flag: TRUE if Hessian KKT test is satisfied to tolerances
#         otherwise FALSE. Note that decisions are sensitive to scaling and tolerances
# 
# NOTE: Does not do projections for bounds and masks! 
# NOTE: Assumes we are MINIMIZING within kktc, even if outer problem is a maximization.
#################################################################
## Post-processing -- Kuhn Karush Tucker conditions
#  Ref. pg 77, Gill, Murray and Wright (1981) Practical Optimization, Academic Press
#   print(control)
   if (is.null(control$ktrace)) trace<-FALSE else trace<-control$ktrace
   if (is.null(control$kkttol)) kkttol<-1e-3 else kkttol<-control$kkttol
   if (is.null(control$kkt2tol)) kkt2tol<-1e-6 else kkt2tol<-control$kkt2tol
   if (trace) { cat("kkttol=",kkttol,"   kkt2tol=",kkt2tol,"  trace=",trace,"\n") }
   if (trace) { cat("fval =",fval,"\n") }
   npar<-dim(nHes)[1]
   if (trace) { 
      cat("KKT condition testing\n") 
      cat("Number of parameters =",npar,"\n")
   }
   kkt1<-NA
   kkt2<-NA
   # test gradient
   gmax<-max(abs(ngr)) # need not worry about sign for maximizing
   if (trace) {
      cat("max abs gradient element =",gmax,"  test tol = ",kkttol*(1.0+abs(fval)),"\n")
   }
   kkt1<-(gmax <= kkttol*(1.0+abs(fval)) ) # ?? Is this sensible?
   if (trace) {cat("KKT1 result = ",kkt1,"\n") }
   if (maxfn) {
      nHes<- -nHes
      if (trace) cat("Maximizing: use negative Hessian\n")
   }
   hev<- try(eigen(nHes)$values, silent=TRUE) # 091215 use try in case of trouble, 
                                              # 20100711 silent
   if (trace) {
      cat("Hessian eigenvalues:\n")
      print(hev)
   }
   if (class(hev) != "try-error") {
       # now look at Hessian
       negeig<-(hev[npar] <= (-1)*kkt2tol*(1.0+abs(fval))) # 20100711 kkt2tol
       evratio<-hev[npar-nbm]/hev[1]
       # If non-positive definite, then there will be zero eigenvalues (from the projection)
       # in the place of the "last" eigenvalue and we'll have singularity.
       # WARNING: Could have a weak minimum if semi-definite.
       kkt2<- (evratio > kkt2tol) && (! negeig) 
       if (trace) {
          cat("KKT2 result = ",kkt2,"\n") 
       }
       ans<-list(gmax,evratio,kkt1,kkt2,hev)
       names(ans)<-c("gmax","evratio","kkt1","kkt2","hev")
       return(ans)
   } else {
       warning("Eigenvalue failure")
       if(trace) cat("Eigenvalue calculation has failed!\n") # JN 111207 added \n
   } # end kkt test
} ## end of kktc
