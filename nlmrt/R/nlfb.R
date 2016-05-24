nlfb <-function(start, resfn, jacfn=NULL, trace=FALSE, 
            lower=-Inf, upper=Inf, maskidx=NULL, control=list(), ...){
#
#  A simplified and hopefully robust alternative to finding the 
#  nonlinear least squares minimizer that causes 'formula' to 
#  give a minimal residual sum of squares. 
#
#  nls.mn is particularly intended to allow for the resolution of 
#  very ill-conditioned or else near zero-residual problems for 
#  which the regular nls() function is ill-suited. It may also be
#  a useful pre-processor for nls().
#  
#  J C Nash  2012-3-4   nashjc _at_  uottawa.ca
#
#  start MUST be a vector where all the elements are named:
#     e.g., start=c(b1=200, b2=50, b3=0.3)
#  trace -- TRUE for console output 
#  lower is a vector of lower bounds
#  upper is a vector of upper bounds
#  maskidx is an index vector of positions of parameters that are fixed.
#  control is a list of control parameters. These are:
#     ...
#  
#  ... will need to contain data for other variables that appear in
#  the functions
#
#  This variant uses a qr solution without forming the sum of squares and cross
#  products t(J)%*%J 
# 
# Function to display SS and point
showprms<-function(SS, pnum){
    pnames<-names(pnum)
    npar<-length(pnum)
    cat("lamda:",lamda," SS=",SS," at")
    for (i in 1:npar){
       cat(" ",pnames[i],"=",pnum[i])
    }
    cat(" ",feval,"/",jeval)
    cat("\n")
}
#########
# ensure params in vector
pnames<-names(start)
start<-as.numeric(start)
names(start)<-pnames # ?? needed
# bounds
npar<-length(start) # number of parameters
if (length(lower)==1) lower<-rep(lower,npar)
if (length(upper)==1) upper<-rep(upper,npar) 
# ?? more tests on bounds
if (length(lower)!=npar) stop("Wrong length: lower")
if (length(upper)!=npar) stop("Wrong length: upper")
if (any(start<lower) || any(start>upper)) stop("Infeasible start")
if (trace) {
   cat("lower:")
   print(lower)
   cat("upper:")
   print(upper)
}
# Should make this more informative??
# controls
   ctrl<-list(
    watch=FALSE, # monitor progress
    phi=1, # the phi parameter
    lamda=0.0001, # lamda (spelled lamda in JNWMS)
    offset=100, # to determine if paramters changed
    laminc=10,
    lamdec=4, # use decreased_lamda<-lamda*lamdec/laminc
    femax=10000,
    jemax=5000,
    rofftest = TRUE,
    smallsstest = TRUE,
    ndstep=1e-7 # numerical jacobian step
   )
   epstol<-(.Machine$double.eps)*ctrl$offset
   ncontrol <- names(control)
   nctrl <- names(ctrl)
   for (onename in ncontrol) {
      if (!(onename %in% nctrl)) {
         if (trace) cat("control ",onename," is not in default set\n")
      }
      ctrl[onename]<-control[onename]
   }
   # print(ctrl)
   phiroot<-sqrt(ctrl$phi)
   lamda<-ctrl$lamda
   offset<-ctrl$offset
   laminc<-ctrl$laminc
   lamdec<-ctrl$lamdec # save typing
   watch<-ctrl$watch
   femax<-ctrl$femax
   jemax<-ctrl$jemax
# Then see which ones are parameters (get their positions in the set xx
    pnum<-start # may simplify later??
    pnames<-names(pnum)
    bdmsk<-rep(1,npar) # set all params free for now
    if (length(maskidx)>0 && trace) {
       cat("The following parameters are masked:")
       print(maskidx)
    }
    bdmsk[maskidx]<-0 # fixed parameters
# 20120607 -- put in if needed
    if (is.null(jacfn)){
       if (trace) cat("Using default jacobian approximation\n")
       numjac<-TRUE
       myjac<-function(pars, rfn=resfn, bdmsk=bdmsk, resbest=resbest, ...){
         npar<-length(pars)
         jacmat<-matrix(0, ncol=npar, nrow=length(resbest))
         tpars<-pars
         for (j in 1:npar){
            if (bdmsk[j]!=0) {
               step<-ctrl$ndstep*(abs(pars[j])+ctrl$ndstep)
##               cat("step ",j," = ",step,"\n")
               tpars[j]<-tpars[j]+step
               jacmat[,j]<-(rfn(tpars,...)-resbest)/step
               tpars[j]<-pars[j]
##               cat("Column ",j," of jacobian :\n")
##               print(jacmat[,j])
            }
         }
         jacmat
       } # end myjac numeric
    } else { 
       numjac<-FALSE
    }
    resbest<-resfn(pnum, ...)
#    cat("resbest:")
#    print(resbest)
    ssbest<-crossprod(resbest)
    ss0 <- ssbest # reference value for rofftest
    ssminval <- ssbest*epstol^4
    if (trace) cat("ssminval =",ssminval,"\n")
    feval<-1
    pbest<-pnum
    feval<-1 # number of function evaluations
    jeval<-0 # number of Jacobian evaluations
    if (trace) {
       cat("Start:")
       showprms(ssbest,pnum)
       if (watch) tmp<-readline("Continue")
    }
    if (length(maskidx) == npar) stop("All parameters are masked") # Should we return?
    ssquares<-.Machine$double.xmax # make it big
    newjac<-TRUE # set control for computing Jacobian
    eqcount<-0
    roffstop <- FALSE
    smallstop <- FALSE
    while ((! roffstop) && (eqcount < npar) 
             && (feval <= femax) && (jeval <= jemax)
             && (! smallstop) ) {
       if (newjac) {
          bdmsk<-rep(1,npar)
          bdmsk[maskidx]<-0
          bdmsk[which(pnum-lower<epstol*(abs(lower)+epstol))]<- -3 
          bdmsk[which(upper-pnum<epstol*(abs(upper)+epstol))]<- -1
          if (trace && watch) {
            cat("bdmsk:")
            print(bdmsk)
          }
          if (numjac) Jac<-myjac(pbest, rfn=resfn, bdmsk=bdmsk, resbest=resbest, ...)
          else Jac<-jacfn(pbest, ...)
          jeval<-jeval+1 # count Jacobians
          if (any(is.na(Jac))) stop("NaN in Jacobian")
          JTJ<-crossprod(Jac)
          gjty<-t(Jac)%*%resbest # raw gradient
          for (i in 1:npar){
             bmi<-bdmsk[i]
             if (bmi==0) {
                gjty[i]<-0 # masked
                Jac[,i]<-0
             }
             if (bmi<0) {
                if((2+bmi)*gjty[i] > 0) { # free parameter
                   bdmsk[i]<-1
                   if (trace) cat("freeing parameter ",i," now at ",pnum[i],"\n")
                } else {
                   gjty[i]<-0 # active bound
                   Jac[,i]<-0
                   if (trace) cat("active bound ",i," at ",pnum[i],"\n") 
                }
             } # bmi
          } # end for loop
          if (npar == 1) dee <- diag(as.matrix(sqrt(diag(crossprod(Jac)))))
          else dee <- diag(sqrt(diag(crossprod(Jac))))  # to append to Jacobian
       } # end newjac
       lamroot<-sqrt(lamda)
       JJ<-rbind(Jac,lamroot*dee, lamroot*phiroot*diag(npar)) # build the matrix
       if (trace && watch) {
         cat("JJ\n")
         print(JJ)
       }
       JQR<-qr(JJ)# ??try
       rplus<-c(resbest, rep(0,2*npar))
       roff <- max(abs(as.numeric(crossprod(qr.Q(JQR), rplus))))/ss0
       if (trace) cat("roff =", roff,"  converged = ",(roff <= sqrt(epstol)),"\n")
       if (ctrl$rofftest && (roff <= sqrt(epstol))) roffstop <- TRUE
#        tmp <- readline('cont')
       delta<-try(qr.coef(JQR,-rplus)) # Note appended rows of y)
       if (class(delta)=="try-error") {
          if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
          lamda<-laminc*lamda
          newjac<-FALSE # increasing lamda -- don't re-evaluate
          if (trace) cat(" Equation solve failure\n")
          feval<-feval+1 # count as a function evaluation to force stop
       } else { # solution OK
          gproj<-crossprod(delta,gjty)
          gangle <- gproj/sqrt(crossprod(gjty) * crossprod(delta))
          gangle <- 180 * acos(sign(gangle)*min(1, abs(gangle)))/pi
          if (trace) cat("gradient projection = ",gproj," g-delta-angle=",gangle,"\n")
          if (is.na(gproj) || (gproj >= 0) ) { # uphill direction -- should NOT be possible
            if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
            lamda<-laminc*lamda
            newjac<-FALSE # increasing lamda -- don't re-evaluate
            if (trace) cat(" Uphill search direction\n")
          } else { # downhill
            delta[maskidx]<-0
            delta<-as.numeric(delta)
            if (trace && watch) {
              cat("delta:")
              print(delta)
            }
            step<-rep(1,npar)
            for (i in 1:npar){
                bd<-bdmsk[i]
                da<-delta[i]
#                if (trace) cat(i," bdmsk=",bd,"  delta=",da,"\n")
                if (bd==0 || ((bd==-3) && (da<0)) ||((bd==-1)&& (da>0))) {
                  delta[i]<-0
                } else {
                  if (delta[i]>0) step[i]<-(upper[i]-pbest[i])/delta[i]
                  if (delta[i]<0) step[i]<-(lower[i]-pbest[i])/delta[i] # positive
                }
            }
            stepsize<-min(1,step[which(delta!=0)])
            if (trace) cat("Stepsize=",stepsize,"\n")
            if (stepsize<.Machine$double.eps) {
              if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
              lamda<-laminc*lamda
              newjac<-FALSE # increasing lamda -- don't re-evaluate
              if (trace) cat(" Step too small or not possible\n")
            } else { # continue
            pnum<-pbest+stepsize*delta # adjust (note POSITIVE here, but not in nlsmn0
            names(pnum)<-pnames # NOT inherited through %*% !!!
            eqcount<-length(which((offset+pbest)==(offset+pnum)))
            if (eqcount<npar) {
#              for (i in 1:npar){
#                joe<-paste(pnames[[i]],"<-",pnum[[i]])
#                eval(parse(text=joe))
#              }
              feval<-feval+1 # count evaluations
              resid<-resfn(pnum, ...)
              ssquares<-as.numeric(crossprod(resid))
              if (is.na(ssquares)) ssquares<-.Machine$double.xmax
              if (ssquares>=ssbest) {
                if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
                lamda<-laminc*lamda
                newjac<-FALSE # increasing lamda -- don't re-evaluate
                if(trace) showprms(ssquares, pnum)
              } else {
                lamda<-lamdec*lamda/laminc
                if (trace) {
                  cat("<<")
                  showprms(ssquares, pnum)
                }
                ssbest<-ssquares
                if (ctrl$smallsstest) { smallstop<-(ssbest <= ssminval) }
                resbest<-resid
                pbest<-pnum
                newjac<-TRUE
              } # reduced sumsquares
            } else {# end if equcount
  	            if (trace) cat("No parameter change\n")
            }
            } # end stepsize not too small
          } # end downhill
        } # solution OK
        if (watch) tmp<-readline("Cycle")
     } # end main while loop 
    pnum<-as.vector(pnum)
    names(pnum) <- pnames
    result <- list(resid = resbest, jacobian = Jac, feval = feval, 
        jeval = jeval, coefficients = pnum, ssquares = ssbest, lower=lower, upper=upper, maskidx=maskidx)
    class(result) <- "nlmrt"
    result
}

