nlsmn0b <-function(formula, start, trace=FALSE, data=NULL, 
            lower=-Inf, upper=Inf, masked=NULL, control=list(), ...){
#
#  A simplified and hopefully robust alternative to finding the 
#  nonlinear least squares minimizer that causes 'formula' to 
#  give a minimal residual sum of squares. 
#
#  nlsmn0b is particularly intended to allow for the resolution of 
#  very ill-conditioned or else near zero-residual problems for 
#  which the regular nls() function is ill-suited. It may also be
#  a useful pre-processor for nls().
#  
#  J C Nash  2012-3-4   nashjc _at_  uottawa.ca
#
#  formula looks like "y~b1/(1+b2*exp(-b3*T))"
#  start MUST be a vector where all the elements are named:
#     e.g., start=c(b1=200, b2=50, b3=0.3)
#  trace -- TRUE for console output 
#  data is a data frame containing data for variables used in the formula
#     that are NOT the parameters. This data may also be defined in the 
#     parent frame i.e., "global" to this function
#  lower is a vector of lower bounds
#  upper is a vector of upper bounds
#  masked is a character vector of names of parameters that are fixed.
#  control is a list of control parameters. These are:
#     ...
#  
#  ... will need to contain data for other variables that appear in
#  the formula and are defined in a parent frame (Not sure how needed??)
#
#  This variant uses a traditional solve() approach
# 
# to display SS and point
showpoint<-function(SS, pnum){
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
# get data from data frame if exists
if (! is.null(data)){
    for (dfn in names(data)) {
       cmd<-paste(dfn,"<-data$",dfn,"")
       eval(parse(text=cmd))
    }
} else stop("'data' must be a list or an environment")
# ensure params in vector
pnames<-names(start)
start<-as.numeric(start)
names(start)<-pnames
# bounds
npar<-length(start) # number of parameters
if (length(lower)==1) lower<-rep(lower,npar)
if (length(upper)==1) upper<-rep(upper,npar) 
# ?? more tests on bounds
if (length(lower)!=npar) stop("Wrong length: lower")
if (length(upper)!=npar) stop("Wrong length: upper")
if (any(start<lower) || any(start>upper)) stop("Infeasible start")
if (trace) {
   cat("formula: ")
   print(formula)
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
    jemax=5000,
    femax=10000
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
   phiroot<-sqrt(ctrl$phi)
   lamda<-ctrl$lamda
   offset<-ctrl$offset
   laminc<-ctrl$laminc
   lamdec<-ctrl$lamdec # save typing
   watch<-ctrl$watch
   jemax<-ctrl$jemax
   femax<-ctrl$femax
#  First get all the variable names:
    vn <- all.vars(parse(text=formula))
# Then see which ones are parameters (get their positions in the set xx
    pnum<-start # may simplify later??
    pnames<-names(pnum)
    bdmsk<-rep(1,npar) # set all params free for now
    mskdx<-which(pnames %in% masked) # NOTE: %in% not == or order gives trouble
    bdmsk[mskdx]<-0 # fixed parameters
    if (trace) {
      parpos  <- match(pnames, vn)
      datvar<-vn[-parpos] # NOT the parameters
      for (dvn in datvar){
         cat("Data variable ",dvn,":")
         print(eval(parse(text=dvn)))
      }
    }
   if (is.character(formula)){
       es<-formula
    } else {
       tstr<-as.character(formula) # note ordering of terms!
       es<-paste(tstr[[2]],"~",tstr[[3]],'')
    }
# Now separate the sides
    parts<-strsplit(as.character(es), "~")[[1]]
    if (length(parts)!=2) stop("Model expression is incorrect!")
    lhs<-parts[1]
    rhs<-parts[2]
# And build the residual at the parameters
    resexp<-paste(rhs,"-",lhs, collapse=" ")
    for (i in 1:npar){ # put parameters in separate variables
       joe<-paste(pnames[[i]],"<-",pnum[[i]])
       eval(parse(text=joe))
    }
    gradexp<-deriv(parse(text=resexp), names(start)) # gradient expression
    resbest<-with(data, eval(parse(text=resexp)))
    ssbest<-crossprod(resbest)
    feval<-1
    pbest<-pnum
    feval<-1 # number of function evaluations
    jeval<-0 # number of Jacobian evaluations
    if (trace) {
       cat("Start:")
       showpoint(ssbest,pnum)
    }
    ssquares<-.Machine$double.xmax # make it big
    newjac<-TRUE # set control for computing Jacobian
       bdmsk<-rep(1,npar)
    eqcount<-0
    while ((eqcount < npar) && (feval<=femax) && (jeval<=jemax) ) {
       if (trace && watch) {
          cat("bdmsk:")
          print(bdmsk)
       }
       if (newjac) {
          bdmsk<-rep(1,npar) # set the constraint indicators
          bdmsk[mskdx]<-0
          bdmsk[which(pnum-lower<epstol*(abs(lower)+epstol))]<- -3 
          bdmsk[which(upper-pnum<epstol*(abs(upper)+epstol))]<- -1
          J0<-with(data, eval(gradexp))
          Jac<-attr(J0,"gradient")
          jeval<-jeval+1 # count Jacobians
          if (any(is.na(Jac))) stop("NaN in Jacobian")
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
                   if (trace) print(Jac)
                }
             } # end for loop
          } # end for loop
          JTJ<-crossprod(Jac)
          gjty[mskdx]<-0 # only masked ones
          dde<-diag(diag(JTJ))+phiroot*diag(npar) # to append to Jacobian
          gjty<-as.numeric(gjty)
          if (trace && watch) {
             cat("gjty:")
             print(gjty)
             cat("JTJ\n")
             print(JTJ)
             cat("dde")
             print(dde)
             cat("bdmsk:")
             print(bdmsk)
          }
       } # end newjac
       Happ<-JTJ+lamda*dde # build the crossprods matrix
       if (trace && watch) {
         cat("Happ\n")
         print(Happ)
       }
#       tmp<-readline("more")
       delta<-try(solve(Happ,-gjty))
       if (class(delta)=="try-error") {
          if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
          lamda<-laminc*lamda
          newjac<-FALSE # increasing lamda -- don't re-evaluate
          if (trace) cat(" Equation solve failure\n")
          print(Happ)
          print(gjty)
          feval<-feval+1 # count as a function evaluation
          print(pnum)
          tmp<-readline("continue")
       } else { # solution OK
          gproj<-crossprod(delta,gjty)
          gangle<-gproj/sqrt(crossprod(gjty)*crossprod(delta))
          if (trace) cat("gradient projection0 = ",gproj," gangle=",gangle,"\n")
          if (gproj >= 0) { # uphill direction (??test)
            if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
            lamda<-laminc*lamda
            newjac<-FALSE # increasing lamda -- don't re-evaluate
            if (trace) cat(" Uphill search direction\n")
          } else { # downhill
            delta[mskdx]<-0
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
#            if (trace && watch) {
#              cat("step:")
#              print(step)
#            }
            stepsize<-min(1,step[which(delta!=0)])
            if (trace) cat("Stepsize=",stepsize,"\n")
            if (stepsize<.Machine$double.eps) {
              if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
              lamda<-laminc*lamda
              newjac<-FALSE # increasing lamda -- don't re-evaluate
              if (trace) cat(" Stepsize too small\n")
            } else { # continue
            pnum<-pbest+stepsize*delta # adjust (note POSITIVE here, but not in nlsmn0
            names(pnum)<-pnames # NOT inherited through %*% !!!
            eqcount<-length(which((offset+pbest)==(offset+pnum)))
            if (eqcount<npar) {
              for (i in 1:npar){
                joe<-paste(pnames[[i]],"<-",pnum[[i]])
                eval(parse(text=joe))
              }
              feval<-feval+1 # count evaluations
              resid<-with(data, eval(parse(text=resexp)))
              print(pnum)
              ssquares<-as.numeric(crossprod(resid))
              print(ssquares)
              if (is.na(ssquares) ) ssquares<-.Machine$double.xmax
              if (ssquares>=ssbest) {
                if (lamda<1000*.Machine$double.eps) lamda<-1000*.Machine$double.eps
                lamda<-laminc*lamda
                newjac<-FALSE # increasing lamda -- don't re-evaluate
                if(trace) showpoint(ssquares, pnum)
              } else {
                lamda<-lamdec*lamda/laminc
                if (trace) {
                  cat("<<")
                  showpoint(ssquares, pnum)
                }
                ssbest<-ssquares
                resbest<-resid
                pbest<-pnum
                newjac<-TRUE
              } # reduced sumsquares
#              if (watch) tmp<-readline()
            } else {# end if equcount
              if (trace) cat("No parameter change\n")
            }
            } # end nonzero step
          } # end downhill
        } # solution OK
        if (watch) tmp<-readline()
     } # end main while loop 
    pnum<-as.vector(pnum)
    result<-list(resid=resbest, jacobian=Jac, feval=feval, jeval=jeval, coeffs=pnum, ssquares=ssbest)
}
