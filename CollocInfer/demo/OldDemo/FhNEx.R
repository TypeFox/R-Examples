##########################################
##       FitzHugh-Nagumo Example       ###
##########################################

library("CollocInfer")

#########################################################
####   Data generation from a solution of the ODE #######
#########################################################


data(FhNdata)

#  define the functions using function make.fhn

#  the functions that make.fhn defines:
#        fn       =  fhn.fun,
#        fn.ode  = fhn.fun.ode,
#        dfdx    = fhn.dfdx,
#        dfdp    = fhn.dfdp,
#        d2fdx2  = fhn.d2fdx2,
#        d2fdxdp = fhn.d2fdxdp

FhN = make.fhn()



# We can also remove the 'R' part of the data (ie, assume we did not measure
# it) and set the corresponding coefficients to zero

V.FhNdata = FhNdata
V.FhNdata[,2] = NA

# Now we want to set up a basis expansion; this makes use of
# routines in the fda library

knots  = seq(0,20,0.5)
norder = 3
nbasis = length(knots) + norder - 2
range  = c(0,20)

FhNbasis = create.bspline.basis(range=range,nbasis=nbasis,
	norder=norder,breaks=knots)
	
# We'll start off by creating a smooth of each state variable to get initial
# values for the coefficients.  

fdnames=list(NULL,c('V','R'),NULL) 
bfdPar = fdPar(FhNbasis,lambda=1,int2Lfd(1))
DEfd0 = smooth.basis(FhNtimes,FhNdata,bfdPar,fdnames=fdnames)$fd

#  plot the solution

par(ask=FALSE)
plotfit.fd(FhNdata, FhNtimes, DEfd0)

#  extract the coefficients and assign variable names

coefs0 = DEfd0$coefs
colnames(coefs0) = FhNvarnames

#  In the case of only V being measured, we can get coefficients for R,
#  so we set these to zero as starting values

V.coefs0 = coefs0
V.coefs0[,2] = 0

#  Set a value for lambda

lambda = 1000

# Optimize the coefficients  in coefs0, that is, run the lowest level or
# inner optimization loop, maximizing function J(c|theta,lambda)

Ires1	= Smooth.LS(FhN,FhNdata,FhNtimes,FhNpars,coefs0,FhNbasis,lambda,
                  in.meth='nlminb')

# We have some handy plotting routines 

plotout1 = CollocInferPlots(Ires1$coefs,FhNpars,Ires1$lik,Ires1$proc,FhNtimes,FhNdata)


# Now we can do the profiling, running the outer loop to optimize
#  parameter values.  You may wish to turn off output buffering in the
#  Misc menu in the command window at this point.

Ores1 = Profile.LS(FhN,FhNdata,FhNtimes,FhNpars,coefs=Ires1$coefs,
                   basisvals=FhNbasis,lambda=lambda)

# And look at the result

Ores1$pars

plotout2 = CollocInferPlots(Ores1$coefs,Ores1$pars,Ores1$lik,Ores1$proc,FhNtimes,FhNdata)

#  Repeat these analyses for the data with only V measured

Ores1.V = Profile.LS(FhN,V.FhNdata,FhNtimes,FhNpars,coefs=V.coefs0,
                     basisvals=FhNbasis,lambda=lambda)

Ores1.V$pars

plotout2.V = CollocInferPlots(Ores1.V$coefs,Ores1.V$pars,Ores1.V$lik,Ores1.V$proc,FhNtimes,V.FhNdata)



# If we only coded the FitzHugh Nagumo function with no derivatives we would just
# have make.fhn()$fn, in which case we default to estimating the derivatives by
# finite differencing.  The first argument provides access only to the 
# function evaluating the right hand side.

fhn.fun <- function(times, y, p, more) {
        r = y
        r[, "V"] = p["c"] * (y[, "V"] - y[, "V"]^3/3 + y[, "R"])
        r[, "R"] = -(y[, "V"] - p["a"] + p["b"] * y[, "R"])/p["c"]
        return(r)
    }


Ires1a = Smooth.LS(fhn.fun,FhNdata,FhNtimes,FhNpars,coefs0,FhNbasis,lambda,
                   in.meth='nlminb')

# Now we can do the profiling

Ores1a = Profile.LS(fhn.fun,FhNdata,FhNtimes,FhNpars,coefs=Ires1a$coefs,
                    basisvals=FhNbasis,lambda=lambda,out.meth='nlminb')
  
Ores1a$pars

plotout1a = CollocInferPlots(Ores1a$coefs,Ores1a$pars,Ores1a$lik,Ores1a$proc,FhNtimes,V.FhNdata)
  
#### Option 2:  set-up functions for proc and lik objects and then call
# optimizers.  This is a longer form of the analysis, in which the
#  lik and proc objects defining the first and second terms of 
#  J(c|theta,lambda) are specified.  This is done in function LS.setup.
#  Then we can optimize initial parameter values using function
#  ParsMatchOpt

profile.obj = LS.setup(pars=FhNpars,fn=FhN,lambda=lambda,times=FhNtimes,
                       coefs=coefs1,basisvals=FhNbasis)
lik  = profile.obj$lik
proc = profile.obj$proc

# Now we can get initial parameter estimates from "gradient matching"
#  using function ParsMatchOpt

Pres  = ParsMatchOpt(FhNpars,coefs1,proc)
pars1 = Pres$pars

pars1

# Smoothing can be done more specifically with 'inneropt'

Ires2 = inneropt(FhNdata,times=FhNtimes,pars1,coefs1,lik,proc)

# And we can also do profiling

Ores2 = outeropt(FhNdata,FhNtimes,pars1,coefs1,lik,proc)

Ores2$pars

DEfd2 = fd(Ores2$coefs,FhNbasis,fdnames)

plotfit.fd(FhNdata, FhNtimes, DEfd2)

# Option 3:  Set everything up manually.  This is the really long way to
# do the analysis, in which the error sum of squares fit is defined
# manually for both terms of J

## lik object

lik2 = make.SSElik()                      # we are using squared error
lik2$bvals = eval.basis(FhNtimes,FhNbasis)  # values of the basis at times
lik2$more = make.id()                     # identity transformation
lik2$more$weights = c(1,0)                # only use data for V

## proc object

qpts = knots[1:(length(knots)-1)]+diff(knots)/2  # Collocation points at 
                                                 # midpointsbetween knots

proc2 = make.SSEproc()                    # we are using squared error
proc2$bvals = list(bvals = eval.basis(qpts,FhNbasis),    # values and derivative of
                  dbvals = eval.basis(qpts,FhNbasis,1))  # basis at collocation points
proc2$more = make.fhn()                   # FitzHugh-Nagumo right hand side
proc2$more$names = FhNvarnames            # State variable names
proc2$more$parnames = FhNparnames         # Parameter names
proc2$more$qpts = qpts                    # Collocation points
proc2$more$weights  = c(1000,1000)        # Weight relative to observations of
                                          # each collocation point.
                                          
# We can also specify optimization methods

in.meth = 'nlminb'            # Inner Optimization
control.in = list(rel.tol=1e-12,iter.max=1000,eval.max=2000,trace=0)  # Optimization control

out.meth = 'optim'            # Outer Optimization
control.out = list(trace=6,maxit=100,reltol=1e-8,meth='BFGS') # Optimization control
                                                                                   
# We can now fix the smooth or 'V' and try and choose 'R' to make the differential
# equation match as well as possible.

Vres0 = FitMatchOpt(coefs=V.coefs0,which=2,pars=FhNpars,proc2)
V.coefs1 = Vres0$coefs

# And we can call the same inner optimization as above, this time with our
# own control parameters and method

Ires3 = inneropt(V.FhNdata,FhNtimes,FhNpars,V.coefs1,lik2,proc2,
                 in.meth=in.meth,control.in=control.in)

# And we can also do profiling with specified controls

Ores3.V = outeropt(V.FhNdata,FhNtimes,FhNpars,coefs=V.coefs1,lik=lik2,proc=proc2,
                   in.meth=in.meth,out.meth=out.meth,
                   control.in=control.in,control.out=control.out)

Ores3.V$pars

DEfd3.V = fd(Ores3.V$coefs,FhNbasis,fdnames)

plot(DEfd3.V['V'])
points(FhNtimes, V.FhNdata[,1])

plot(DEfd3.V['R'])

# We can also use finite differencing to calculate derivatives of the right hand side of the
# FitzHugh-Nagumo equations, this is defined by modifying the proc object.

proc2a = make.SSEproc()                    # we are using squared error
proc2a$bvals = list(bvals = eval.basis(qpts,FhNbasis),    # values and derivative of
                  dbvals = eval.basis(qpts,FhNbasis,1))  # basis at collocation points
proc2a$more = make.findif.ode()                   # FitzHugh-Nagumo right hand side
proc2a$more$names = FhNvarnames            # State variable names
proc2a$more$parnames = FhNparnames         # Parameter names
proc2a$more$qpts = qpts                    # Collocation points
proc2a$more$weights  = c(1000,1000)        # Weight relative to observations of
                                           # each collocation point.
proc2a$more$more = list(fn = make.fhn()$fn, eps=1e-6) # Tell findif that the function
                                                      # to difference is fhn and the
                                                      # difference stepsize


Ires3a = inneropt(V.FhNdata,FhNtimes,FhNpars,Ores3.V$coefs,lik2,proc2a,
  in.meth=in.meth,control.in=control.in)

# And we can also do profiling with specified controls

Ores3a = outeropt(V.FhNdata, FhNtimes, FhNpars, coefs=Ores3.V$coefs, 
                  lik=lik2, proc=proc2a, 
                  in.meth=in.meth, out.meth=out.meth, 
                  control.in=control.in, control.out=control.out)
