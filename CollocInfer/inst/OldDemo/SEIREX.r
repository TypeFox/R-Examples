##########
# This demonstration illustrates the use of log transformations with
# CollocInfer. We use the SEIR equations with a seasonally varying
# infection rate for this purpose.
#
# The equations are given as
#
# Sdot = mu - beta(t)*S*(I+i) - nu*S        (Susceptibles)
# Edot = beta(t)*S*(I+i) - (sigma+nu)*E     (Exposed)
# Idot = sigma*E - (gamma+nu)*I              (Infectious)
#
#
# Here beta(t) - the infection rate - is parameterized by a sinusoidal function
# plus a constant.
#
# Other parameters are
# i - a visiting process
# nu - death rate
# sigma - the rate of movement from Exposed to Infectious.
# gamma - the rate of recovery from infection.
#
# It is generally more stable to solve the equations for the log states rather
# than the states themselves. CollocInfer contains a number of useful tools
# that let you make this transition without needing to re-code all of your
# differential equations.

library('CollocInfer')
            
            
#### labels for the state process

SEIRnames = c('S','E','I')

#### define the right side evaluation function

SEIRodefn = make.SEIR()$fn.ode


N = 5e6

####  define true parameter values and names

pars = c(7/5,0.02*N,0,0.02,365/8,365/5,17*(365/5)/N, 0, 0.08*17*(365/5)/N)
parnames = c('rho','mu','i','nu','sigma','gamma','b0','b1','b2')
names(pars) = parnames


#### A couple of functions to define the infection rate

beta.fun = function(t,p,more){
    return( p['b0'] + p['b1']*sin(2*pi*t) + p['b2']*cos(2*pi*t) )
}

beta.dfdp = function(t,p,more){
    dfdp =  cbind(rep(1,length(t)), sin(2*pi*t), cos(2*pi*t)) 
    colnames(dfdp) = c('b0','b1','b2')
    return(dfdp)
}

betamore = list(beta.fun=beta.fun,
                beta.dfdp=beta.dfdp,
                beta.ind=c('b0','b1','b2'))


#### Now create ODE paths


y0 = N*c(0.9,0,0.1)


y0 =N*c(.056,.00014,.000085)
names(y0) = SEIRnames
times = seq(0,3,1/52)

y = lsoda(y0,times=times,func=make.SEIR()$ode.fn,parms=list(p=pars,more=betamore))


data = y[,4,drop=FALSE]*7/5



# and we're going to need some basis values

rr = range(times)
knots = seq(rr[1],rr[2],1/52)
norder = 3
nbasis = length(knots)+norder-2

bbasis = create.bspline.basis(range=rr,norder=norder,nbasis=nbasis,breaks=knots)
qpts = 0.5*(knots[1:(length(knots)-1)] + knots[2:length(knots)])

bvals.obs = Matrix(eval.basis(times,bbasis),sparse=TRUE)

bvals = list(bvals = Matrix(eval.basis(qpts,bbasis),sparse=TRUE),
            dbvals = Matrix(eval.basis(qpts,bbasis,1),sparse=TRUE))


### Proc object, the first of these is just the standard squared error deviation
# from the right hand side of a differential equation

sproc = make.SSEproc()
sproc$bvals = bvals
sproc$more = make.SEIR()
sproc$more$more = betamore
sproc$more$qpts = qpts
sproc$more$weights = matrix(1,length(qpts),3)%*%diag(c(1e2,1e0,1e0))
sproc$more$names = SEIRnames
sproc$more$parnames = parnames


### However, ODEs are often much more numerically stable if represented on a
# log scale. The make.logtrans() function will take the right hand side functions
# and derivatives defined for any differential equation and provide the equivalent
# system for the log state. Note that this does affect the way  you represent
# your state when considering the observations. 

lsproc = sproc
lsproc$more = make.logtrans()
lsproc$more$more = make.SEIR()
lsproc$more$more$more = betamore
lsproc$more$qpts = qpts
lsproc$more$weights = matrix(1,length(qpts),3)%*%diag(c(1e2,1e0,1e0))
lsproc$more$names = SEIRnames
lsproc$more$parnames = parnames

### Lik objects, this is the standard squared error. 

slik = make.SSElik()
slik$bvals = eval.basis(times,bbasis)
slik$more = make.genlin()
slik$more$weights = rep(1,length(times))
slik$more$more = list(mat = matrix(0,1,3),sub=matrix(c(1,3,1),1,3))
slik$more$names = SEIRnames
slik$more$parnames = parnames

# Log transform transformation. For this, we note that we have represented
# the trajectory on the log scale and will need to transform back. 

lslik = make.logstate.lik()
lslik$bvals = slik$bvals
lslik$more$weights = slik$more$weights
lslik$more = slik
lslik$more$parnames = parnames


# Numerically things work much better on the log scale

dfd = data2fd(log(data/pars[1]),times,bbasis)

coefs = matrix(0,nbasis,3)
coefs[,3] = dfd$coefs
 
res = FitMatchOpt(coefs=coefs,which=1:2,proc=lsproc,pars=pars,meth='nlminb')

res2 = inneropt(data=data,times=times,pars=pars,proc=lsproc,lik=lslik,coefs=res$coefs)

res3 = outeropt(data=data,times=times,pars=pars,proc=lsproc,lik=lslik,coefs=res2$coefs,
  active=c('rho','i','b0','b1','b2'))


## Or just log the observations, this is faster. 

spars = pars
spars[1] = 1

tdata = log(data/pars[1])

dfd = data2fd(tdata,times,bbasis)

coefs = matrix(0,nbasis,3)
coefs[,3] = dfd$coefs
 
res = FitMatchOpt(coefs=coefs,which=1:2,proc=lsproc,pars=spars)

res2 = inneropt(data=tdata,times=times,pars=spars,proc=lsproc,lik=slik,coefs=res$coefs)

res3 = outeropt(data=tdata,times=times,pars=spars,proc=lsproc,lik=slik,coefs=res2$coefs,
  active=c('i','b0','b1','b2'),out.meth='nlminb')

