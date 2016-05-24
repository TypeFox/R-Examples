#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                        Bayesian inference on a standardized effect
#-----------------------------------------------------------------------------------------
.ws14 = proto(

  create = function(.,h,...) {
  
    # Don't activate till it is finished
    # gmessage("Not available yet")
    # return()
    
    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Bayesian inference\non a standardized effect") %in% names(.ws$nb)) return()
    
   .$priorparam1 = gedit("0",width=15,handler=.$updatePlot)
    enabled(.$priorparam1) = FALSE
   .$priorparam2 = gedit("1",width=15,handler=.$updatePlot)
    enabled(.$priorparam2) = FALSE

    add(.ws$nb,group <- ggroup(horizontal=FALSE),label=.$translate("Bayesian inference\non a standardized effect"))
    tmp = gframe(.$translate("Gaussian prior"), horizontal=FALSE,container = group,expand=TRUE)
    priorGroup = glayout(container=tmp)
    priorGroup[2,2,anchor=c(-1,0)] = glabel(.$translate("Mean"))
    priorGroup[2,3] = .$priorparam1
    priorGroup[3,2,anchor=c(-1,0)] = glabel(.$translate("Standard dev."))
    priorGroup[3,3] = .$priorparam2
    priorGroup[4,2] = ""
    visible(priorGroup) = TRUE

   .$xbar  = gedit("",width=15,handler=.$updatePlot)
   .$n     = gedit("",width=15,coerce.with=as.numeric,handler=.$updatePlot)
   .$s     = gedit("",width=15,handler=.$updatePlot)
   .$refmean = gedit("",width=15,handler=.$updatePlot)
   .$delta = glabel("")

    tmp = gframe(.$translate("Observed data"), container = group,expand=TRUE)
    dataGroup = glayout(container=tmp)
    dataGroup[2,2] = "Ref. mean"
    dataGroup[2,3] = .$refmean
    dataGroup[3,2,anchor=c(-1,0)] = glabel(.$translate("Sample size"))
    dataGroup[3,3] = .$n
    dataGroup[4,2,anchor=c(-1,0)] = glabel(.$translate("Mean"))
    dataGroup[4,3] = .$xbar
    dataGroup[5,2,anchor=c(-1,0)]=glabel(.$translate("Standard dev."))
    dataGroup[5,3] = .$s
    dataGroup[6,2,anchor=c(-1,0)]=glabel(.$translate("Standard effect"))
    dataGroup[6,3] = .$delta
    dataGroup[7,2] = ""
    visible(dataGroup)=TRUE

   .$level     = gedit("0.95",width=5,coerce.with=as.numeric,handler=.$updatePlot)
   .$postdelta = glabel("")
   .$interval  = glabel("")

    tmp = gframe(.$translate("Posterior statistics"),container=group,expand=TRUE)
    statsGroup = glayout(container=tmp)
    statsGroup[2,2,anchor=c(-1,0)] = glabel(.$translate("Level"))
    statsGroup[2,3] = .$level
    statsGroup[3,2] = glabel(.$translate("Credibility"))
    statsGroup[3,3:4] = .$interval
    statsGroup[4,2] = glabel(.$translate("Standard effect"))
    statsGroup[4,3] = .$postdelta
    statsGroup[5,2] = ""
    visible(statsGroup)=TRUE

   .$value = gedit("0.2",width=10,handler=.$updatePlot)
   .$refmean = gedit("",width=10,handler=.$updatePlot)
   .$priorprob = gedit("1/2",width=5)
   .$bf = glabel("")
   .$postprob = glabel("")
   .$op = gdroplist(c("H : \u03B4 >","H : \u03B4 <","H : \u03B4 ="),handler=.$updatePlot)

    tmp = gframe(.$translate("Hypothesis test"), container = group,expand=TRUE)
    testGroup = glayout(container=tmp)
    testGroup[2,2] = .$op
    testGroup[2,3] = .$value
    testGroup[3,2,anchor=c(-1,0)] = glabel(.$translate("Prior Pr(H)"))
    testGroup[3,3] = .$priorprob
    testGroup[4,2] = glabel(.$translate("Bayes factor"))
    testGroup[4,3] = .$bf
    testGroup[5,2] = glabel(.$translate("Pr(H|D)"))
    testGroup[5,3] = .$postprob
    testGroup[6,2] = ""
    visible(testGroup) = TRUE

    addSpring(group)

    buttonGroup=ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Compute"),container=buttonGroup, handler=.$updatePlot)

  },
  
  updatePlot = function(.,h,...) {
    
    # Check parameters
    if(any(is.na(c(svalue(.$priorparam1),svalue(.$priorparam2))))) {
      gmessage(.$translate("Please specify all prior parameters."))
      return()
    }
    
    if(any(is.na(c(svalue(.$xbar),svalue(.$s),svalue(.$n))))) {
      gmessage(.$translate("Please specify observed data statistics."))
      return()
    }

    if(is.na(svalue(.$refmean))) {
      gmessage(.$translate("Please provide a reference mean in the data section."))
      return()
    }

    if(is.na(svalue(.$value))) {
      gmessage(.$translate("Please provide a normative value in the test section."))
      return()
    }

    # Get prior parameters
    mu.lambda = eval(parse(text=svalue(.$priorparam1)))
    sigma.lambda = eval(parse(text=svalue(.$priorparam2)))
    
    # Get data
    xbar = eval(parse(text=svalue(.$xbar)))
    n = eval(parse(text=svalue(.$n)))
    s = eval(parse(text=svalue(.$s)))
    
    # Get test info
    conf = svalue(.$level)
    m0 = eval(parse(text=svalue(.$refmean)))
    obs.d = (xbar - m0)/s
    svalue(.$delta) = round(obs.d,3)
    op = svalue(.$op)
    if(is.na(svalue(.$priorprob))) svalue(.$priorprob) ="1/2"
    prior.prob = eval(parse(text=svalue(.$priorprob)))
    
    # Compute posterior statistics
    nu0 = max(prior.n2-1,0)
    post.n1 = n + prior.n1
    nu.n = nu0 + n
    post.mean = (n*xbar + prior.n1*prior.mean)/post.n1
    scale = sqrt(((n-1)*(s**2) + ((xbar-prior.mean)**2)/((1/n)+(1/prior.n1)))/nu.n)
    post.d = (post.mean - m0)/scale
    svalue(.$post.delta) = round(post.d,3)
    
    # Plots: Prior, likelihood, posterior
    lambda.scale = sqrt(1/(prior.n1 + n))
    left.lim  = .$qlambdaprime((1-.995)/2,nu.n,post.d,lambda.scale)
    right.lim = .$qlambdaprime((1+.995)/2,nu.n,post.d,lambda.scale)
    x = seq(left.lim,right.lim,len=500)
    post = .$dlambdaprime(x,nu.n,post.d,lambda.scale)
    q1 = .$qlambdaprime((1-conf)/2,nu.n,post.d,lambda.scale)
    q2 = .$qlambdaprime((1+conf)/2,nu.n,post.d,lambda.scale)
    prior = dnorm(x,prior.mean,prior.sd)
    lh = dt(d.obs*sqrt(n),n-1,ncp=x)*sqrt(n)

    svalue(.$bf) = ""
    svalue(.$postprob) = ""
    
    # Default two-sided test
    BF = (dt(sqrt(prior.n)*(m0-prior.mean)/prior.sd,nu0)/(prior.sd/sqrt(prior.n)))/(dt(sqrt(post.n)*(m0-post.mean)/scale,nu.n)/(scale/sqrt(post.n)))
    post.prob = prior.prob * BF/(prior.prob * BF + 1 - prior.prob)
    
    if(op != "H : \u03B4 =") {
      priorH = pt(sqrt(prior.n)*(m0-prior.mean)/prior.sd,nu0)
      priorA = 1 - priorH
      prior.odds = priorH/priorA
      postH = pt(sqrt(post.n)*(m0-post.mean)/scale,nu.n)
      postA = 1 - postH
      post.odds = postH/postA
      BF = post.odds/prior.odds
      post.prob = postH
      if(op=="H : \u03B4 >") { BF = 1/BF     ; post.prob = postA }
    }        

    svalue(.$bf) = round(BF,4)
    svalue(.$postprob) = round(post.prob,4)

    m = max(c(prior, lh, post))
      
    svalue(.$interval) = paste("[",round(q1,4),";",round(q2,4),"]")
    svalue(.$postmean) = round(post.mean,4)
    svalue(.$postsd)   = round(post.sd,4)

    plot(x, post, type = "l", xlab=.$translate("Standard effect"),ylab = .$translate("Density"), lty = 1, lwd = 2, main = .$translate("Posterior probability"),ylim = c(0,m))
    z = seq(q1,q2,len=500)

    # One-sided regions
    if(!is.null(m0) && (op != "H : \u03B4 =")) {  
      if(op == "H : \u03B4 >") z = seq(m0,max(x),len=500)
      else                     z = seq(min(x),m0,len=500)        
    }
    if(svalue(.$sdFixed)==.$translate("Known")) dz = dnorm(z,post.mean,post.sd)
    else                           dz = dt(sqrt(post.n)*(z-post.mean)/scale,nu.n)/(scale/sqrt(post.n))
    polygon(c(z,max(z),min(z)),c(dz,0,0),density=-1,col="lightgrey",lwd=2)

    lines(x, lh, lty = 2, lwd = 2, col = "red")
    lines(x, prior, lty = 3, lwd = 2, col = "darkgreen")
    if(is.null(m0)) text(post.mean,max(post)/3,cex=1.3,paste(round(conf*100),"%",sep=""))
    legend("topright", .$translate(c("Prior", "Likelihood", "Posterior")), cex=.8,lty = c(3,2,1), lwd = rep(2,3), col = c("darkgreen", "red", "black"),inset=0.01,bg="white")

  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },
  # Lambda-prime distribution (Lecoutre, 1999)
  rlambdaprime = function(.,n,nu,ncp,scale) { rnorm(n,0,scale) + ncp*sqrt(rchisq(n,nu)/nu) },
  plambdaprime = function(.,x,nu,ncp,scale) { 1-pt(ncp/scale,nu,x/scale) },
  qlambdaprime = function(.,p,nu,ncp,scale) {

    t = ncp/scale
    k = exp( ((log(2)-log(nu))/2) + lgamma((nu+1)/2) - lgamma(nu/2) )
    M = k*t
    V = 1 + t**2 - M**2

    N = 100
    from = M - 4.5*sqrt(V)
    to   = M + 4.5*sqrt(V)

    # By interpolation
    for(iter in 1:3) {
    
      range = seq(from,to,len=N)
      index = which.min(abs(p - plambdaprime(range,nu,ncp,scale)))
      from = range[index-1]
      to = range[index+1]
    }
    
    range[index]
  },
  dlambdaprime = function(.,x,nu,ncp,scale) {
  
    # Numerical derivatives
    dx = .0001
    dF1 = plambdaprime(x-dx/2,nu,ncp,scale)
    dF2 = plambdaprime(x+dx/2,nu,ncp,scale)
    
    (dF2-dF1)/dx
  },
  #---------------------------------------------------------------------------------------
  #  SLOT            INITIAL VALUE                           CONTENT
  #---------------------------------------------------------------------------------------
  priorparam1     =  NULL,                     #
  priorparam2     =  NULL,                     #
  refmean         =  NULL,                     #
  xbar            =  NULL,                     #
  n               =  NULL,                     #
  s               =  NULL,                     #
  delta           =  NULL,                     #
  level           =  NULL,                     #
  postdelta       =  NULL,                     #
  postsd          =  NULL,                     #
  interval        =  NULL,                     #
  value           =  NULL,                     #
  priorprob       =  NULL,                     #
  bf              =  NULL,                     #
  postprob        =  NULL,                     #
  op              =  NULL                      #
)


