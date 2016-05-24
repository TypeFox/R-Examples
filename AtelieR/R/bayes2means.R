#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#               Bayesian inference on the mean (informative prior)
#-----------------------------------------------------------------------------------------
.ws15 = proto(

  create = function(.,h,...) {
  
    # Don't activate till it is finished
    gmessage("Not available yet")
    return()

    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Bayesian inference\non a sample mean") %in% names(.ws$nb)) return()
    
   .$priorparam1 = gedit("",width=15,handler=.$updatePlot)
   .$priorparam2 = gedit("",width=15,handler=.$updatePlot)
   .$priorn      = gedit("2",width=15,coerce.with=as.numeric,handler=.$updatePlot)

    add(.ws$nb,group <- ggroup(horizontal=FALSE),label=.$translate("Bayesian inference\non a sample mean"))
    tmp = gframe(.$translate("Gaussian prior"), horizontal=FALSE,container = group,expand=TRUE)
    priorGroup = glayout(container=tmp)
    priorGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Prior size"))
    priorGroup[2,3]=.$priorn
    priorGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Mean"))
    priorGroup[3,3]=.$priorparam1
    priorGroup[4,2,anchor=c(-1,0)]=glabel(.$translate("Standard dev."))
    priorGroup[4,3]=.$priorparam2
    visible(priorGroup)=TRUE

   .$xbar = gedit("",width=15,handler=.$updatePlot)
   .$n = gedit("",width=15,coerce.with=as.numeric,handler=.$updatePlot)
   .$s = gedit("",width=15,handler=.$updatePlot)
    enabled(.$s) = FALSE
   .$sdFixed = gradio(.$translate(c("Known","Estimated")),horizontal=TRUE,handler=.$onCheck)

    tmp = gframe(.$translate("Observed data"), container = group,expand=TRUE)
    dataGroup = glayout(container=tmp)
    dataGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Sample size"))
    dataGroup[2,3]=.$n
    dataGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Mean"))
    dataGroup[3,3]=.$xbar
    dataGroup[4,2,anchor=c(-1,0)]=glabel(.$translate("Standard dev."))
    dataGroup[4,3]=.$s
    dataGroup[5,3]=.$sdFixed
    visible(dataGroup)=TRUE

   .$level = gedit("0.95",width=5,coerce.with=as.numeric,handler=.$updatePlot)
   .$postmean = glabel("")
   .$postsd = glabel("")
   .$interval = glabel("")

    tmp = gframe(.$translate("Posterior statistics"), container = group,expand=TRUE)
    statsGroup = glayout(container=tmp)
    statsGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Level"))
    statsGroup[2,3]=.$level
    statsGroup[3,2]=glabel(.$translate("Credibility"))
    statsGroup[3,3:4]=.$interval
    statsGroup[4,2]=glabel(.$translate("Mean"))
    statsGroup[4,3]=.$postmean
    statsGroup[5,2]=glabel(.$translate("Standard dev."))
    statsGroup[5,3]=.$postsd
    visible(statsGroup)=TRUE

   .$value = gedit("",width=10,handler=.$updatePlot)
   .$priorprob = gedit("0.5",width=5,handler=.$updatePlot)
   .$bf = glabel("")
   .$postprob = glabel("")
   .$op = gdroplist(c("H : mu <","H : mu =","H : mu >"),handler=.$updatePlot)

    tmp = gframe(.$translate("Hypothesis test"), container = group,expand=TRUE)
    testGroup = glayout(container=tmp)
    testGroup[2,2]=.$op
    testGroup[2,3]=.$value
    testGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Prior Pr(H)"))
    testGroup[3,3]=.$priorprob
    testGroup[4,2]=glabel(.$translate("Bayes factor"))
    testGroup[4,3]=.$bf
    testGroup[5,2]=glabel(.$translate("Pr(H|D)"))
    testGroup[5,3]=.$postprob
    visible(testGroup)=TRUE

    addSpring(group)

    buttonGroup=ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Compute"),container=buttonGroup, handler=.$updatePlot)

  },
  
  onCheck = function(.,h,...) {
    
    if(svalue(.$sdFixed)==.$translate("Known")) {
      svalue(.$s) = ""
      enabled(.$s) = FALSE
    }
    else {
      enabled(.$s) = TRUE
    }
  },
  
  updatePlot = function(.,h,...) {
    
    # Check parameters
    if(any(is.na(c(svalue(.$priorn),svalue(.$priorparam1),svalue(.$priorparam2))))) {
      gmessage(.$translate("Please specify prior parameters."))
      return()
    }
    
    if(any(is.na(c(svalue(.$xbar),svalue(.$n))))) {
      gmessage(.$translate("Please specify observed data."))
      return()
    }

    # Get input info
    xbar = eval(parse(text=svalue(.$xbar)))
    n = svalue(.$n)
    s = eval(parse(text=svalue(.$s)))
    conf = svalue(.$level)
    prior.mean = eval(parse(text=svalue(.$priorparam1)))
    prior.sd = eval(parse(text=svalue(.$priorparam2)))
    prior.var = prior.sd^2
    prior.precision = 1/prior.var
    prior.n = svalue(.$priorn)
    m0 = eval(parse(text=svalue(.$value)))
    op = svalue(.$op)
    prior.prob = eval(parse(text=svalue(.$priorprob)))
    post.n = n+prior.n
        
    # Case 1: Variance known
    if(svalue(.$sdFixed)==.$translate("Known")) {
      
      post.var = prior.var / (n+prior.n)
      post.sd = sqrt(post.var)
      post.mean = (xbar * n + prior.mean*prior.n)/(n+prior.n)
      x = seq(post.mean-3*post.sd,post.mean+3*post.sd,len=500)
      prior = dnorm(x,prior.mean,prior.sd/sqrt(prior.n))
    
      post = dnorm(x,post.mean,post.sd)
      q1 = qnorm((1-conf)/2,post.mean,post.sd)
      q2 = qnorm((1+conf)/2, post.mean, post.sd)   
      lh = dnorm(xbar,x,prior.sd/sqrt(n))

      svalue(.$bf) = ""
      svalue(.$postprob) = ""
      
      # Bayesian test (if some normative value has been provided)
      if(!is.null(m0) && !is.null(prior.prob)) {
      
        # Default two-sided test
        BF = dnorm(m0,post.mean,post.sd)/dnorm(m0,prior.mean,prior.sd/sqrt(prior.n))
        post.prob = prior.prob * BF/(prior.prob * BF + 1 - prior.prob)
        
        # One-sided test
        if(op != "H : mu =") {
          priorH = pnorm(m0, prior.mean, prior.sd/sqrt(prior.n))
          priorA = 1 - priorH
          prior.odds = priorH/priorA
          postH = pnorm(m0, post.mean, post.sd)
          postA = 1 - postH
          post.odds = postH/postA
          BF = post.odds/prior.odds
          post.prob = postH
          if(op=="H : mu >") { BF = 1/BF     ; post.prob = postA }
        }
        svalue(.$bf) = round(BF,4)
        svalue(.$postprob) = round(post.prob,4)
      }
    }
    
    # Case 2: Variance unknown
    else {
    
      # Check parameters
      if(is.null(s)) {
        gmessage(.$translate("Please provide the observed standard deviation."))
        return()
      }
      nu0 = prior.n - 1

      if(nu0<1) {
        gmessage(.$translate("The prior size is too low (0 df.)."))
        return()
      }
    
      nu.n = nu0 + n
      post.mean = (n*xbar + prior.n*prior.mean)/post.n
      scale = sqrt((nu0*prior.var + (n-1)*s**2 + (prior.n*n)*((xbar-prior.mean)**2)/post.n)/nu.n)
      post.sd = (scale/sqrt(post.n)) * sqrt((post.n-1)/(post.n-3))
      
      left.lim  = qt((1-.995)/2,nu.n)*(scale/sqrt(post.n)) + post.mean
      right.lim = qt((1+.995)/2,nu.n)*(scale/sqrt(post.n)) + post.mean      
      x = seq(left.lim,right.lim,len=500)
      post = dt(sqrt(post.n)*(x-post.mean)/scale,nu.n)/(scale/sqrt(post.n))
      q1 = qt((1-conf)/2,nu.n)*(scale/sqrt(post.n)) + post.mean
      q2 = qt((1+conf)/2,nu.n)*(scale/sqrt(post.n)) + post.mean
      prior = dt(sqrt(prior.n)*(x-prior.mean)/prior.sd,nu0)/(prior.sd/sqrt(prior.n))
      lh = dt(sqrt(n)*(xbar-x)/s,n-1)/(scale/sqrt(post.n))

      svalue(.$bf) = ""
      svalue(.$postprob) = ""
      
      # Bayesian test (if some normative value has been provided)
      if(!is.null(m0) && !is.null(prior.prob)) {
      
        # Default two-sided test
        BF = (dt(sqrt(prior.n)*(m0-prior.mean)/prior.sd,nu0)/(prior.sd/sqrt(prior.n)))/(dt(sqrt(post.n)*(m0-post.mean)/scale,nu.n)/(scale/sqrt(post.n)))
        post.prob = prior.prob * BF/(prior.prob * BF + 1 - prior.prob)
        
        if(op != "H : mu =") {
          priorH = pt(sqrt(prior.n)*(m0-prior.mean)/prior.sd,nu0)
          priorA = 1 - priorH
          prior.odds = priorH/priorA
          postH = pt(sqrt(post.n)*(m0-post.mean)/scale,nu.n)
          postA = 1 - postH
          post.odds = postH/postA
          BF = post.odds/prior.odds
          post.prob = postH
          if(op=="H : mu >") { BF = 1/BF     ; post.prob = postA }
        }        

      svalue(.$bf) = round(BF,4)
      svalue(.$postprob) = round(post.prob,4)
      }
    }

    m = max(c(prior, lh, post))
      
    svalue(.$interval) = paste("[",round(q1,4),";",round(q2,4),"]")
    svalue(.$postmean) = round(post.mean,4)
    svalue(.$postsd)   = round(post.sd,4)

    plot(x, post, type = "l", xlab=.$translate("Mean"),ylab = .$translate("Density"), lty = 1, lwd = 2, main = .$translate("Posterior probability"),ylim = c(0,m))
    z = seq(q1,q2,len=500)

    # One-sided regions
    if(!is.null(m0) && (op != "H : mu =")) {  
      if(op == "H : mu >") z = seq(m0,max(x),len=500)
      else                z = seq(min(x),m0,len=500)        
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

  #---------------------------------------------------------------------------------------
  #  SLOT            INITIAL VALUE                           CONTENT
  #---------------------------------------------------------------------------------------
  priorparam1     =  NULL,                     #
  priorparam2     =  NULL,                     #
  priorn          =  NULL,                     #
  xbar            =  NULL,                     #
  n               =  NULL,                     #
  s               =  NULL,                     #
  sdFixed         =  NULL,                     #
  level           =  NULL,                     #
  postmean        =  NULL,                     #
  postsd          =  NULL,                     #
  interval        =  NULL,                     #
  value           =  NULL,                     #
  priorprob       =  NULL,                     #
  bf              =  NULL,                     #
  postprob        =  NULL,                     #
  op              =  NULL                      #
)


