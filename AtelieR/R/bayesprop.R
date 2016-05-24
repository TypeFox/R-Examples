#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                           Bayesian inference on a proportion
#-----------------------------------------------------------------------------------------
.ws6 = proto(

  create = function(.,h,...) {

    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Bayesian inference\non a proportion") %in% names(.ws$nb)) return()
    
   .$betaparam1 = gedit("1",width=15,coerce.with=as.numeric,handler=.$updatePlot)
   .$betaparam2 = gedit("1",width=15,coerce.with=as.numeric,handler=.$updatePlot)

    add(.ws$nb,group <- ggroup(horizontal=FALSE),label=.$translate("Bayesian inference\non a proportion"))
    tmp = gframe(.$translate("Beta prior"), container = group,horizontal=FALSE)
    priorGroup = glayout(container=tmp)
    priorGroup[2,2,anchor=c(-1,0)]=glabel("a=")
    priorGroup[2,3,anchor=c(-1,0)]=.$betaparam1
    priorGroup[3,2,anchor=c(-1,0),anchor=c(-1,0)]=glabel("b=")
    priorGroup[3,3,anchor=c(-1,0)]=.$betaparam2
    visible(priorGroup)=TRUE

   .$success = gedit("6",width=15,handler=.$updatePlot)
   .$Ntot = gedit("10",width=15,handler=.$updatePlot)

    tmp = gframe(.$translate("Observed data"), container = group)
    dataGroup = glayout(container=tmp)
    dataGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Successes"))
    dataGroup[2,3]=.$success
    dataGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Total"))
    dataGroup[3,3]=.$Ntot
    visible(dataGroup)=TRUE

   .$level = gedit("0.95",width=5,coerce.with=as.numeric,handler=.$updatePlot)
   .$postmean = glabel("")
   .$postmode = glabel("")
   .$postsd = glabel("")
   .$postconf = gedit("0.05")
   .$interval = glabel("")

    tmp = gframe(.$translate("Posterior statistics"), container = group)
    statsGroup = glayout(container=tmp)
    statsGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Level"))
    statsGroup[2,3]=.$level
    statsGroup[3,2]=glabel(.$translate("Credibility"))
    statsGroup[3,3:4]=.$interval
    statsGroup[4,2]=glabel(.$translate("Mean"))
    statsGroup[4,3]=.$postmean
    statsGroup[5,2]=glabel(.$translate("Mode"))
    statsGroup[5,3]=.$postmode
    statsGroup[6,2]=glabel(.$translate("Standard deviation"))
    statsGroup[6,3]=.$postsd
    visible(statsGroup)=TRUE

   .$value = gedit("1/2",width=10,handler=.$updatePlot)
   .$priorprob = gedit("0.5",width=5,handler=.$updatePlot)
   .$bf = glabel("")
   .$postprob = glabel("")
   .$op = gdroplist(c("H : pi <","H : pi =","H : pi >"),handler=.$updatePlot)
    svalue(.$op,index=TRUE) = 2

    tmp = gframe(.$translate(.$translate("Hypothesis test")), container = group)
    testGroup = glayout(container=tmp)
    testGroup[2,2,anchor=c(-1,0)]=.$op
    testGroup[2,3]=.$value
    testGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Prior Pr(H)"))
    testGroup[3,3]=.$priorprob
    testGroup[4,2]=glabel(.$translate("Bayes Factor"))
    testGroup[4,3]=.$bf
    testGroup[5,2]=glabel(.$translate("Pr(H|D)"))
    testGroup[5,3]=.$postprob
    visible(testGroup)=TRUE

    addSpring(group)

    buttonGroup=ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Compute"),container=buttonGroup, handler=.$updatePlot)
  
  },
  
  updatePlot = function(.,h,...) {

    # Vérification des paramètres
    if(any(is.na(c(svalue(.$betaparam1),svalue(.$betaparam2))))) {
      gmessage(.$translate("Please specify prior parameters."))
      return()
    }
    
    # Vérification des paramètres
    if(any(is.na(c(svalue(.$success),svalue(.$Ntot))))) {
      gmessage(.$translate("Please specify observed data."))
      return()
    }
      
    a = svalue(.$betaparam1)
    b = svalue(.$betaparam2)
    s = eval(parse(text=svalue(.$success)))
    N = eval(parse(text=svalue(.$Ntot)))
    f = N-s
    p0 = eval(parse(text=svalue(.$value)))
    conf = svalue(.$level)
    priorprobH0 = eval(parse(text=svalue(.$priorprob)))
    op = svalue(.$op)

    p = seq(0.005, 0.995, length = 500)
    prior = dbeta(p, a, b)
    post = dbeta(p, a + s, b + f)
    lh = dbeta(p, s + 1, f + 1)
    m = max(c(prior, lh, post))
    q1 = qbeta((1-conf)/2,a + s, b + f)
    q2 = qbeta((1+conf)/2,a + s, b + f)
    
    postmode = (a+s-1)/(a+s+b+f-2)
    postmean = (a+s)/(a+s+b+f)
    postvar = postmean*(1-postmean)/(a+s+b+f+1)
    postsd = sqrt(postvar)

    svalue(.$interval) = paste("[",round(q1,4),";",round(q2,4),"]")
    svalue(.$postmode) = round(postmode,4)
    svalue(.$postmean) = round(postmean,4)
    svalue(.$postsd) = round(postsd,4)
    
    # Bayesian test
    svalue(.$bf) = ""
    svalue(.$postprob) = ""
    if(!is.null(p0) && !is.null(priorprobH0)) {
      if(op == "H : pi =") {
        # Savage-Dickey ratio
        bf = dbeta(p0,a+s,b+f)/dbeta(p0,a,b)
        post.prob = priorprobH0 * bf/(priorprobH0 * bf + 1 - priorprobH0)
      }
      else  {
        priorH = pbeta(p0,a,b)
        priorA = 1 - priorH
        prior.odds = priorH/priorA
        postH = pbeta(p0,a+s,b+f)
        postA = 1 - postH
        post.odds = postH/postA
        bf = post.odds/prior.odds
        post.prob = postH
        if(op=="H : pi >") { bf = 1/bf     ; post.prob = postA }
      }
      svalue(.$bf) = round(bf,4)
      svalue(.$postprob) = round(post.prob,4)
    }

    plot(p, post, type = "l", xlab=.$translate("Probability of success"),ylab = .$translate("Density"), lty = 1, lwd = 2, main = .$translate("Posterior probability"), xlim=c(0,1),ylim = c(0, max(m,5)))
    z = seq(q1,q2,len=500)

    # One-sided regions
    if(!is.null(p0) && (op != "H : pi =")) {  
      if(op == "H : pi >") z = seq(p0,max(p),len=500)
      else                 z = seq(min(p),p0,len=500)        
    }
    
    dz = dbeta(z,a+s,b+f)
    polygon(c(z,max(z),min(z)),c(dz,0,0),density=-1,col="lightgrey",lwd=2)

    lines(p, lh, lty = 2, lwd = 2, col = "red")
    lines(p, prior, lty = 3, lwd = 2, col = "darkgreen")
    if(op == "H : pi =") text((q2+q1)/2,max(post)/3,cex=1.3,paste(round(conf*100),"%",sep=""))
    legend("topright", .$translate(c("Prior", "Likelihood", "Posterior")), cex=.8,lty = c(3,2,1), lwd = rep(2,3), col = c("darkgreen", "red", "black"),inset=0.01,bg="white")
  
  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },
  #---------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                                    CONTENT
  #---------------------------------------------------------------------------------------
  success      =  NULL,                #
  Ntot         =  NULL,                #
  level        =  NULL,                #
  postmean     =  NULL,                #
  postmode     =  NULL,                #
  postsd       =  NULL,                #
  postconf     =  NULL,                #
  interval     =  NULL,                #
  value        =  NULL,                #
  priorprob    =  NULL,                #
  bf           =  NULL,                #
  postprob     =  NULL                 #
)


