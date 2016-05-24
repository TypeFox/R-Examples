#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                       Bayesian inference on a proportion
#-----------------------------------------------------------------------------------------

.ws9 = proto(

  create = function(.,h,...) {

    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Bayesian inference\non several proportions") %in% names(.ws$nb)) return()
    
   .$betaparam1 = gedit("1",width=15,coerce.with=as.numeric,handler=.$updatePlot)
   .$betaparam2 = gedit("1",width=15,coerce.with=as.numeric,handler=.$updatePlot)
    enabled(.$betaparam1) = FALSE
    enabled(.$betaparam2) = FALSE

    add(.ws$nb,group <- ggroup(horizontal=FALSE),label=.$translate("Bayesian inference\non several proportions"))
    add(group,tmp <- gframe(.$translate("Beta prior"),horizontal=FALSE))
    priorGroup = glayout(container=tmp)
    priorGroup[2,2,anchor=c(-1,0)]=glabel("a=")
    priorGroup[2,3,anchor=c(-1,0)]=.$betaparam1
    priorGroup[3,2,anchor=c(-1,0),anchor=c(-1,0)]=glabel("b=")
    priorGroup[3,3,anchor=c(-1,0)]=.$betaparam2
    priorGroup[4,1]=""
    visible(priorGroup)=TRUE

   .$success = gedit("",width=20,handler=.$updatePlot)
   .$Ntot = gedit("",width=20,handler=.$updatePlot)

    add(group,tmp <- gframe(.$translate("Observed data")))
    dataGroup = glayout(container=tmp)
    dataGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Successes"))
    dataGroup[2,3]=.$success
    dataGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Total"))
    dataGroup[3,3]=.$Ntot
    dataGroup[4,1]=""
    visible(dataGroup)=TRUE

   .$priorprob = gedit("")
   .$handler.ID['changePriorProb'] = addhandlerchanged(.$priorprob,.$updatePlot)
   .$model = gedit("",width=15,handler=.$updatePlot)
   .$bf = glabel("")
   .$possible = glabel("")
   .$postprob = glabel("")
   .$testAll = gcheckbox(.$translate("Test all"),checked=FALSE,handler=.$onTestAll)
   .$bestOne = glabel("")

    add(group,tmp <- gframe(.$translate("Hypothesis test")))
    testGroup = glayout(container=tmp)
    testGroup[2,2,anchor=c(-1,0)] = glabel(.$translate("Possible models"))
    testGroup[2,3] = .$possible
    testGroup[3,2,anchor=c(-1,0)] = glabel(.$translate("Target model"))
    testGroup[3,3] = .$model
    testGroup[4,2,anchor=c(-1,0)] = glabel(.$translate("Prior Pr(M)"))
    testGroup[4,3] =.$priorprob
    testGroup[5,3] =.$testAll
    testGroup[6,2] = glabel(.$translate("Best model"))
    testGroup[6,3] =.$bestOne
    testGroup[7,2] = glabel(.$translate("Bayes factor"))
    testGroup[7,3] =.$bf
    testGroup[8,2] = glabel(.$translate("Pr(M|D)"))
    testGroup[8,3] =.$postprob
    testGroup[9,1] = ""
    visible(testGroup)=TRUE

    add(group, tmp <- gframe(.$translate("Model posterior estimates")),expand=TRUE)
    est.tab = cbind(Group=rep("",10),Saturated=rep("",10),Target=rep("",10),Averaged=rep("",10))
    add(tmp,.$postestim <- gtable(est.tab),expand=TRUE)
    names(.$postestim) = .$translate(names(.$postestim))

    buttonGroup=ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Compute"),container=buttonGroup, handler=.$updatePlot)
  
  },
  onTestAll = function(.,h,...) {
  
    if(svalue(.$testAll)) enabled(.$model) = FALSE
    else                  enabled(.$model) = TRUE
    svalue(.$model) = ""
    
  },
  updatePlot = function(.,h,...) {

    # Check parameters
    if(any(is.na(c(svalue(.$betaparam1),svalue(.$betaparam2))))) {
      gmessage(.$translate("Please specify prior parameters."))
      return()
    }
    
    # Check data input
    if(any(is.na(c(svalue(.$success),svalue(.$Ntot))))) {
      gmessage(.$translate("Please specify observed data."))
      return()
    }
    
    a = svalue(.$betaparam1)
    b = svalue(.$betaparam2)
    s = unlist(strsplit(svalue(.$success)," "))
    s = as.numeric(s[s != ""])
    N = unlist(strsplit(svalue(.$Ntot)," "))
    N = as.numeric(N[N != ""])
    f = N - s

    if( (length(s)<2) || (length(f)<2) || (length(s) != length(f)) ) {
      gmessage(.$translate("Please specify at least two counts in each field.\nMake sure there are as many values in both fields."))
      return()
    }

    if( (length(a)>1) || (length(b)>1) ) {
      gmessage(.$translate("Please specify a single prior parameter value in each field."))
      return()
    }
    
    K = length(s)

    # Check model definition
    if(!nchar(svalue(.$model))) {
      svalue(.$model) = paste(1:length(s),collapse=" ")
      # Return or the analysis is performed twice!
      return()
    }
    m = .$getModel()
    if(length(m) != K) {
      gmessage(.$translate("Incorrect number of groups in model definition."))
      return()
    }
    
   .ws$setStatus(.$translate("Estimation in progress..."))
    svalue(.$bestOne) = ""

    # Only two groups -> two models
    if(K==2) models = cbind(c(1,1),c(1,2))
    
    # More than two groups
    else     models = setparts(K)
    nmodels = ncol(models)
    svalue(.$possible) = paste(nmodels)

    p0 = (sum(s)+a)/(sum(s+f)+a+b)

    # Test all possible models
    if(svalue(.$testAll)) {
    
      # Equal prior probs in this case
      blockhandler(.$priorprob,.$handler.ID['changePriorProb'])
      svalue(.$priorprob) = paste("1/",nmodels,sep="")
      unblockhandler(.$priorprob,.$handler.ID['changePriorProb'])

      priorprob = 1/nmodels

      results = matrix(0,0,K+1)
      group.names = paste(1:K)
      for(i in 1:ncol(models)) {
        test = .$testModel(s,f,models[,i])
        est = (test$s+a)/(test$s+test$f+a+b)
        results = rbind(results,c(est,test$bf))
        model.name = tapply(group.names,models[,i],paste,collapse=",")
        model.name = paste(sapply(model.name,function(x) paste("(",x,")",sep="")),collapse=",")
        rownames(results)[i] = model.name
      }
      results = cbind(results,results[,K+1]/sum(results[,K+1]))
      colnames(results) = c(paste("G",1:K,sep=""),"BF","Prob")
      
      # Best model
      best.one = which.max(results[,"Prob"])
      test = .$testModel(s,f,models[,best.one])
      svalue(.$bestOne) = rownames(results)[best.one]
      best.bf = results[best.one,"BF"]
      post.prob = results[best.one,"Prob"]
      best.estim = results[best.one,1:K]
      svalue(.$bf) = round(best.bf,4)
      svalue(.$postprob) = round(post.prob,4)

      # A posteriori model averaged estimations (Viana, 1991)
      psat = (s+a)/(s+f+a+b)
      post.estim = colSums(results[,"Prob"]*results[,1:K])

     .$postestim[,] = matrix("",10,4)
     .$postestim[1:K,1] = paste(1:K)
     .$postestim[1:K,2] = round(psat,4)
     .$postestim[1:K,3] = round(best.estim,4)
     .$postestim[1:K,4] = round(post.estim,4)
    }
    else {

      priorprob = eval(parse(text=svalue(.$priorprob)))

      if(is.null(priorprob)) {
        blockhandler(.$priorprob,.$handler.ID['changePriorProb'])
        svalue(.$priorprob) = "1/2"
        unblockhandler(.$priorprob,.$handler.ID['changePriorProb'])
        priorprob = .5
      }
      
      # Convert model symbols in integers if letters were supplied
      if(all(m %in% letters)) m = match(tolower(m),letters)
    
      test = .$testModel(s,f,m)
      post.prob = priorprob * test$bf/(priorprob * test$bf + 1 - priorprob)
      
      svalue(.$bf) = paste(round(test$bf,4))
      svalue(.$postprob) = paste(round(post.prob,4))

      # A posteriori model averaged estimations (Viana, 1991)
      psat = (s+a)/(s+f+a+b)
      p1 = (test$s+a)/(test$s+test$f+a+b)
      post.estim = post.prob*p1 + (1-post.prob)*p0

     .$postestim[,] = matrix("",10,4)
     .$postestim[1:K,1] = paste(1:K)
     .$postestim[1:K,2] = round(psat,4)
     .$postestim[1:K,3] = round(p1,4)
     .$postestim[1:K,4] = round(post.estim,4)
    }
    
    # Plot
    stats = .$postestim[1:K,]
    plot(stats[,1],stats[,2],xlab=.$translate("Groups"),ylab=.$translate("Estimated probabilities"),ylim=c(0,1),main="",type="n",xaxt="n")
    axis(1,1:K,paste(1:K))
    abline(h=p0,col="lightgrey",lty=2,lwd=2)
    points(stats[,1],stats[,2],cex=1.3)             # Saturated
    for(k in 1:K) lines(rbind(c(k,qbeta(.025,test$s[k]+a,test$f[k]+b)),c(k,qbeta(.975,test$s[k]+a,test$f[k]+b))))
    points(stats[,1],stats[,3],pch=19,col="red")    # Target
    points(stats[,1],stats[,4],pch=19,col="blue")   # Averaged
    legend("topleft",pch=c(1,19,19),col=c("black","red","blue"),legend=.$translate(c("Saturated","Target","Averaged")),inset=.01)
    
   .ws$setStatus(.$translate("Ready."))
  },
  getModel = function(.) {
  
    m = unlist(strsplit(svalue(.$model)," "))
    m[m != ""]
  },
  testModel = function(.,s0,f0,mod) {
  
    # Test
    s = tapply(s0,mod,sum)
    f = tapply(f0,mod,sum)
    bf = .$BF10(s,s+f)
    
    # Counts
    s = s[mod]
    f = f[mod]

    list(bf=bf,s=s,f=f)
  },
  # Bayes factor
  BF10 = function(.,k,n) {

    N = sum(n)
    K = sum(k)
      
    log.bf10 = log(N+1) + lchoose(N,K) - sum(lchoose(n,k)) - sum(log(n+1))
    exp(log.bf10)
  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },
  #---------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                                    CONTENT
  #---------------------------------------------------------------------------------------
  handler.ID   = list(),               # IDs of various handlers that may have to be blocked
  betaparam1   =  NULL,                # First parameter of the beta prior (default 1)
  betaparam2   =  NULL,                # Second parameter of the beta prior (default 1)
  success      =  NULL,                #
  Ntot         =  NULL,                #
  priorprob    =  NULL,                #
  model        =  NULL,                #
  possible     =  NULL,                #
  testAll      =  NULL,                #
  bestOne      =  NULL,                #
  bf           =  NULL,                #
  postprob     =  NULL,                #
  postestim    =  NULL                 #
)


