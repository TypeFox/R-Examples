#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                       Bayesian inference on a proportion
#-----------------------------------------------------------------------------------------

.ws16 = proto(

  create = function(.,h,...) {

    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Bayesian inference\non a categorical distribution") %in% names(.ws$nb)) return()
    
   .$param1 = gedit("1",width=15,coerce.with=as.numeric,handler=.$compute)
    enabled(.$param1) = FALSE

    add(.ws$nb,group <- ggroup(horizontal=FALSE),label=.$translate("Bayesian inference\non a categorical distribution"))
    add(group,tmp <- gframe(.$translate("Dirichlet prior"),horizontal=FALSE))
    priorGroup = glayout(container=tmp)
    priorGroup[2,2,anchor=c(-1,0)]=glabel("Parameter")
    priorGroup[2,3,anchor=c(-1,0)]=.$param1
    priorGroup[4,1]=""
    visible(priorGroup)=TRUE

   .$success = gedit("5 4 6 11 15 21",width=22,handler=.$compute)

    add(group,tmp <- gframe(.$translate("Observed data")))
    dataGroup = glayout(container=tmp)
    dataGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Counts"))
    dataGroup[2,3]=.$success
    dataGroup[4,1]=""
    visible(dataGroup)=TRUE

   .$priorprob = gedit("",width=5)
   .$handler.ID['changePriorProb'] = addhandlerchanged(.$priorprob,.$compute)

   .$model = gedit("",width=16)
   .$handler.ID['changeModel'] = addhandlerchanged(.$model,.$compute)

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
    est.tab = cbind(Levels=rep("",20),Saturated=rep("",20),Target=rep("",20),Averaged=rep("",20))
    add(tmp,.$postestim <- gtable(est.tab),expand=TRUE)
    names(.$postestim) = .$translate(names(.$postestim))

    buttonGroup=ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Compute"),container=buttonGroup, handler=.$compute)
  
  },
  onTestAll = function(.,h,...) {
  
    if(svalue(.$testAll)) enabled(.$model) = FALSE
    else                  enabled(.$model) = TRUE
    svalue(.$model) = ""
  },
  compute = function(.,h,...) {

    # Check parameters
    if(any(is.na(svalue(.$param1)))) {
      gmessage(.$translate("Please specify prior parameters."))
      return()
    }
    
    # Check data input
    if(any(is.na(svalue(.$success)))) {
      gmessage(.$translate("Please specify observed data."))
      return()
    }
    
    a = svalue(.$param1)
   .$s = unlist(strsplit(svalue(.$success)," "))
   .$s = as.numeric(s[s != ""])
   .$N = sum(s)
   .$K = length(.$s)

    if(.$K<2) {
      gmessage(.$translate("Please specify at least two counts."))
      return()
    }

    # Saturated model by default
    if(!nchar(svalue(.$model))) {
      blockhandler(.$model,.$handler.ID['changeModel'])
      svalue(.$model) = paste(1:.$K,collapse=" ")
      unblockhandler(.$model,.$handler.ID['changeModel'])
    }
    m = .$getModel()
    if(length(m) != .$K) {
      gmessage(.$translate("Incorrect number of groups in model definition."))
      return()
    }
    
   .ws$setStatus(.$translate("Estimation in progress..."))
    svalue(.$bestOne) = ""

    if(.$K==2) models = cbind(c(1,1),c(1,2))
    else     models = setparts(.$K)
    nmodels = ncol(models)
    svalue(.$possible) = paste(nmodels)

    # Reference models (null and saturated)
   .$p0 = .$testModel(rep("1",.$K),a)$estimates
   .$ps = .$testModel(1:.$K,a)$estimates

    if(svalue(.$testAll)) {
    
      # Equal prior probs in this case
      blockhandler(.$priorprob,.$handler.ID['changePriorProb'])
      svalue(.$priorprob) = paste("1/",nmodels,sep="")
      unblockhandler(.$priorprob,.$handler.ID['changePriorProb'])
      priorprob = 1/nmodels

      # We start with the constant model
     .$avestimates = rep(1/.$K,.$K)
     .$p1 = rep(1/.$K,.$K)
      sum.bf = 1
      best.bf = 1
      best.model = paste("(",paste(1:.$K,collapse=","),")",sep="")
      best.model.index = 1
      group.names = paste(1:.$K)

      for(i in 2:nmodels) {
      
        # Display model number
        svalue(.$possible) = paste(i,"/",nmodels)
        Sys.sleep(.001)
        
        # Test the model
        test = .$testModel(models[,i],a)
        sum.bf = sum.bf + test$bf
       .$avestimates = .$avestimates + test$bf*test$estimates
        if(test$bf> best.bf) {
          best.bf = test$bf
          model.name = tapply(group.names,models[,i],paste,collapse=",")
          model.name = paste(sapply(model.name,function(x) paste("(",x,")",sep="")),collapse=",")
          best.model = model.name
          best.model.index = i
         .$p1 = test$estimates
        }
      }
     .$avestimates = .$avestimates / sum.bf
      
      # Best model
      svalue(.$bf) = round(best.bf,4)
      svalue(.$bestOne) = best.model
      post.prob = best.bf/sum.bf
      svalue(.$postprob) = round(post.prob,4)
     .$groups = models[,best.model.index]
      
      # Estimates
     .$postestim[,] = matrix("",20,4)
     .$postestim[1:.$K,1] = paste(1:.$K)
     .$postestim[1:.$K,2] = round(.$ps,4)
     .$postestim[1:.$K,3] = round(.$p1,4)
     .$postestim[1:.$K,4] = round(.$avestimates,4)
    } 
    
    # Test a target model against the null
    else {

      # Check prior prob on the target model
      priorprob = eval(parse(text=svalue(.$priorprob)))
      if(is.null(priorprob)) {
        blockhandler(.$priorprob,.$handler.ID['changePriorProb'])
        svalue(.$priorprob) = "1/2"
        unblockhandler(.$priorprob,.$handler.ID['changePriorProb'])
        priorprob = .5
      }
      
      # Convert model symbols in integers if letters were supplied
      if(all(m %in% letters)) m = match(tolower(m),letters)
    
      test = .$testModel(m,a)
      post.prob = priorprob * test$bf/(priorprob * test$bf + 1 - priorprob)
     .$groups = as.numeric(factor(m))
      
      svalue(.$bf) = paste(round(test$bf,4))
      svalue(.$postprob) = paste(round(post.prob,4))

      # A posteriori model averaged estimations
     .$p1 = test$estimates
     .$avestimates = post.prob*.$p1 + (1-post.prob)*.$p0

     .$postestim[,] = matrix("",20,4)
     .$postestim[1:.$K,1] = paste(1:.$K)
     .$postestim[1:.$K,2] = round(.$ps,4)
     .$postestim[1:.$K,3] = round(.$p1,4)
     .$postestim[1:.$K,4] = round(.$avestimates,4)
    }
    
   .$updatePlot(h,...)
   .ws$setStatus(.$translate("Ready."))
   
  },
  updatePlot = function(.,h,...) {
      
    # Saturated model
    plot(1:.$K,.$ps,xlab=.$translate("Response"),ylab=.$translate("Estimated probabilities"),cex=1.5,lwd=2,ylim=c(0,1),main="")

    # Best model
    points((1:.$K)-.02,.$p1,type="h",col="red",lwd=3)
    
    # Averaged model
    points((1:.$K)+.02,.$avestimates,type="h",col="blue",lwd=3)
    
    # Uniform distribution
    abline(h=1/.$K,lty=2,lwd=2,col="lightgray")
    
    legend("topleft",inset=.01,legend=.$translate(c("Saturated","Target","Averaged")),pch=c(1,NA,NA),lty=c(0,1,1),lwd=c(2,3,3),col=c("black","red","blue"))
  },
  getModel = function(.) {
  
    m = unlist(strsplit(svalue(.$model)," "))
    m[m != ""]
  },
  testModel = function(.,mod,a) {
  
    lBeta = function(z)      sum(lgamma(z)) - lgamma(sum(z))
    lPsi  = function(nn,aa)  lBeta(nn+aa) - lBeta(aa)
    lComb = function(nn)     lgamma(sum(nn)+1) - sum(lgamma(nn+1))
    
    # Grouped counts
    nn = tapply(.$s,mod,sum)
    C  = length(nn)
    aa = rep(a,C)
    N  = sum(nn)
    T = table(mod)
    
    # Bayes factor: O'Hagan & Forster (2005) p. 350
    l0 = -N * log(.$K)
    l1 = lPsi(nn,aa) - sum(nn*log(T))
    bf = exp(l1-l0)

    # Estimates (posterior means)
    estimates = (nn/T) + aa
    estimates = estimates[mod]
    estimates = estimates/sum(estimates)

    list(estimates=estimates,bf=bf)
  },

  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },
  #---------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                                    CONTENT
  #---------------------------------------------------------------------------------------
  handler.ID   = list(),               # IDs of various handlers that may have to be blocked
  param1       =  NULL,                #
  success      =  NULL,                # Data as string
  s            =  NULL,                # Data in numeric vector format
  N            =  NULL,                # Total count
  K            =  NULL,                # Number of categories
  priorprob    =  NULL,                #
  model        =  NULL,                #
  p0           =  NULL,                # Probabilities from a a uniform model
  ps           =  NULL,                # Probability posterior estimates from the saturated model
  p1           =  NULL,                # Probability posterior estimates from the best model
  possible     =  NULL,                #
  testAll      =  NULL,                #
  bestOne      =  NULL,                #
  bf           =  NULL,                #
  postprob     =  NULL,                #
  postestim    =  NULL                 #
)


