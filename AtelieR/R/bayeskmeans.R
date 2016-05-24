#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                        Bayesian inference on several means
#-----------------------------------------------------------------------------------------

.ws13 = proto(

  create = function(.,h,...) {

    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Bayesian inference\non several means") %in% names(.ws$nb)) return()
    
    add(.ws$nb,group <- ggroup(horizontal=FALSE),label=.$translate("Bayesian inference\non several means"))

   .$means = gedit("",width=25,handler=.$updatePlot)
   .$sds   = gedit("",width=25,handler=.$updatePlot)
   .$Ntot  = gedit("",width=25,handler=.$updatePlot)

    add(group,tmp <- gframe(.$translate("Observed data")))
    dataGroup = glayout(container=tmp)
    dataGroup[2,2,anchor=c(-1,0)] = glabel(.$translate("Means"))
    dataGroup[2,3] = .$means
    dataGroup[3,2,anchor=c(-1,0)] = glabel(.$translate("Standard dev."))
    dataGroup[3,3] = .$sds
    dataGroup[4,2,anchor=c(-1,0)] = glabel(.$translate("Sample sizes"))
    dataGroup[4,3] = .$Ntot
    dataGroup[5,1]= ""
    visible(dataGroup) = TRUE

   .$priorprob = gedit("",width=5)
   .$handler.ID['changePriorProb'] = addhandlerchanged(.$priorprob,.$updatePlot)
   .$model = gedit("",width=15)
   .$handler.ID['changeModel'] = addhandlerchanged(.$model,.$updatePlot)
   
   .$bf = glabel("")
   .$bartlett = glabel("")
   .$possible = glabel("")
   .$postprob = glabel("")
   .$testAll = gcheckbox(.$translate("Test all"),checked=FALSE,handler=.$onTestAll)
   .$bestOne = glabel("")
   .$doBMA = gcheckbox(.$translate("Average"),handler=.$updatePlot)

    add(group,tmp <- gframe(.$translate("Hypothesis test"),horizontal=FALSE))
    testGroup = glayout(container=tmp)
    testGroup[2,2,anchor=c(-1,0)] = glabel(.$translate("Homogeneity of variances"))
    testGroup[2,3] = .$bartlett
    testGroup[3,2,anchor=c(-1,0)] = glabel(.$translate("Possible models"))
    testGroup[3,3] = .$possible
    testGroup[4,2,anchor=c(-1,0)] = glabel(.$translate("Target model"))
    testGroup[4,3] = .$model
    testGroup[5,2,anchor=c(-1,0)] = glabel(.$translate("Prior Pr(M)"))
    testGroup[5,3] =.$priorprob
    testGroup[6,2] =.$testAll
    testGroup[6,3] =.$doBMA
    # testGroup[10,1] = ""
    visible(testGroup)=TRUE

    buttonGroup1 = ggroup(container=tmp)
    addSpring(buttonGroup1)
    gbutton(.$translate("Compute"),container=buttonGroup1, handler=.$updatePlot)

    add(group, tmp <- gframe(.$translate("Model and posterior estimates"),horizontal=FALSE),expand=TRUE)
    resGroup = glayout(container=tmp)
    resGroup[2,2] = glabel(.$translate("Best model"))
    resGroup[2,3] =.$bestOne
    resGroup[3,2] = glabel(.$translate("BIC"))
    resGroup[3,3] =.$bf
    resGroup[4,2] = glabel(.$translate("Pr(M|D)"))
    resGroup[4,3] =.$postprob

    est.tab = cbind(Group=rep("",10),Means=rep("",10),Inf95=rep("",10),Sup95=rep("",10))
    colnames(est.tab) = .$translate(colnames(est.tab))
    add(tmp,.$postestim <- gtable(est.tab),expand=TRUE)

    buttonGroup2 = ggroup(container=tmp)
    addSpring(buttonGroup2)
    gbutton(.$translate("Details"),container=buttonGroup2, handler=.$printModels)
  
  },
  onTestAll = function(.,h,...) {
  
    if(svalue(.$testAll)) { 
  
      enabled(.$model)     = FALSE
      enabled(.$priorprob) = FALSE
    }
    else { 
      enabled(.$model)     = TRUE
      enabled(.$priorprob) = TRUE
      enabled(.$doBMA)     = FALSE
    }
    svalue(.$priorprob)  = ""
  },
  updatePlot = function(.,h,...) {

    require(partitions)
    
    if(svalue(.$doBMA) && !svalue(.$testAll)) {
      svalue(.$testAll) = TRUE
      # Return or the analysis is performed twice!
      return()
    }
    
    # Any data provided?
    if(any(is.na(c(svalue(.$means),svalue(.$sds))))) {
      gmessage(.$translate("Please specify observed data."))
      return()
    }
    
    mn = unlist(strsplit(svalue(.$means)," "))
    mn = as.numeric(mn[mn != ""])
    sd = unlist(strsplit(svalue(.$sds)," "))
    sd = as.numeric(sd[sd != ""])
    N = unlist(strsplit(svalue(.$Ntot)," "))
    N = as.numeric(N[N != ""])

    # Check data input
    if( (length(mn)<2) || (length(sd)<2) || (length(N)<2) || (var(c(length(mn),length(sd),length(N)))!=0) ) {
      gmessage(.$translate("Please specify at least two groups.\nMake sure there are as many values in all three fields."))
      return()
    }
        
    # Number of groups and observations
    K = length(mn)
    No = sum(N)

    # Bartlett test for homogeneity of variances
    S2 = sd^2
    Sp2 = sum((N-1)*S2)/(No-K)
    A = (No-K)*log(Sp2) - sum((N-1)*log(S2))
    B = 1+(sum(1/(N-1))-(1/(No-K)))/(3*(K-1))
    K2 = A/B
    svalue(.$bartlett) = paste("K2=",round(K2,3),"p<",round(1-pchisq(K2,K-1),4))

    # Default model is the saturated one
    if(!nchar(svalue(.$model))) {
      blockhandler(.$model,.$handler.ID['changeModel'])
      svalue(.$model) = paste(1:length(mn),collapse=" ")
      unblockhandler(.$model,.$handler.ID['changeModel'])
    }

    # Check model definition
    m = .$getModel()
    if(length(m) != K) {
      gmessage(.$translate("Incorrect number of groups in model definition."))
      return()
    }
    
   .ws$setStatus(.$translate("Estimation in progress..."))
    svalue(.$bestOne) = ""

    # Generate all possible model definitions
    if(K==2) models = cbind(c(1,1),c(1,2))
    else     models = setparts(K)
    nmodels = ncol(models)
    svalue(.$possible) = paste(nmodels)

    # Get model prior prob
    priorprob = eval(parse(text=svalue(.$priorprob)))
    m0   = .$testModel(mn,sd,N,rep(1,K))
    msat = .$testModel(mn,sd,N,1:K)

    if(svalue(.$testAll)) {
    
      # Actually not used by BIC approx. but just for displaying
      blockhandler(.$priorprob,.$handler.ID['changePriorProb'])
      svalue(.$priorprob) = paste("1/",nmodels,sep="")
      unblockhandler(.$priorprob,.$handler.ID['changePriorProb'])

      # Compute all possible models
      results = matrix(0,0,3*K)
      bics = vector()
      for(i in 1:ncol(models)) {
        test = .$testModel(mn,sd,N,models[,i])
        results = rbind(results,c(test$means,test$IC.inf,test$IC.sup))
        bics = c(bics,test$bic)
        rownames(results)[i] = test$name
      }
      colnames(results) = c(paste("m",1:K,sep=""),paste("Inf",1:K,sep=""),paste("Sup",1:K,sep=""))

      # Compute approximate model posterior probs
      modProbs = exp(-0.5*(bics - min(bics)))
      modProbs = modProbs / sum(modProbs)
      
      # Store model details
     .$modelDetails = data.frame(Model=row.names(results),BIC=bics, Probs=round(modProbs,4))

      # Bayesian Model Averaging
      if(svalue(.$doBMA)) {
      
        # A posteriori model averaged estimations (Neath & Cavanaugh, 2006)
        post.estim = colSums(results*modProbs)
        target = list(name=.$translate("Averaged"),means=post.estim[1:K],IC.inf=post.estim[(K+1):(2*K)],IC.sup=post.estim[(2*K+1):(3*K)])
        svalue(.$bestOne) = .$translate("Averaged")
        svalue(.$bf) = ""
        svalue(.$postprob) = ""
      }
      
      # Best model selection
      else {

        best.one = which.max(modProbs)
        target = .$testModel(mn,sd,N,models[,best.one])
        svalue(.$bestOne) = target$name
        svalue(.$bf) = round(target$bic,4)
        svalue(.$postprob) = round(modProbs[best.one],4)
      }

     .$postestim[,] = matrix("",10,4)
     .$postestim[1:K,1] = paste(1:K)
     .$postestim[1:K,2] = round(target$means,4)
     .$postestim[1:K,3] = round(target$IC.inf,4)
     .$postestim[1:K,4] = round(target$IC.sup,4)
    }
    
    # Theoretical model
    else {
    
      # Get prior prob. on the target model
      if(is.null(priorprob)) {
        blockhandler(.$priorprob,.$handler.ID['changePriorProb'])
        svalue(.$priorprob) = "1/2"
        unblockhandler(.$priorprob,.$handler.ID['changePriorProb'])
        priorprob = 1/2
      }
      target = .$testModel(mn,sd,N,m)

      # Compute posterior prob
      target.bf = exp(-.5*(target$bic-m0$bic))
      post.prob = priorprob * target.bf/(priorprob * target.bf + 1 - priorprob)
      
      svalue(.$bf) = paste(round(target$bic,4))
      svalue(.$postprob) = paste(round(post.prob,4))

      # A posteriori model averaged estimations (Stein, 1955)
      post.estim = post.prob*target$means + (1-post.prob)*m0$means

     .$postestim[,] = matrix("",10,4)
     .$postestim[1:K,1] = paste(1:K)
     .$postestim[1:K,2] = round(target$means,4)
     .$postestim[1:K,3] = round(target$IC.inf,4)
     .$postestim[1:K,4] = round(target$IC.sup,4)
    }
    
    # Plot
    ymin = min(msat$means) - 3*msat$error/sqrt(min(N))
    ymax = max(msat$means) + 3*msat$error/sqrt(min(N))
    plot(1:K,msat$means,xlab=.$translate("Group"),ylab=.$translate("Estimated means"),main="",type="n",xaxt="n",ylim=c(ymin,ymax))
    axis(1,1:K,paste(1:K))
    abline(h=m0$means[1],col="lightgrey",lty=2,lwd=2)
    
    points(1:K,msat$means,cex=1.5)
    for(k in 1:K) lines(rbind(c(k,target$IC.inf[k]),c(k,target$IC.sup[k])),col="red")
    points(1:K,target$means,pch=19,col="red")
    legend("topleft",pch=c(1,19),pt.cex=1.2,col=c("black","red"),legend=.$translate(c("Data","Model")),inset=.01)
    
   .ws$setStatus(.$translate("Ready."))
  },
  getModel = function(.) {
  
    m = unlist(strsplit(svalue(.$model)," "))
    m[m != ""]
  },
  testModel = function(.,m0,s0,N0,mod,conf=.95) {
  
    # Dimension
    N = sum(N0)
    K = length(unique(mod))

    # Sums
    sumx = N0*m0
    sumx2 = (N0-1)*s0**2
    
    # Constrained means
    Nt  = tapply(N0,mod,sum)
    mt  = tapply(sumx,mod,sum)/Nt

    # new sum of squares wrt the new means mt
    sumx2 = sumx2 + N0*(m0-mt[mod])**2
    s2t = sum(sumx2)/N
    corr.s2t = (N/(N-K)) * s2t
    
    # Credibility interval
    IC.inf = mt + qt((1-conf)/2,N-K)*sqrt(corr.s2t/Nt)
    IC.sup = mt + qt((1+conf)/2,N-K)*sqrt(corr.s2t/Nt)
    
    # BIC
    t = length(mt) + 1
    bic = N + N*log(2*pi*s2t) + t*log(N)
    prob = exp(-.5 * bic)
    
    # Model name
    group.names=paste(1:length(m0))
    model.name = tapply(group.names,mod,paste,collapse=",")
    model.name = paste(sapply(model.name,function(x) paste("(",x,")",sep="")),collapse=",")
    
    list(model=mod,name=model.name,means=mt[mod],error=sqrt(corr.s2t),bic=bic,prob=prob,IC.inf=IC.inf[mod],IC.sup=IC.sup[mod])
  },

  ### Print model detais (formula, BIC, Post. prob)
  printModels = function(.,...) {
  
    k = order(.$modelDetails$BIC)
    cat(.$translate("Model ranking:"),"\n")
    print(.$modelDetails[k,],row.names=FALSE)
    cat("\n")
    galert(.$translate("Model ranking has been printed to the console."))
  },
  
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },

  #---------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                                    CONTENT
  #---------------------------------------------------------------------------------------
  handler.ID   = list(),               # IDs of various handlers that may have to be blocked
  means        =  NULL,                #
  sds          =  NULL,                #
  Ntot         =  NULL,                #
  priorprob    =  NULL,                #
  model        =  NULL,                #
  possible     =  NULL,                #
  testAll      =  NULL,                #
  bestOne      =  NULL,                #
  bf           =  NULL,                #
  postprob     =  NULL,                #
  postestim    =  NULL,                #
  modelDetails = data.frame(),         # Stores a list of all possible models with BICs and approximate post. probs.
  detailWindow = NULL
)


