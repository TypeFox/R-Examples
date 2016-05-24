#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                     Bayesian inference on several proportions
#-----------------------------------------------------------------------------------------

.ws10 = proto(

  create = function(.,h,...) {

    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Bayesian inference\non a contingency table") %in% names(.ws$nb)) return()
    group <- ggroup(horizontal=FALSE,cont=.ws$nb,label=.$translate("Bayesian inference\non a contingency table"))
    
    tmp <- gframe(.$translate("Observed data"),cont=group,expand=TRUE)
   .$counts <- gtext("",cont=tmp,expand=TRUE,font.attr=c(family="monospace"))

   .$alpha = gedit("1",width=3,coerce.with=as.numeric,handler=.$compute)

   .$priorprob = gedit("",width=5)
   .$handler.ID['changePriorProb'] = addhandlerchanged(.$priorprob,.$compute)

   .$model = gedit("",width=15)
   .$handler.ID['changeModel'] = addhandlerchanged(.$model,.$compute)

   .$bf = glabel("")
   .$possible = glabel("")
   .$postprob = glabel("")
   .$testAll = gcheckbox(.$translate("Test all"),checked=FALSE,handler=.$onTestAll)
   .$bestOne = glabel("")

    add(group,tmp <- gframe(.$translate("Model validation")))
    testGroup = glayout(container=tmp)
    testGroup[2,2,anchor=c(-1,0)] = glabel(.$translate("Prior parameter"))
    testGroup[2,3] = .$alpha
    testGroup[3,2,anchor=c(-1,0)] = glabel(.$translate("Possible models"))
    testGroup[3,3] = .$possible
    testGroup[4,2,anchor=c(-1,0)] = glabel(.$translate("Target model"))
    testGroup[4,3] = .$model
    testGroup[5,2,anchor=c(-1,0)] = glabel(.$translate("Prior Pr(M)"))
    testGroup[5,3] =.$priorprob
    testGroup[6,3] =.$testAll
    testGroup[7,2] = glabel(.$translate("Best model"))
    testGroup[7,3] =.$bestOne
    testGroup[8,2] = glabel(.$translate("Bayes factor"))
    testGroup[8,3] =.$bf
    testGroup[9,2] = glabel(.$translate("Pr(M|D)"))
    testGroup[9,3] =.$postprob
    testGroup[10,1] = ""
    visible(testGroup)=TRUE

    add(group, tmp <- gframe(.$translate("Model posterior estimates")),expand=TRUE)
    add(tmp, .$estimNb <- gnotebook(), expand=TRUE)
    add(.$estimNb, tmp <- ggroup(),label=.$translate("Model"),expand=TRUE)
    add(tmp,.$postestim <- gtext("",font.attr=c(family="monospace")),expand=TRUE)
    add(.$estimNb, tmp <- ggroup(),label=.$translate("Averaged"),expand=TRUE)
    add(tmp,.$avestim <- gtext("",font.attr=c(family="monospace")),expand=TRUE)
    svalue(.$estimNb,index=TRUE) = 1

    buttonGroup=ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Compute"),container=buttonGroup, handler=.$compute)
  
    # A associer en dernier pour ne pas déclencher pendant la construction
    addHandlerChanged(.$estimNb,.$updatePlot)
    addHandlerChanged(.$testAll,.$compute)

  },
  compute = function(.,h,...) {

    # Check parameters
    if(is.na(svalue(.$alpha))) {
      gmessage(.$translate("Please specify a parameter value for the Dirichlet prior."))
      return()
    }
    a = svalue(.$alpha)
    
    # Check data
   .$getData()
    if(is.null(.$n)) {
      gmessage(.$translate("Please specify observed data."))
      return()
    }
    
    if(any(.$n<0)) {
      gmessage(.$translate("Error: All rows in the table should\ncontain the same number of values."))
      return()
    }

    I = nrow(.$n)
    if(I<2) {
      gmessage(.$translate("Please provide at least two count rows."))
      return()
    }

    # Saturated model by default
    if(!nchar(svalue(.$model))) {
      blockhandler(.$model,.$handler.ID['changeModel'])
      svalue(.$model) = paste(1:I,collapse=" ")
      unblockhandler(.$model,.$handler.ID['changeModel'])
    }
    
    m = .$getModel()
    if( (length(m) != I) && !svalue(.$testAll)) {
      gmessage(.$translate("Incorrect number of groups in model definition."))
      return()
    }
    
    # Number of possible models
    if(I==2) models = cbind(c(1,1),c(1,2))
    else     models = setparts(I)
    nmodels = ncol(models)
    svalue(.$possible) = paste(nmodels)
    
   .ws$setStatus(.$translate("Estimation in progress..."))
    svalue(.$bestOne) = ""

    # Reference models (null and saturated)
   .$p0 = .$testModel(rep("1",I),a)$estimates
   .$ps = .$testModel(1:I,a)$estimates

    # Test all models
    if(svalue(.$testAll)) {
    
      # Equal prior probs in this case
      blockhandler(.$priorprob,.$handler.ID['changePriorProb'])
      svalue(.$priorprob) = paste("1/",nmodels,sep="")
      unblockhandler(.$priorprob,.$handler.ID['changePriorProb'])
      priorprob = 1/nmodels

     .$avestimates = matrix(0,nrow(.$n),ncol(.$n))
      group.names = paste(1:I)
      sum.bf = 0
      best.bf = -1
      best.model = ""
            
      for(i in 1:nmodels) {
        svalue(.$possible) = paste(i,"/",nmodels)
        Sys.sleep(.001)
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
      svalue(.$postestim) = capture.output(round(.$p1,3))
      svalue(.$avestim)   = capture.output(round(.$avestimates,3))
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
      svalue(.$postestim) = capture.output(round(.$p1,3))
      svalue(.$avestim)   = capture.output(round(.$avestimates,3))
    }
    
   .$updatePlot(h,...)
   .ws$setStatus(.$translate("Ready."))
   
  },
  
  updatePlot = function(.,h,...) {
      
    # Graphique de base : modèle homogène
    C = ncol(.$n)
    responses = paste(1:C)
    matplot(t(.$p0),type="l",col="grey",lwd=2,xlab=.$translate("Response"),ylab=.$translate("Estimated probabilities"),ylim=c(0,1),main="",xaxt="n")
    axis(1,1:C,responses)

    # Estimation modèle
    estimType = if(is.null(h$pageno)) svalue(.$estimNb) else h$pageno
    if(estimType == 1) {
      matlines(t(.$p1),col="red",lwd=2,lty=.$groups)
      legend("topleft",lty=1,lwd=2,col=c("red","gray"),legend=c(ifelse(svalue(.$testAll),.$translate("Best"),.$translate("Target")),.$translate("Homogeneous")),inset=.01)
    }
    else {
      matlines(t(.$avestimates),lwd=2,col="blue",lty=.$groups)
      legend("topleft",lty=1,lwd=2,col=c("blue","gray"),legend=.$translate(c("Averaged","Homogeneous")),inset=.01)
    }
    
    # Saturé
    matpoints(t(.$ps),type="p",pch=21,bg="white",col=1,cex=1.5)
    matpoints(t(.$ps),col=1,cex=.5)
  },
  
  getData = function(.,h,...) {
  
    data = svalue(.$counts)
    for(ch in c(",",";","\t")) data = gsub(ch," ",data)
    
    # Strip double spaces
    while(length(grep("  ",data))) data = gsub("  "," ",data)
    
    # Read lines
    rows = unlist(strsplit(data,split="\n"))
    rows = rows[rows != ""]
    
    if(!length(rows)) return(-1)
    
    rows = sapply(rows,strsplit,split=" ")
    l = sapply(rows,length)

    if(var(l)) return(-1)

   .$n = do.call("rbind",lapply(rows,as.numeric))
    row.names(.$n) = paste("l",1:nrow(.$n),sep="")
    
  },
  getModel = function(.) {
  
    m = unlist(strsplit(svalue(.$model)," "))
    m[m != ""]
  },

  testModel = function(.,mod,a) {
  
    stopifnot(is.matrix(.$n))
    
    C  = ncol(.$n)
    Kj = colSums(.$n)
    n2 = rowsum(.$n,mod)
    A  = matrix(a,nrow(.$n),ncol(.$n))
    A2 = matrix(a,nrow(n2),ncol(n2))
    
    lBeta = function(z)     rowSums(lgamma(z)) - lgamma(rowSums(z))
    lPsi  = function(nn,aa) lBeta(nn+aa) - lBeta(aa)
    
    lnum = sum(lPsi(n2,A2))
    lden = lPsi(matrix(Kj,ncol=C),matrix(a,1,C))

    bf = exp(lnum - lden)
    estimates = n2[mod,] + A
    estimates = estimates/rowSums(estimates)
    row.names(estimates) = paste(1:nrow(estimates)," ")

    list(estimates=estimates,bf=bf)
  },

  onTestAll = function(.,h,...) {
  
    if(svalue(.$testAll)) enabled(.$model) = FALSE
    else                  enabled(.$model) = TRUE    
  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },
  #---------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                                    CONTENT
  #---------------------------------------------------------------------------------------
  handler.ID   = list(),               # IDs of various handlers that may have to be blocked
  alpha        =  NULL,                #
  counts       =  NULL,                #
  n            =  NULL,                #
  priorprob    =  NULL,                #
  model        =  NULL,                #
  m            =  NULL,                #
  possible     =  NULL,                #
  testAll      =  NULL,                #
  bestOne      =  NULL,                #
  bf           =  NULL,                #
  postprob     =  NULL,                #
  estimNb      =  NULL,                #
  postestim    =  NULL,                #
  avestim      =  NULL,                #
  p0           =  NULL,                #
  p1           =  NULL,                #
  ps           =  NULL,                #
  avestimates  =  NULL                 #
)


