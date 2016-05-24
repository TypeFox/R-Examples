#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                         Distribution of a sample variance
#-----------------------------------------------------------------------------------------

.ws5 = proto(

  create = function(.,h,...) {

    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Distribution of\na sample variance") %in% names(.ws$nb)) return()
    
   .$sampleSize   = gradio(c(10,20,50,100),handler=.$updatePlot,coerce.with=as.numeric)
   .$nvar         = gradio(c(5,10,50,500),handler=.$updatePlot,coerce.with=as.numeric)
   .$compType     = gradio(.$translate(c("With true mean","With sample mean","Corrected variance")),handler=.$updatePlot)
   .$displayNorm  = gcheckbox(.$translate("gaussian dist."),handler=.$updatePlot)
   .$popDistr     = gcheckbox(.$translate("Chi-2 dist."),handler=.$updatePlot)
   .$param1       = gedit("0",width=5,coerce.with=as.numeric)
   .$param2       = gedit("1",width=5,coerce.with=as.numeric)
   .$var          = glabel("")
   .$truevar      = glabel("")

    # Construction de l'interface
    add(.ws$nb,group <- ggroup(horizontal=FALSE),label=.$translate("Distribution of\na sample variance"))

    tmp = gframe(.$translate("Gaussian distribution of scores"), container = group)
    distribGroup = glayout(container=tmp)
    distribGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Mean"))
    distribGroup[2,3]=.$param1
    distribGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Standard deviation"))
    distribGroup[3,3]=.$param2
    visible(distribGroup)=TRUE

    tmp = gframe(.$translate("Variance computation"), container = group, horizontal=FALSE)
    add(tmp,.$compType)

    sizeGroup = ggroup(cont = group,expand=TRUE)
    tmp = gframe(.$translate("Samples"), container = sizeGroup,expand=TRUE)
    add(tmp, .$nvar)
    tmp = gframe(.$translate("Obs. by sample"), container = sizeGroup,expand=TRUE)
    add(tmp, .$sampleSize)

    # Options d'affichage
    tmp = gframe(.$translate("Display"), container = group, horizontal=FALSE)
    add(tmp,.$popDistr)
    add(tmp,.$displayNorm)

    # Statistiques descriptives
    tmp = gframe(.$translate("Descriptive statistics"), container = group, horizontal=FALSE)
    resGroup = glayout(container=tmp)
    resGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Expected mean"))
    resGroup[2,3]=.$truevar
    resGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Observed mean"))
    resGroup[3,3]=.$var
    visible(resGroup)=TRUE

    addSpring(group)

    buttonGroup=ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Sample"),container=buttonGroup, handler=.$updatePlot)
  },
  
  updatePlot = function(.,h,...) {

    # Vérification des paramètres
    if(any(is.na(c(svalue(.$param1),svalue(.$param2))))) {
      gmessage("Please provide parameter values.")
      return()
    }
    
    nobs     = as.numeric(svalue(.$sampleSize))
    nsamp    = as.numeric(svalue(.$nvar))
    compType = svalue(.$compType,index=T)
    plotNorm = svalue(.$displayNorm)
    plotChi2 = svalue(.$popDistr)
    param1   = svalue(.$param1)
    param2   = svalue(.$param2)
    
    # Paramètres vrais
    y = matrix(rnorm(nobs*nsamp,param1,param2),nobs,nsamp)
    true.mean = param1
    true.var = param2**2
    
    # Statistiques d'échantillon
    xbar = colMeans(y)
    
	  # Moyenne vraie
    if(compType==1) { 
      sample.vars = colMeans((y-true.mean)**2)
      nu = nobs
      xlab=expression(S[mu]**2)
    }
	  # Moyenne d'échantillon
    else if(compType==2) {
      sample.vars = rowMeans((t(y)-xbar)**2)
      nu = nobs-1
      xlab=expression(S[bar(X)]**2)
    }
	  # Moyenne corrigée
    else {
      sample.vars = (nobs/(nobs-1))*rowMeans((t(y)-xbar)**2)
      nu = nobs-1
      xlab=expression(S[bar(X)]**2)
    }

    # Affichage graphique  
    title = .$translate("Distribution of sample variances")
    ylab = .$translate("Densities")

    avg.var = mean(sample.vars)
    xlim = c(0,max(sample.vars))
    q = seq(xlim[1],xlim[2],len=100)*nu/true.var
    hist(sample.vars,freq=FALSE,main=title,xlab=xlab,ylab=ylab,xlim=xlim)
    abline(v=true.var,col="green",lwd=2)
    abline(v=avg.var,col="green",lwd=2,lty=2)

    # Changement d'échelle pour s'ajuster à la distribution des variances
    if(plotChi2) lines(q*true.var/nu,nu*dchisq(q,nu)/true.var,lwd=2,col="red")
    if(plotNorm) lines(q*true.var/nu,nu*dnorm(q,nu,sqrt(2*nu))/true.var,lwd=2,col="blue")
    
    # Affichage des stats descriptives
    svalue(.$truevar) = paste(true.var)
    svalue(.$var) = paste(round(avg.var,3))

  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },

  #---------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                                    CONTENT
  #---------------------------------------------------------------------------------------
  sampleSize     =  NULL,                       #
  nvar           =  NULL,                       #
  compType       =  NULL,                       #
  displayNorm    =  NULL,                       #
  popDistr       =  NULL,                       #
  param1         =  NULL,                       #
  param2         =  NULL,                       #
  var            =  NULL,                       #
  truevar        =  NULL                       #
)


