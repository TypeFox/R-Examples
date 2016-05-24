#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                    Construction of the gaussian distribution
#-----------------------------------------------------------------------------------------

.ws1 = proto(

  create = function(.,h,...) {
  
    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Construction of the\ngaussian distribution") %in% names(.ws$nb)) return()
    
    names(.$availDists) = .$translate(names(availDists))

   .$paramNames$Binomial = .$translate(paramNames$Binomial)
   .$paramNames$Uniform  = .$translate(paramNames$Uniform)
   .$paramNames$Gaussian = .$translate(paramNames$Gaussian)
   .$paramNames$Gamma    = .$translate(paramNames$Gamma)
   
    names(.$paramNames)  = .$translate(names(paramNames))

   .$distribution = gdroplist(names(.$availDists),horizontal=FALSE,handler=.$updateParamNames)
   .$sampleSize = gradio(c(500, 1000, 5000, 50000),handler=.$updatePlot,coerce.with=as.numeric)
   .$nvar = gradio(c(1, 2, 10, 50),handler=.$updatePlot,coerce.with=as.numeric)
   .$displayWhat = gradio(.$translate(c("Counts","Proportions","Densities")),handler=.$updatePlot)
   .$displayFunc = gcheckbox(.$translate("theoretical"),handler=.$updatePlot)
   .$displayNorm = gcheckbox(.$translate("gaussian"),handler=.$updatePlot)
   .$param1 = gedit("0",width=5,coerce.with=as.numeric)
   .$paramLabel1 = glabel(.$translate("Left boundary"))
   .$param2 = gedit("1",width=5,coerce.with=as.numeric)
   .$paramLabel2 = glabel(.$translate("Right boundary"))
   .$cutpoints = gedit("",handler=.$updatePlot)

    add(.ws$nb, group <- ggroup(horizontal=FALSE),label=.$translate("Construction of the\ngaussian distribution"))
    tmp = gframe(.$translate("Distribution"), container = group)
    distribGroup = glayout(container=tmp)
    distribGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Family"))
    distribGroup[2,3]=.$distribution
    distribGroup[3,2,anchor=c(-1,0)]=.$paramLabel1
    distribGroup[3,3]=.$param1
    distribGroup[4,2,anchor=c(-1,0)]=.$paramLabel2
    distribGroup[4,3]=.$param2
    visible( distribGroup)=TRUE

    sizeGroup = ggroup(cont = group,expand=TRUE)
    tmp = gframe(.$translate("Sample size"), container =  sizeGroup,expand=TRUE)
    add(tmp, .$sampleSize)
    tmp = gframe(.$translate("Variables"), container =  sizeGroup,expand=TRUE)
    add(tmp, .$nvar)

    histGroup = ggroup(cont = group,expand=TRUE)
    tmp = gframe(.$translate("Display"), container = histGroup,expand=TRUE)
    add(tmp,.$displayWhat)
    tmp = gframe(.$translate("Family"), container = histGroup, horizontal=FALSE,expand=TRUE)
    add(tmp,.$displayFunc)
    add(tmp,.$displayNorm)

    tmp = gframe(.$translate("Cutpoints"), container = group)
    add(tmp,.$cutpoints,expand=TRUE)

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
    
    rfunc = paste("r",.$availDists[svalue(.$distribution)],sep="")
    dfunc = paste("d",.$availDists[svalue(.$distribution)],sep="")
    
    # Génération des données
    nvar = as.numeric(svalue(.$nvar))
    nobs = as.numeric(svalue(.$sampleSize))
    
    y = do.call(rfunc, list(nobs*nvar,svalue(.$param1),svalue(.$param2)))
    x = rowSums(matrix(y,nobs,nvar))
    
    # Définition des coupures
    if(nchar(svalue(.$cutpoints))) { 
      breaks = unlist(strsplit(svalue(.$cutpoints)," "))
      breaks = breaks[breaks!=""]
      if(!length(breaks)) breaks="sturges"
      else breaks = as.numeric(breaks)
    } 
    else { breaks = "sturges" }

    # Affichage de l'histogramme empirique (distributions continues)
    if(svalue(.$distribution)!=.$translate("Binomial")) { 
      hh = hist(x,breaks=breaks,plot=FALSE)
      if(svalue(.$displayWhat)==.$translate("Proportions")) {
        hh$counts = hh$counts / nobs
      }
      xlab = ifelse(nvar==1,.$translate("Variable values"),.$translate("Values of summed variables"))
      plot(hh,freq=svalue(.$displayWhat)!=.$translate("Densities"),main = paste(.$translate("Distribution"),":",svalue(.$distribution)),xlab=xlab,ylab=svalue(.$displayWhat))

      # Affichage de la loi théorique
      if(svalue(.$displayFunc) && (nvar==1)) {
        z = seq(min(x),max(x),len=100)
        lines(z,do.call(dfunc, list(z,svalue(.$param1),svalue(.$param2))),lwd=2,col="red")
      }
    }
    
    # Loi binomiale
    else {
      # Affichage de l'histogramme empirique
      nn = table(x)
      if(svalue(.$displayWhat) %in% .$translate(c("Proportions","Densities"))) {
        nn = nn / nobs
      }
      xlab = ifelse(svalue(nvar)==1,.$translate("Variable values"),.$translate("Values of summed variables"))
      res = plot(nn,main = paste(.$translate("Distribution"),svalue(.$distribution)),xlab=xlab,ylab=svalue(.$displayWhat))
	
	  # Affichage des probabilités théoriques
      if(svalue(.$displayFunc) && (nvar==1)) {
        z = min(x):max(x)
	    lines(z+.1,dbinom(z,svalue(.$param1),svalue(.$param2)),col="red",lwd=2,type="h")
      }
	
    }

    if(svalue(.$displayNorm)) {
      z = seq(min(x),max(x),len=100)
      lines(z,dnorm(z,mean(x),sd(x)),lwd=2,col="blue")
    }

  },
  
  updateParamNames = function(.,h,...) {
  
    newDist = svalue(.$distribution)
    
    # Warning: a droplist may be temporarily set to NULL when changed
    if(is.null(newDist)) return()
    
    svalue(.$paramLabel1) = .$paramNames[[newDist]][1]
    svalue(.$paramLabel2) = .$paramNames[[newDist]][2]
    
  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },
  #---------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                                    CONTENT
  #---------------------------------------------------------------------------------------
  availDists     = c(Uniform = "unif", Binomial = "binom", Gaussian = "norm", Gamma = "gamma"),
  paramNames     = list(Uniform=c("Left boundary","Right boundary"),Binomial=c("Size","Probability"),Gaussian=c("Mean","Standard dev."),Gamma=c("Shape","Scale")),
  distribution   = NULL,
  sampleSize     = NULL,
  nvar           = NULL,
  displayWhat    = NULL,
  displayFunc    = NULL,
  displayNorm    = NULL,
  param1         = NULL,
  paramLabel1    = NULL,
  param2         = NULL,
  paramLabel2    = NULL,
  cutpoints      = NULL
)


