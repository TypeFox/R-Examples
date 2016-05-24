#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                Effect of change of origin and scale on the gaussian
#-----------------------------------------------------------------------------------------

.ws2 = proto(

  create = function(.,h,...) {
   
    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Change of origin\nand scale") %in% names(.ws$nb)) return()
    
   .$nobs = gradio(c(50, 500, 50000),handler=.$updatePlot,coerce.with=as.numeric)
   .$param1 = gedit("100",width=5,coerce.with=as.numeric)
   .$param2 = gedit("15",width=5,coerce.with=as.numeric)
   .$add = gedit("0",width=5)
   .$mult = gedit("1",width=5)
   .$standard = gcheckbox(.$translate("Standardize the sample"),handler=.$updatePlot)
   .$xbar = glabel("")
   .$s = glabel("")

    # Construction de l'interface
    add(.ws$nb,group <- ggroup(horizontal=FALSE),label=.$translate("Change of origin\nand scale"))

    # Distribution
    tmp = gframe(.$translate("Population parameters"), container = group)
    distribGroup = glayout(container= tmp)
    distribGroup[2,2:3]=glabel(.$translate("Gaussian distribution"))
    distribGroup[3,2,anchor=c(-1,0)]=glabel(.$translate("Mean"))
    distribGroup[3,3]=.$param1
    distribGroup[4,2,anchor=c(-1,0)]=glabel(.$translate("Standard dev."))
    distribGroup[4,3]=.$param2
    visible(distribGroup)=TRUE

    # Effectifs          
    tmp = gframe(.$translate("Sample size"), container = group)
    add( tmp,.$nobs)

    # Transformation
    tmp = gframe(.$translate("Transformation"), container = group, horizontal=FALSE)
    transfoGroup = glayout(container= tmp)
    transfoGroup[2,2:3]=glabel(" X' = aX + b ")
    transfoGroup[3,2,anchor=c(-1,0)]=glabel(" a = ")
    transfoGroup[3,3,expand=TRUE]=.$mult
    transfoGroup[4,2,anchor=c(-1,0)]=glabel(" b = ")
    transfoGroup[4,3,expand=TRUE]=.$add
    visible(transfoGroup)=TRUE
    add( tmp,.$standard)

    # Statistiques descriptives
    tmp = gframe(.$translate("Descriptive statistics"), container = group, horizontal=FALSE)
    resultGroup = glayout(container= tmp)
    resultGroup[2,2] = glabel(.$translate("Mean"))
    resultGroup[2,3] = .$xbar
    resultGroup[3,2] = glabel(.$translate("Standard deviation"))
    resultGroup[3,3] = .$s

    addSpring(group)

    # Boutons de commande
    buttonGroup=ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Sample"),container=buttonGroup, handler=.$updatePlot)

  },
  
  updatePlot = function(.,h,...) {
  
    # Vérification des paramètres
    if(any(is.na(c(svalue(.$param1),svalue(.$param2))))) {
      gmessage(.$translate("Please provide parameter values."))
      return()
    }

    if(is.na(svalue(.$add)))  svalue(.$add)  = 0
    if(is.na(svalue(.$mult))) svalue(.$mult) = 1

    nobs = svalue(.$nobs)
    param1 = svalue(.$param1)
    param2 = svalue(.$param2)
    add = eval(parse(text=svalue(.$add)))
    mult = eval(parse(text=svalue(.$mult)))

    # Génération des données
    y = rnorm(nobs,param1,param2)
    m1 = mean(y)
    sd1 = sd(y)
    
    # Transformation des données
    if(svalue(.$standard)) y = as.numeric(scale(y))
    y = y * mult + add
    m2 = mean(y)
    sd2 = sd(y)
    
    # Affichage des stats descriptives
    svalue(.$xbar)=paste(round(m2,3))
    svalue(.$s)=paste(round(sd2,3))

    # Paramètres transformés de population
    new.mu = param1*mult+add
    new.sigma = param2*abs(mult)
    
    # Représentation graphique
    hist(y,xlab=.$translate("Variable"),ylab=.$translate("Density"),main=.$translate("Gaussian distribution"),freq=FALSE)
    if(svalue(.$standard)) curve(dnorm(x,0,1),from=-4,to=4,add=TRUE,lwd=2,col="blue")
    else curve(dnorm(x,param1*mult+add,param2*abs(mult)),from=min(y),to=max(y),add=TRUE,lwd=2,col="blue")
    
    # Afficher moyenne vraie transformée en rouge
    abline(v=new.mu,lwd=2,col="red")
    
    # Afficher écart-type vrai transformé en vert
    lines(rbind(c(new.mu,dnorm(new.mu+new.sigma,new.mu,new.sigma)),c(new.mu+new.sigma,dnorm(new.mu+new.sigma,new.mu,new.sigma))),col="green",lwd=2)

  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },
  #---------------------------------------------------------------------------------------
  #  SLOT         INITIAL VALUE                  CONTENT
  #---------------------------------------------------------------------------------------
  nobs            = NULL,                    #
  param1          = NULL,                    #
  param2          = NULL,                    #
  add             = NULL,                    #
  mult            = NULL,                    #
  standard        = NULL,                    #
  xbar            = NULL,                    #
  s               = NULL                     #
)


