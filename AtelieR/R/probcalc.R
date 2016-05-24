#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                           Probability calculator
#-----------------------------------------------------------------------------------------

.ws4 = proto(

  create = function(.,h,...) {

    # Don't create if no main interface exists
    if(inherits(try(is.environment(.ws)),"try-error")) return()
    
    # Don't create if already opened
    if(.$translate("Probability\ncalculator") %in% names(.ws$nb)) return()
    
    names(.$availDists) = .$translate(names(availDists))
    
   .$paramNames$Uniform               = .$translate(.$paramNames$Uniform)
   .$paramNames$Binomial              = .$translate(.$paramNames$Binomial)
   .$paramNames$Gaussian              = .$translate(.$paramNames$Gaussian)
   .$paramNames$Gamma                 = .$translate(.$paramNames$Gamma)
   .$paramNames$"Inverse Gamma"       = .$translate(.$paramNames$"Inverse Gamma")
   .$paramNames$"Chi-2"               = .$translate(.$paramNames$"Chi-2")
   .$paramNames$"Inverse Chi-2"       = .$translate(.$paramNames$"Inverse Chi-2")
   .$paramNames$"Non central Chi-2"   = .$translate(.$paramNames$"Non central Chi-2")
   .$paramNames$Beta                  = .$translate(.$paramNames$Beta)
   .$paramNames$Poisson               = .$translate(.$paramNames$Poisson)
   .$paramNames$Student               = .$translate(.$paramNames$Student)
   .$paramNames$"Generalized Student" = .$translate(.$paramNames$"Generalized Student")
   .$paramNames$"Non central Student" = .$translate(.$paramNames$"Non central Student")
   .$paramNames$Fisher                = .$translate(.$paramNames$Fisher)
   .$paramNames$"Non central Fisher"  = .$translate(.$paramNames$"Non central Fisher")
   .$paramNames$"Lambda prime"        = .$translate(.$paramNames$"Lambda prime")

    names(.$paramNames) = .$translate(names(paramNames))

   .$distribution = gdroplist(names(.$availDists),horizontal=FALSE,handler=.$initOptions)
   .$param1 = gedit("0",width=15)
   .$param2 = gedit("1",width=15)
   .$param3 = gedit("",width=15)
    enabled(.$param3) = FALSE
   .$paramLabel1 = glabel(.$translate("Mean"))
   .$paramLabel2 = glabel(.$translate("Standard dev."))
   .$paramLabel3 = glabel("")
   
    # Construction de l'interface
    add(.ws$nb, group <- ggroup(horizontal=FALSE),label=.$translate("Probability\ncalculator"))

    tmp = gframe(.$translate("Distribution"),container=group)
    distribGroup = glayout(container=tmp)
    distribGroup[2,2,anchor=c(-1,0)]=glabel(.$translate("Family"))
    distribGroup[2,3]=.$distribution
    distribGroup[3,2,anchor=c(-1,0)]=.$paramLabel1
    distribGroup[3,3]=.$param1
    distribGroup[4,2,anchor=c(-1,0)]=.$paramLabel2
    distribGroup[4,3]=.$param2
    distribGroup[5,2,anchor=c(-1,0)]=.$paramLabel3
    distribGroup[5,3]=.$param3
    visible(distribGroup)=TRUE

   .$calcWhat = gradio(.$translate(c("Quantile => probability","Probability => quantile")),handler=.$updatePlot)
    tmp = gframe(.$translate("Computation"),container=group)
    add(tmp,.$calcWhat)

   .$side = gradio(.$translate(c("Left-sided","Right-sided")),handler=.$updatePlot)
    tmp = gframe(.$translate("Cumulative probability"),container=group)
    add(tmp,.$side)

   .$value  = gedit(width=15,handler=.$updatePlot)
   .$resultx = glabel("")
   .$resultp = glabel("")

    tmp = gframe(.$translate("Value or expression to compute"),container=group)
    add(tmp,.$value,expand=TRUE)

    tmp = gframe(.$translate("Result"),container=group)
    resultGroup = glayout(container=tmp)
    resultGroup[2,2] = " x ="
    resultGroup[2,3,expand=TRUE,anchor=c(-1,0)] = .$resultx
    resultGroup[3,2] = " p ="
    resultGroup[3,3,expand=TRUE,anchor=c(-1,0)] = .$resultp

    addSpring(group)

    # buttons
    buttonGroup = ggroup(container=group)
    addSpring(buttonGroup)
    gbutton(.$translate("Compute"),container=buttonGroup, handler=.$updatePlot)
  
  },
  
  updatePlot = function(.,h,...) {

    param1 = eval(parse(text=svalue(.$param1)))
    param2 = eval(parse(text=svalue(.$param2)))
    param3 = eval(parse(text=svalue(.$param3)))
    
    if(is.null(param1)) {
      gmessage(.$translate("Please provide parameter values."))
      return()
    }
    
    distrib = svalue(.$distribution)
    is1P = distrib %in% .$translate(c("Student","Chi-2","Poisson"))
    is2P = distrib %in% .$translate(c("Gaussian","Inverse Chi-2","Non central Chi-2","Fisher","Binomial","Gamma","Inverse Gamma","Beta"))
    is3P = distrib %in% .$translate(c("Generalized Student","Non central Student","Non central Fisher","Lambda prime"))
    
    if(is2P && is.null(param2)) {
        gmessage(.$translate("Some parameter values are missing."))
        return()
    }
      
    if(is3P && is.null(param3)) {
        gmessage(.$translate("Some parameter values are missing."))
        return()
    }
      
    value = eval(parse(text=svalue(.$value)))
    
    if(is.null(value)) {
      gmessage(.$translate("Please fill in the value field."))
      return()
    }

    if( (distrib %in% .$translate(c("Binomial","Poisson","Chi-2","Inverse Chi-2","Fisher","Gamma","Inverse Gamma","Beta")))  && (value < 0))  {
      gmessage(.$translate("This distribution is not defined for negative values."))
      return()  
    }
    
    isDiscrete = distrib %in% c("Binomial","Poisson")
    is01 = function(x) (x>=0)&&(x<=1)
    isInteger = function(x) abs(x)==round(x)
    probf = .$availDists
    
    p = svalue(.$calcWhat,index=T) == 2 # "Probabilité =>  quantile"
    right = svalue(.$side)==.$translate("Right-sided")

    if(p && !is01(value)) {
      gmessage(.$translate("A probability lies between 0 and 1."))
      return()
    }
    
    # Some extra distributions for Bayesian stats
    
    # Non-central Chi-square
    dncchisq = function(x,nu,ncp) dchisq(x,nu,ncp)
    qncchisq = function(p,nu,ncp) qchisq(p,nu,ncp)
    pncchisq = function(q,nu,ncp) pchisq(q,nu,ncp)
    rncchisq = function(n,nu,ncp) rchisq(n,nu,ncp)

    # Non-central F
    dncf = function(x,nu1,nu2,ncp) df(x,nu1,nu2,ncp)
    qncf = function(p,nu1,nu2,ncp) qf(p,nu1,nu2,ncp)
    pncf = function(q,nu1,nu2,ncp) pf(q,nu1,nu2,ncp)
    rncf = function(n,nu1,nu2,ncp) rf(n,nu1,nu2,ncp)
    
    # Inverse Gamma
    dinvgamma = function(z,a,b) exp(a*log(b) - lgamma(a) -(a+1)*log(z) -(b/z))
    qinvgamma = function(p,a,b) ifelse(((1 - p) <= .Machine$double.eps),Inf,1/qgamma(1-p,a,b))
    pinvgamma = function(z,a,b) 1-pgamma(1/z,a,b)
    rinvgamma = function(n,a,b) 1/rgamma(n=n,shape=a,rate=b)
    
    # Scaled Inversed Chi-Square
    dscaledInvChi2 = function(z,nu,s2) dinvgamma(z,nu/2,nu*s2/2)
    pscaledInvChi2 = function(p,nu,s2) pinvgamma(p,nu/2,nu*s2/2)
    qscaledInvChi2 = function(z,nu,s2) qinvgamma(z,nu/2,nu*s2/2)
    rscaledInvChi2 = function(n,nu,s2) rinvgamma(n,nu/2,nu*s2/2)
    
    # Non standard (3P) Student
    dnst = function(x,nu,pos,scale) dt((x-pos)/scale,nu)/scale
    qnst = function(p,nu,pos,scale) qt(p,nu)*scale + pos
    pnst = function(q,nu,pos,scale) pt((q-pos)/scale,nu)
    rnst = function(n,nu,pos,scale) rt(n,nu)*scale + pos
    
    # Non central scaled Student
    dnct = function(x,nu,ncp,scale) dt(x/scale,nu,ncp/scale)/scale
    qnct = function(p,nu,ncp,scale) qt(p,nu,ncp/scale)*scale
    pnct = function(q,nu,ncp,scale) pt(q/scale,nu,ncp/scale)
    rnct = function(n,nu,ncp,scale) rt(n,nu,ncp/scale)*scale
    
    # Lambda-prime distribution (Lecoutre, 1999)
    rlambdaprime = function(n,nu,ncp,scale) rnorm(n,0,scale) + ncp*sqrt(rchisq(n,nu)/nu)
    plambdaprime = function(x,nu,ncp,scale) 1-pt(ncp/scale,nu,x/scale)
    qlambdaprime = function(p,nu,ncp,scale) {

      t = ncp/scale
      k = exp( ((log(2)-log(nu))/2) + lgamma((nu+1)/2) - lgamma(nu/2) )
      M = k*t
      V = 1 + t**2 - M**2

      N = 100
      from = M - 4.5*sqrt(V)
      to   = M + 4.5*sqrt(V)

      for(iter in 1:3) {
      
        range = seq(from,to,len=N)
        index = which.min(abs(p - plambdaprime(range,nu,ncp,scale)))
        from = range[index-1]
        to = range[index+1]
      }
      
      range[index]
    }
    dlambdaprime = function(x,nu,ncp,scale) {
    
      # Numerical derivatives
      dx = .0001
      dF1 = plambdaprime(x-dx/2,nu,ncp,scale)
      dF2 = plambdaprime(x+dx/2,nu,ncp,scale)
      
      (dF2-dF1)/dx
    }

    dfunction = eval(parse(text=paste("d",probf[distrib],sep="")))
    pfunction = eval(parse(text=paste("p",probf[distrib],sep="")))
    qfunction = eval(parse(text=paste("q",probf[distrib],sep="")))
    rfunction = eval(parse(text=paste("r",probf[distrib],sep="")))
    
    # Chosen distribution has two parameters
    if(is2P) {
    
      # Check parameter values
      if(distrib==.$translate("Binomial")) {
        stopifnot(isInteger(param1) && is01(param2)) }
      if(distrib==.$translate("Fisher")) {
        stopifnot(isInteger(param1) && isInteger(param2)) }
      
      # prob. to quantile
      if(p) {
      
        prob = value
        
        # Continuous distribution
        if(!isDiscrete) {
          if(right) value = qfunction(1-prob,param1,param2)
          else      value = qfunction(prob,param1,param2)
        }
        
        # Discrete distribution
        else {
          if(right) value = qfunction(1-prob,param1,param2)
          else      value = qfunction(prob,param1,param2)
        }
      }
      
      # Quantile to prob.
      else {
      
        # Continuous distribution
        if(!isDiscrete) {
          if(right) prob = 1-pfunction(value,param1,param2)
          else      prob = pfunction(value,param1,param2)
        }
        
        # Discrete distribution
        else {
          if(right)  prob = 1-pfunction(value-1,param1,param2)
          else       prob = pfunction(value,param1,param2)
        }
      }
      
      dens = dfunction(value,param1,param2)
    }
    
    # One-parameter distributions
    else if(is1P) {
    
      svalue(.$param2)=""
      
      # Check parameter values
      if(distrib==.$translate("Student")) {
        stopifnot(param1>0)
      }
      if(distrib==.$translate("Chi-2")) {
        stopifnot(param1>0) 
      }
      
      # Prob. to quantile
      if(p) {
      
        prob = value
        
        # Continuous distribution
        if(!isDiscrete) {
          if(right) value = qfunction(1-prob,param1)
          else      value = qfunction(prob,param1)
        }
        
        # Discrete distribution
        else {
          if(right) value = qfunction(1-prob,param1)
          else      value = qfunction(prob,param1)
        }
      }
      
      # Quantile to prob.
      else {
      
        # Continuous distribution
        if(!isDiscrete) {
          if(right) prob = 1-pfunction(value,param1)
          else      prob = pfunction(value,param1)
        }
        
        # Discrete distribution
        else {
          if(right) prob = 1-pfunction(value-1,param1) 
          else      prob = pfunction(value,param1)
        }
      }
      
      dens = dfunction(value,param1)
    }
    
    # Chosen distribution has 3 parameters
    else if(is3P) {
          
      # prob. to quantile
      if(p) {
      
        prob = value
        
        # Continuous distribution
        if(!isDiscrete) {
          if(right) value = qfunction(1-prob,param1,param2,param3)
          else      value = qfunction(prob,param1,param2,param3)
        }
        
        # Discrete distribution
        else {
          if(right) value = qfunction(1-prob,param1,param2,param3)
          else      value = qfunction(prob,param1,param2,param3)
        }
      }
      
      # Quantile to prob.
      else {
      
        # Continuous distribution
        if(!isDiscrete) {
          if(right) prob = 1-pfunction(value,param1,param2,param3)
          else      prob = pfunction(value,param1,param2,param3)
        }
        
        # Discrete distribution
        else {
          if(right)  prob = 1-pfunction(value-1,param1,param2,param3)
          else       prob = pfunction(value,param1,param2,param3)
        }
      }
      
      dens = dfunction(value,param1,param2,param3)
    }
    
    # Result
    svalue(.$resultx) = paste(value)
    svalue(.$resultp) = paste(prob)
    
    # Construction du graphique
    xlab="X"
    title = paste(.$translate("Distribution"),distrib)
    ylab = expression(f(X==x))
    from = 0
    
    # Two-parameter distribution
    if(is2P) {
    
      # Continuous distribution
      if(!isDiscrete) { 
        from = ifelse(distrib==.$translate("Gaussian"),param1-4*param2,0)
        from = min(from,value,na.rm=TRUE)
        to = ifelse(distrib==.$translate("Gaussian"),param1+4*param2,max(rfunction(1000,param1,param2),na.rm=TRUE))
        to = max(to,value,na.rm=TRUE)
        curve(dfunction(x,param1,param2),n=1000,from=from,to=to,lwd=2,main=title,xlab=xlab,ylab=ylab)
        if(!right) {
          z = c(seq(from,value,len=1000),value,from)
          dz = c(dfunction(z[1:1000],param1,param2),0,0)
          polygon(z,dz,density=-1,col="red",lwd=2)
        }
        else {
          z = c(seq(value,to,len=1000),to,value)
          dz = c(dfunction(z[1:1000],param1,param2),0,0)
          polygon(z,dz,density=-1,col="red",lwd=2)
        }
      }
      
      # Discrete distribution
      else {
        from = 0
        to = ifelse(distrib==.$translate("Binomial"),param1,max(rfunction(1000,param1,param2)))
        z = 0:to
        plot(z,dfunction(z,param1,param2),type="h",lwd=2,main=title,xlab=xlab,ylab=ylab)
        
        if(!right) {
          for(i in 0:(value-1)) { lines(rbind(c(i,0),c(i,dfunction(i,param1,param2))),lwd=2,col="red") }
          lines(rbind(c(value,0),c(value,prob-pfunction(value-1,param1,param2))),lwd=2,col="red")
        }
        else {
          for(i in param1:(value+1)) {
          lines(rbind(c(i,0),c(i,dfunction(i,param1,param2))),lwd=2,col="red") }
          lines(rbind(c(value,0),c(value,prob-1+pfunction(value,param1,param2))),lwd=2,col="red")
        }
      }
    } 
    
    # One parameter distributions
    else if(is1P) {
    
      # Continuous distribution
      if(!isDiscrete) {
      
        from = ifelse(distrib==.$translate("Student"),min(rfunction(1000,param1)),0)
        from = min(from,value,na.rm=TRUE)
        to = max(rfunction(1000,param1),na.rm=TRUE)
        to = max(to,value,na.rm=TRUE)
        curve(dfunction(x,param1),n=1000,from=from,to=to,lwd=2,main=title,xlab=xlab,ylab=ylab)
        
        if(!right) {
          z = c(seq(from,value,len=1000),value,from)
          dz = c(dfunction(z[1:1000],param1),0,0)
          polygon(z,dz,density=-1,col="red",lwd=2)
        }
        else {
          z = c(seq(value,to,len=1000),to,value)
          dz = c(dfunction(z[1:1000],param1),0,0)
          polygon(z,dz,density=-1,col="red",lwd=2)
        }
      }
      
      # Discrete distribution
      else {
        from = 0
        to = max(rfunction(1000,param1),na.rm=TRUE)
        z = 0:to
        plot(z,dfunction(z,param1),type="h",lwd=2,main=title,xlab=xlab,ylab=ylab)
        
        if(!right) {
          for(i in 0:(value-1)) {
            lines(rbind(c(i,0),c(i,dfunction(i,param1))),lwd=2,col="red")
          }
          lines(rbind(c(value,0),c(value,prob-pfunction(value-1,param1))),lwd=2,col="red")
        }
        
        else {
          for(i in to:(value+1)) {
            lines(rbind(c(i,0),c(i,dfunction(i,param1))),lwd=2,col="red")
          }
          lines(rbind(c(value,0),c(value,prob-1+pfunction(value,param1))),lwd=2,col="red")
        }
      }
    }
    
    # Three-parameter distribution
    else if(is3P) {
    
      # Continuous distribution
      if(!isDiscrete) { 
        sample = rfunction(1000,param1,param2,param3)
        from = min(sample,value,na.rm=TRUE)
        to = max(sample,value,na.rm=TRUE)
        curve(dfunction(x,param1,param2,param3),n=1000,from=from,to=to,lwd=2,main=title,xlab=xlab,ylab=ylab)
        if(!right) {
          z = c(seq(from,value,len=1000),value,from)
          dz = c(dfunction(z[1:1000],param1,param2,param3),0,0)
          polygon(z,dz,density=-1,col="red",lwd=2)
        }
        else {
          z = c(seq(value,to,len=1000),to,value)
          dz = c(dfunction(z[1:1000],param1,param2,param3),0,0)
          polygon(z,dz,density=-1,col="red",lwd=2)
        }
      }
      
      # Discrete distribution
      else {
        from = 0
        to = max(rfunction(1000,param1,param2,param3),value,na.rm=TRUE)
        z = 0:to
        plot(z,dfunction(z,param1,param2,param3),type="h",lwd=2,main=title,xlab=xlab,ylab=ylab)
        
        if(!right) {
          for(i in 0:(value-1)) { lines(rbind(c(i,0),c(i,dfunction(i,param1,param2,param3))),lwd=2,col="red") }
          lines(rbind(c(value,0),c(value,prob-pfunction(value-1,param1,param2,param3))),lwd=2,col="red")
        }
        else {
          for(i in param1:(value+1)) {
          lines(rbind(c(i,0),c(i,dfunction(i,param1,param2,param3))),lwd=2,col="red") }
          lines(rbind(c(value,0),c(value,prob-1+pfunction(value,param1,param2,param3))),lwd=2,col="red")
        }
      }
    } 

  },
  initOptions = function(.,h,...) {

    distrib = svalue(.$distribution)
    param2 = eval(parse(text=svalue(.$param2)))
    
    is1P = distrib %in% .$translate(c("Student","Chi-2","Poisson"))
    is2P = distrib %in% .$translate(c("Gaussian","Uniform","Inverse Chi-2","Non central Chi-2","Fisher","Binomial","Gamma","Inverse Gamma","Beta"))
    is3P = distrib %in% .$translate(c("Generalized Student","Non central Student","Non central Fisher","Lambda prime"))
        
    # Warning: a droplist may be temporarily set to NULL in gWidgets when changed
    if(is.null(distrib)) return()
    
    svalue(.$paramLabel1) = .$paramNames[[distrib]][1]
    svalue(.$paramLabel2) = .$paramNames[[distrib]][2]
    svalue(.$paramLabel3) = .$paramNames[[distrib]][3]
    
    if(is1P) {
      svalue(.$param2) = ""
      enabled(.$param2) = FALSE
      svalue(.$param3) = ""
      enabled(.$param3) = FALSE
    }

    if(is2P) {
      svalue(.$param2) = ""
      enabled(.$param2) = TRUE
      svalue(.$param3) = ""
      enabled(.$param3) = FALSE
    }

    if(is3P) {
      svalue(.$param2) = ""
      enabled(.$param2) = TRUE
      svalue(.$param3) = ""
      enabled(.$param3) = TRUE
    }   
    
    svalue(.$resultx) = ""
    svalue(.$resultp) = ""    
  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-AtelieR")
  },
  #---------------------------------------------------------------------------------------
  #  SLOT                   INITIAL VALUE                          CONTENT
  #---------------------------------------------------------------------------------------
  availDists   = c(Gaussian="norm",Uniform="unif",Binomial="binom",Poisson="pois",                       #     Available distributions
                   Student="t","Generalized Student"="nst","Non central Student"="nct",
                  "Chi-2"="chisq","Inverse Chi-2"="scaledInvChi2","Non central Chi-2"="ncchisq",
                   Fisher="f","Non central Fisher"="ncf",
                   Gamma="gamma","Inverse Gamma"="invgamma",
                   Beta="beta","Lambda prime"="lambdaprime"),       
  paramNames   = list(Uniform               = c("Left boundary", "Right boundary", "       "),          #    Parameter names
                      Binomial              = c("Size         ", "Probability   ", "       "),
                      Gaussian              = c("Mean         ", "Standard dev. ", "       "),
                      Gamma                 = c("Shape        ", "Scale         ", "       "),
                     "Inverse Gamma"        = c("Shape        ", "Scale         ", "       "),
                     "Chi-2"                = c("df           ", "              ", "       "),
                     "Inverse Chi-2"        = c("df           ", "Scale         ", "       "),
                     "Non central Chi-2"    = c("df           ", "Non centrality", "       "),
                      Beta                  = c("alpha        ", "beta          ", "       "),
                      Poisson               = c("Mean         ", "              ", "       "),
                      Student               = c("df           ", "              ", "       "),
                     "Generalized Student"  = c("df           ", "Center        ", "Scale  "),
                     "Non central Student"  = c("df           ", "Non centrality", "Scale  "),
                     "Fisher"               = c("df1          ", "df2           ", "       "),
                     "Non central Fisher"   = c("df1          ", "df2           ", "Non centrality"),
                     "Lambda prime"         = c("df           ", "Non centrality", "Scale  ")),
  distribution = NULL,                                 #    Distribution choisie
  calcWhat     = NULL,                                 #    Type de calcul (de quantile à prob. ou l'inverse)
  side         = NULL,                                 #    Cumul à droite ou à gauche
  param1       = NULL,                                 #    Valeur du paramètre 1 de la loi 
  param2       = NULL,                                 #    Valeur du paramètre 2 de la loi (s'il y en a)   
  param3       = NULL,                                 #    Valeur du paramètre 3 de la loi (s'il y en a)   
  paramLabel1  = NULL,                                 #    Nom du paramètre 1 de la loi
  paramLabel2  = NULL,                                 #    Nom du paramètre 2 de la loi (s'il y en a)
  paramLabel3  = NULL,                                 #    Nom du paramètre 3 de la loi (s'il y en a)
  value        = NULL,                                 #    Valeur numérique de départ (quantile ou probabilité)
  resultx      = NULL,                                 #    Résultat quantile
  resultp      = NULL                                  #    Résultat probabilité
)


