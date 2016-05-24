# NAMESPACE issues: http://r.789695.n4.nabble.com/another-import-puzzle-td3408077.html
# getME() slots: "X" "Z" "Zt" "u" "Gp" "L" "Lambda" "Lambdat" "RX" "RZX" "beta" "theta" "REML" "n_rtrms" "is_REML"
#--------------------------------------------------------------------------------------------------
#
#                      R2STATS: A Graphical User Interface for GLM and GLMM in R
#                       Yvonnick Noel, University of Brittany, Rennes 2, France
#                                            2006-2010
#
#--------------------------------------------------------------------------------------------------
#                                        GLM PROTOTYPE
#--------------------------------------------------------------------------------------------------

r2sGLMM = proto(

  ###------------------------------------- Model constructor --------------------------------------------

  new = function(.,name,class,func,dvField,dv,ivField,iv,data,family,link,weights,constrFactor,subset) {
    .$proto(name=name,class=class,func=func,dvField=dvField,dv=dv,ivField=ivField,iv=iv,data=data,
            family=family,link=link,weights=weights,constrFactor=constrFactor,subset=subset)
  },

  ###---------------------------------- Model estimation method -----------------------------------------
  estimate = function(.) {
  
   .$estimated = FALSE

    family = .$getFamily()
    link   = .$getLink()

    # Is the dependent variable a matrix?
    dv = .$dvField
    if(length(.$dv) > 1) dv = paste("cbind(",.$dvField,")",sep="")
    
    # Estimate
    if( (family=="gaussian") && (link=="identity") ) {
      res = try(eval(parse(text=paste(".$Rmodel <- lmer","(",dv,"~",.$ivField,",data=",.$data,",subset=",.$subset,",weights=",.$weights,")",sep=""))))
    }
    else {
      res = try(eval(parse(text=paste(".$Rmodel <-",.$func,"(",dv,"~",.$ivField,",family=",family,"(link=",link,"),data=",.$data,",subset=",.$subset,",weights=",.$weights,")",sep=""))))
    }
    
    if(inherits(res,"try-error")) return(res)
        
    # Store initial factor names
    indicList = .$getIndicList()
    modelVars = .$getModelVars()
    
    # This allows to get design factors in model definition order
   .$designFactors  = modelVars[ modelVars %in% union(.$getFactorList(),indicList) ]
   
    # Exclude random factors from the IV list
   .$iv = .$getFIV()
  
    # Get explicit labels for factor or indic categories (e.g. Sex0/Sex1 instead of 0/1, -1/1 etc.)
    if(length(.$designFactors)) {
    
      # 'subset' seems to be the only way to get a data.frame when variable selection retains a single factor
      modelFactors = subset(.$getModelData(),select=.$designFactors)
      for(i in 1:ncol(modelFactors)) modelFactors[,i] = paste(.$designFactors[i],modelFactors[,i],sep="")
     .$groupLabels = factor(do.call("paste",modelFactors))
     
      # Set colors for plots
     .$groupFullColors   = r2stats$getColors(nlevels(.$groupLabels),"Set1")
      names(.$groupFullColors) = levels(.$groupLabels)
     .$groupPastelColors = r2stats$getColors(nlevels(.$groupLabels),"Pastel1")
      names(.$groupPastelColors) = levels(.$groupLabels)
    }
    
    else {
     .$groupFullColors   = r2stats$getColors(1,"Set1")
     .$groupPastelColors = r2stats$getColors(1,"Pastel1")
    }

    # If a contrast factor is defined, refit with this constraint
    if(.$constrFactor != .$translate("No factor")) {
      res = try(eval(parse(text=paste(".$Rmodel <-update(.$Rmodel,.~",.$getConstrainedFormula(),")"))))
      if(inherits(res,"try-error")) return(res)
    }

   .$estimated = TRUE
  },
  ###------------------------------------------------ Model summary and print method ------------------------------------------------
  Summary = function(.) {
  
    if(! .$estimated) {
      r2stats$setStatus(.$translate("Status: Ready."))
      return()
    }
    
    model.data   = .$getModelData()  
    familyIndex  = .$getFamilyAsIndex()
    family       =  names(r2stats$linkLists)[familyIndex]
    linkIndex    = .$getLinkAsIndex()
    link         =  r2stats$linkLists[[familyIndex]][linkIndex]
    residuals    = .$Residuals()

    # Subset command
    nobs = length(residuals)
    subset = rep(TRUE,nobs)
    if(.$subset != "NULL") subset = .$getSubset()

    r2stats$setStatus(.$translate("Status: Goodness of fit statistics..."))
    s = summary(.$Rmodel)
    coef.table = coef(s)

    if(nrow(coef.table)) {
      coef.table = round(coef.table,3)
      colnames(coef.table)[1:3]=.$translate(c("Estimate","Std. Error","t value"))
      if(.$hasIntercept()) rownames(coef.table)[1] = .$translate("Intercept")
    }

    res = data.frame(Statistic = .$translate(c("Loglikelihood","No. of parameters","AIC","BIC","Deviance")), 
                      Value=round(c(.$LogLik(),.$df(),.$Aic(),.$Bic(),.$Deviance()),3),row.names=1)
    names(res) = .$translate(names(res))
    
    # Group structure (if any) for additional tests
    if(length(.$designFactors))                        groups = .$groupLabels
    if( !(.$constrFactor %in% .$translate(c("No factor","Constant"))) ) groups = model.data[,.$constrFactor]
    
    r2stats$setStatus(.$translate("Status: Test of model assumptions..."))

    # Test normality and homogeneity of variances if family is gaussian and model not constant
    normtest = NULL
    vartest = NULL
    
    # Gaussian family
    if(familyIndex == 1) {

      # Shapiro-Wilk test for normality
      normtest = .$normalityTest(.$Residuals())
    
      # Levene test for homogeneity of variances
      if(exists("groups") && (.$constrFactor!=.$translate("Constant"))) .$varTest(.$Residuals(),groups)
    }

    r2stats$setStatus(.$translate("Status: Output of numerical results..."))
    
    # Display model name and specifications
    add(r2stats$results,paste(.$translate("Model"),.$name),font.attr=c(style="normal",weights="bold",size="large",col="blue"))
    modelProps = cbind(.$translate(c("Table","Formula","Constraint","Link","Distribution")),
                       c(.$data,.$getFormula(),.$getConstrFactor(),link,family))
    add(r2stats$results,capture.output(prmatrix(modelProps,rowlab=rep("",5),collab=rep("",2),quote=F)),font.attr=c(family="monospace",size="medium"))
    add(r2stats$results,"")

    # Display fit statistics
    add(r2stats$results,.$translate("Goodness of fit"), font.attr=c(style="normal",col="black",weights="bold"))
    add(r2stats$results,"")
    add(r2stats$results,capture.output(res),font.attr=c(family="monospace",size="medium"))
    add(r2stats$results,"")

    # Display test of normality
    if(!is.null(normtest)) { 
      add(r2stats$results,.$translate("Test of normality"),font.attr=c(style="normal",col="black",weights="bold"))
      add(r2stats$results,"")
      add(r2stats$results,capture.output(normtest),font.attr=c(family="monospace",size="medium"))
      add(r2stats$results,"")
    }
    
    # Display test of homogeneity of variance
    if(!is.null(vartest)) {
      add(r2stats$results,.$translate("Test of homogeneity of variances"), font.attr=c(style="normal",col="black",weights="bold"))
      add(r2stats$results,"")
      add(r2stats$results,capture.output(vartest),font.attr=c(family="monospace",size="medium"))
      add(r2stats$results,"")
    }
    
    # Display fitted values averaged by group, and group frequencies
    if(length(.$designFactors)) {
      add(r2stats$results,.$translate("Fitted values (by group)"),font.attr=c(style="normal",col="black",weights="bold"))
      add(r2stats$results,"")
      add(r2stats$results,capture.output(tapply(.$getPrediction(),model.data[,.$designFactors],mean)),font.attr=c(family="monospace",size="medium"))
      add(r2stats$results,"")
      add(r2stats$results,.$translate("Group sizes"),font.attr=c(style="normal",col="black",weights="bold"))
      add(r2stats$results,"")
      add(r2stats$results,capture.output(tapply(.$getPriorWeights(),model.data[,.$designFactors],sum)),font.attr=c(family="monospace",size="medium"))
      add(r2stats$results,"")
    }

    # Display estimates
    add(r2stats$results,.$translate("Fixed effects"),font.attr=c(style="normal",col="black",weights="bold"))
    add(r2stats$results,"")
    if(nrow(coef.table) > 0) {
      add(r2stats$results,capture.output(coef.table),font.attr=c(family="monospace",size="medium"))
      add(r2stats$results,"")
    }
    add(r2stats$results,.$translate("Random effects"),font.attr=c(style="normal",col="black",weights="bold"))
    add(r2stats$results,"")
    varcorr = VarCorr(.$Rmodel)
    if(length(varcorr) > 0) {

      # No more REmat matrix in this version of lme4! Reconstruct it.
      for(l in 1:length(varcorr)) {
        REmat = cbind(Groups="",Name=attr(varcorr[[l]],"dimnames")[[1]],"Std. dev."=round(attr(varcorr[[l]],"stddev"),3),Corr=round(attr(varcorr[[l]],"correlation"),3))
        REmat[1,1] = names(varcorr)[l] # Name of the random factor
        REmat = rbind(REmat,"")         # Add a line for the residual
        REmat[nrow(REmat),3] = round(attr(varcorr, "sc"), 3)
        colnames(REmat) = .$translate(colnames(REmat))
        REmat[nrow(REmat),1] = .$translate("Residual")
        add(r2stats$results,capture.output(print(as.data.frame(REmat),row.names=FALSE)),font.attr=c(family="monospace",size="medium"))
      }
      
      add(r2stats$results,"")
    }  

  },
  ###------------------------------------------------ Model plot methods ------------------------------------------------
  Plot = function(.,h) {
  
    # Model specs
    distr    = .$getFamily()
    liens    = .$getLink()
    varNames = .$getModelVars()
    is.vect  = .$getNumVarList()
    randFact = .$poolRandomFactors()

    # Get R2STATS plotting options
    plotType    = r2stats$getPlotTypeAsIndex()
    selGroup    = r2stats$getSelectedGroup()
    subsetVar   = r2stats$getSubset()
    legend.loc  = r2stats$getLegendLocation()
    legend.cols = r2stats$getLegendCols()
    xlim        = r2stats$getXLim()
    ylim        = r2stats$getYLim()
    addData     = svalue(r2stats$addData)
    addModel    = svalue(r2stats$addModel)
    addCondMeans= svalue(r2stats$addCondMeans)
    addGrid     = svalue(r2stats$addGrid)
    addRefLine  = svalue(r2stats$addRefLine)
    addNoise    = svalue(r2stats$addNoise)
    addSmooth   = svalue(r2stats$addSmooth)
    addRandCurves = svalue(r2stats$addRandCurves)

    # Get data
    current.data = .$getModelData()
    y <- yaxis <- .$getY()

    # Fitted values
    fit  = .$getPrediction()
    ffit = .$getFixedPrediction()

    # TODO: Data subset as defined on the graph tab

    # Set colors
    if(length(.$designFactors)) r2stats$setPlotParams(nlevels(.$groupLabels))
    else                        r2stats$setPlotParams(1)
    
    # Set the graphic device (except when we want to copy to clipboard or save to a file)
    if(h$action != "save") visible(r2stats$plotArea) = TRUE
    
    # 1 - Quantile-quantile plot
    if(plotType==4) {
    
      if(distr=="gaussian") {
      
        if(is.null(.$groupLabels)) {
        
          r2stats$currentPlot = qqmath(~.$Residuals(),main = paste(.$translate("Quantile-quantile plot"),.$name,sep=" - "),
                                        xlab = .$translate("Expected quantiles"),ylab=.$translate("Residuals"),
                                        panel = function(x, ...) {
                                           panel.qqmathline(x, ...)
                                           panel.qqmath(x, ...)
                                      })
        }
        
        else {

          r2stats$currentPlot = qqmath(~.$Residuals(),groups=.$groupLabels,
                                        key = list(space=legend.loc,text=list(levels(.$groupLabels)),col=.$groupFullColors,columns=legend.cols),
                                        panel = "panel.superpose",main = paste(.$translate("Quantile-quantile plot"),.$name,sep=" - "),
                                        xlab = .$translate("Expected quantiles"),ylab=.$translate("Residuals"),
                                        panel.groups = function(x, ...) {
                                           panel.qqmathline(x, ...)
                                           panel.qqmath(x, ...)
                                      })
        }
      }
      
      else if(distr=="Gamma") {
      
        shape = gamma.shape(.$Rmodel)$alpha
        fit = unique(fit)

        if(is.null(.$groupLabels)) {
        
          r2stats$currentPlot = qqmath(~y,distribution = function(x) qgamma(x,shape=shape,rate=shape/fit),
                                        main = paste(.$translate("Quantile-quantile plot"),.$name,sep=" - "),
                                        xlab = .$translate("Expected quantiles"),ylab=.$translate("Observed quantiles"),
                                        panel = function(x, ...) {
                                           panel.qqmathline(x, ...)
                                           panel.qqmath(x, ...)
                                      })
        }
        
        else {
        
          rate = shape/fit
          r2stats$currentPlot = qqmath(~y,groups=.$groupLabels,
                                        key = list(space=legend.loc,text=list(levels(.$groupLabels)),col=.$groupFullColors,columns=legend.cols),
                                        panel = "panel.superpose",main = paste(.$translate("Quantile-quantile plot"),.$name,sep=" - "),
                                        xlab = .$translate("Expected quantiles"),ylab=.$translate("Observed quantiles"),
                                        distribution = function(x) qgamma(x,shape,rate),
                                        prepanel = prepanel.qqmathline,
                                        panel.groups = function(x,group.number,...) {
                                           panel.qqmathline(x,rate=rate[group.number],...)
                                           panel.qqmath(x,...)
                                      })
        }
      }
      
      # No QQ plot for other distributions
      else r2stats$currentPlot = r2stats$emptyPlot()
    }

    # 2 - Regression plot
    else if(plotType==1) {
    
      # No numeric IV
      if(length(is.vect)==0) {
      
        xaxis = .$getFixedLinearTerm()
        
        # No factors neither: A constant model
        if(is.null(.$groupLabels)) {

          xlabel = .$translate("Estimated intercept")
          r2stats$currentPlot = xyplot(yaxis~xaxis,xlab=xlabel,ylab=.$dv[1],main=paste(.$translate("Regression"),.$name,sep=" - "),
                                       panel = function(x,y,...) {
                                         if(addData)  panel.xyplot(x,y,type="p",cex=.8)
                                         if(addModel) panel.xyplot(x,fit,type="p",pch="+",cex=1.3)
                                       })          
        }
        
        # Purely categorical model
        else {

          xlabel = .$translate("Fitted values")
          obsColors = r2stats$getColors(nlevels(randFact),"Set1")
          r2stats$currentPlot = xyplot(y~reorder(.$groupLabels,ffit),groups=randFact,xlab="",ylab=.$dv[1],main=paste(.$translate("Regression"),.$name,sep=" - "),
                                       auto.key=list(space=legend.loc,col=obsColors,points=F),
                                       panel = function(x,y,groups,subscripts,...) { 
                                         if(addGrid)       panel.grid(h=-1,v=-1)
                                         if(addData)       panel.superpose(x,y,groups,subscripts,type="p",jitter.x=addNoise,jitter.y=addNoise,col=obsColors)
                                         if(addRandCurves) panel.superpose(x,fit[subscripts],groups,subscripts,type="a",col=obsColors)
                                         if(addModel)      panel.superpose(x,ffit[subscripts],groups,subscripts,type="a",col.lines=.$groupFullColors)
                                       })
        }
      }
      
      # A single numeric variable: Take it as the x-axis
      else if(length(is.vect)==1) {
      
        xAxisType = svalue(r2stats$graphAxisX)
        if(xAxisType==.$translate("Default")) {
          xaxis = current.data[,is.vect]
          xlabel = is.vect
        }
        else {
          xaxis = .$getFullData()[,xAxisType]
          xlabel = xAxisType      
        }

        # No groups
        if(is.null(.$groupLabels)) {

          r2stats$currentPlot = xyplot(y~xaxis,groups=randFact,xlab=xlabel,ylab=.$dv[1],main=paste(.$translate("Regression"),.$name,sep=" - "),
                                       auto.key=FALSE,
                                       panel = function(x,y,groups,subscripts,...) { 
                                         if(addGrid)       panel.grid(h=-1,v=-1)
                                         if(addRefLine)    panel.abline(a=0,b=1,col="lightgrey",lty=2)
                                         if(addData)       panel.superpose(x,y,groups,subscripts,type="p",jitter.x=addNoise,jitter.y=addNoise)
                                         if(addCondMeans)  panel.average(x,y,fun=mean,horizontal=FALSE,groups,subscripts,col="darkgrey",lty=2)
                                         if(addRandCurves) panel.superpose(x,fit[subscripts],groups,subscripts,type="a")
                                         if(addModel)      panel.superpose(x,ffit[subscripts],groups,subscripts,type="a",col="black")
                                       })
          if(addSmooth)  r2stats$currentPlot = r2stats$currentPlot + layer(panel.xyplot(x,y,type="smooth",col.lines="darkgrey",lty=2))
        }
        
        # With a group structure
        else {
        
          gfcol = tapply(.$groupFullColors[.$groupLabels],randFact,unique)
          gpcol = tapply(.$groupPastelColors[.$groupLabels],randFact,unique)

          r2stats$currentPlot = xyplot(y~xaxis,groups=randFact,xlab=xlabel,ylab=.$dv[1],main=paste(.$translate("Regression"),.$name,sep=" - "),
                                       key = list(space=legend.loc,text=list(levels(.$groupLabels)),col=.$groupFullColors,columns=legend.cols),
                                       panel = function(x,y,groups,subscripts) { 
                                         if(addGrid)       panel.grid(h=-1,v=-1)
                                         if(addRefLine)    panel.abline(a=0,b=1,col="lightgrey",lty=2)
                                         if(addData)       panel.superpose(x,y,groups,subscripts,type="p",jitter.x=addNoise,jitter.y=addNoise,col=gfcol,fill=gpcol)
                                         if(addCondMeans)  panel.superpose(x,y,panel.groups=panel.average,fun=mean,horizontal=FALSE,groups=.$groupLabels,subscripts,type="a",col.lines=.$groupFullColors[.$groupLabels],lty=2)
                                         if(addRandCurves) panel.superpose(x,fit,groups,subscripts,type="a",col=gpcol)
                                         # if(addSmooth)     panel.superpose(x,y,panel.groups=panel.xyplot,type="smooth",groups=.$groupLabels,subscripts,col.lines=.$groupFullColors[.$groupLabels],lty=2)
                                         if(addSmooth)     panel.superpose(x,y,panel.groups=panel.loess,groups=.$groupLabels,subscripts,col.lines=.$groupFullColors[.$groupLabels],lty=2)
                                         if(addModel)      panel.superpose(x,ffit,panel.groups=panel.average,fun=mean,horizontal=FALSE,groups=.$groupLabels,subscripts,type="a",col.lines=.$groupFullColors[.$groupLabels])
                                       })
        }
      }
      
      # Several numeric predictors: Combine them
      else {
      
        # Prepare plot area
        xAxisType = svalue(r2stats$graphAxisX)
        if(xAxisType==.$translate("Default")) {
          xaxis = .$getNumLinearTerm()
          xlabel = .$translate("Predictor combination")        
        }
        
        else {
          xaxis = .$getFullData()[,xAxisType]
          xlabel = xAxisType
        }

        # No groups
        if(is.null(.$groupLabels)) {

          r2stats$currentPlot = xyplot(y~xaxis,groups=randFact,xlab=xlabel,ylab=.$dv[1],main=paste(.$translate("Regression"),.$name,sep=" - "),
                                       auto.key=FALSE,
                                       panel = function(x,y,groups,subscripts,...) { 
                                         if(addGrid)       panel.grid(h=-1,v=-1)
                                         if(addRefLine)    panel.abline(a=0,b=1,col="lightgrey",lty=2)
                                         if(addData)       panel.superpose(x,y,groups,subscripts,type="p",jitter.x=addNoise,jitter.y=addNoise)
                                         if(addCondMeans)  panel.average(x,y,fun=mean,horizontal=FALSE,groups,subscripts,type="p",col="lightgrey")
                                         if(addRandCurves) panel.superpose(x,fit[subscripts],groups,subscripts,type="a")
                                         if(addModel)      panel.superpose(x,ffit[subscripts],groups,subscripts,type="a",col="black")
                                       })
          if(addSmooth)  r2stats$currentPlot = r2stats$currentPlot + layer(panel.xyplot(x,y,type="smooth",col.lines="grey",lty=2))
        }
        
        # With a group structure
        else {
        
          gfcol = tapply(.$groupFullColors[.$groupLabels],randFact,unique)
          gpcol = tapply(.$groupPastelColors[.$groupLabels],randFact,unique)

          r2stats$currentPlot = xyplot(y~xaxis,groups=randFact,xlab=xlabel,ylab=.$dv[1],main=paste(.$translate("Regression"),.$name,sep=" - "),
                                       key = list(space=legend.loc,text=list(levels(.$groupLabels)),col=.$groupFullColors,columns=legend.cols),
                                       panel = function(x,y,groups,subscripts) { 
                                         if(addGrid)       panel.grid(h=-1,v=-1)
                                         if(addRefLine)    panel.abline(a=0,b=1,col="lightgrey",lty=2)
                                         if(addData)       panel.superpose(x,y,groups,subscripts,type="p",jitter.x=addNoise,jitter.y=addNoise,col=gfcol,fill=gpcol)
                                         if(addCondMeans)  panel.superpose(x,y,panel.groups=panel.average,fun=mean,horizontal=FALSE,groups=.$groupLabels,subscripts,type="p",col.lines=.$groupPastelColors[.$groupLabels])
                                         if(addRandCurves) panel.superpose(x,fit,groups,subscripts,type="a",col=gpcol)
                                         if(addSmooth)     panel.superpose(x,y,panel.groups=panel.xyplot,type="smooth",groups=.$groupLabels,subscripts,col.lines=.$groupPastelColors[.$groupLabels],lty=2)
                                         if(addModel)      panel.superpose(x,ffit,groups,subscripts,type="a",col=gfcol)
                                       })
        }
      }
    }

    # 3 - Prediction plot
    else if(plotType==3) {
    
      # No groups
      if(is.null(.$groupLabels)) {
        r2stats$currentPlot = xyplot(y~fit,xlab=.$translate("Fitted values"),ylab=.$dv[1],main=paste(.$translate("Prediction"),.$name,sep=" - "),
                                     panel = function(x,y,...) { 
                                       if(addGrid)    panel.grid(h=-1,v=-1)
                                       panel.abline(a=0,b=1,col="lightgrey",lty=2)
                                       panel.xyplot(x,y,type="p",...)
                                     })
      }
          
      # With a group structure
      else {
        r2stats$currentPlot = xyplot(y~fit,groups=.$groupLabels,xlab=.$translate("Fitted values"),ylab=.$dv[1],main=paste(.$translate("Prediction"),.$name,sep=" - "),
                                     key = list(space=legend.loc,text=list(levels(.$groupLabels)),col=.$groupFullColors,columns=legend.cols),
                                     panel = function(x,y,groups,subscripts,...) { 
                                       if(addGrid)    panel.grid(h=-1,v=-1)
                                       panel.abline(a=0,b=1,col="lightgrey",lty=2)
                                       panel.superpose(x,y,groups,subscripts,type="p",...)
                                     })
      }
    }

    # 4 - TODO: Response distribution plot
    else if(plotType==2) {
      r2stats$currentPlot = r2stats$emptyPlot()
    }
    
    # 5 - Residuals distribution
    else if(plotType==5) {
      # Only relevant for gaussian family
      if(distr=="gaussian") {
      
        resids = .$Residuals()
        h = hist(resids,plot=FALSE,freq=FALSE)
        r2stats$currentPlot = histogram(~resids,type="density",main=.$translate("Histogram"),xlab=.$translate("Residuals"),ylab=.$translate("Density"),
                                        panel = function(x,...) {
                                          panel.histogram(x,breaks=h$breaks,col=NULL,border=.$groupFullColors[1])
                                          panel.mathdensity(dnorm,args = list(mean=mean(x),sd=.$getDispersion()),col="black")
                                       })
      }
      else r2stats$currentPlot = r2stats$emptyPlot()
    }
    
    # 6 - Prediction and residuals
    else if(plotType==6) {
    
      # No groups
      if(is.null(.$groupLabels)) {
        r2stats$currentPlot = xyplot(.$Residuals(type="standard")~fit,xlab=.$translate("Fitted values"),ylab=.$translate("Standardized residuals"),main=paste(.$translate("Prediction and residuals"),.$name,sep=" - "),
                                     panel = function(x,y,...) { 
                                       if(addGrid)    panel.grid(h=-1,v=-1)
                                       panel.abline(a=1.96,b=0,col="lightgrey",lty=2)
                                       panel.abline(a=-1.96,b=0,col="lightgrey",lty=2)
                                       panel.xyplot(x,y,type="p",...)
                                     })
      }
          
      # With a group structure
      else {
        r2stats$currentPlot = xyplot(.$Residuals(type="standard")~fit,groups=.$groupLabels,xlab=.$translate("Fitted values"),ylab=.$translate("Standardized residuals"),
                                     main=paste(.$translate("Prediction and residuals"),.$name,sep=" - "),
                                     key = list(space=legend.loc,text=list(levels(.$groupLabels)),col=.$groupFullColors,columns=legend.cols),
                                     panel = function(x,y,groups,subscripts,...) { 
                                       if(addGrid)    panel.grid(h=-1,v=-1)
                                       panel.abline(a=1.96,b=0,col="lightgrey",lty=2)
                                       panel.abline(a=-1.96,b=0,col="lightgrey",lty=2)
                                       panel.superpose(x,y,groups,subscripts,type="p",...)
                                     })
      }
    }  

    # Is there a legend?
    if(legend.loc == "none") r2stats$currentPlot = update(r2stats$currentPlot,legend=NULL)
    
    # Are there restriction or dilation on the axis limits?
    if(!is.null(xlim)) r2stats$currentPlot = update(r2stats$currentPlot,xlim=xlim)
    if(!is.null(ylim)) r2stats$currentPlot = update(r2stats$currentPlot,ylim=ylim)
    
    # Want noise?
    if(addNoise) r2stats$currentPlot = update(r2stats$currentPlot,jitter.x=TRUE,jitter.y=TRUE)    
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 ACCESSOR METHODS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Get model name
  getName = function(.) {
    return(.$name)
  },
  ### Get model class ("glm","mer",...)
  getClass = function(.) {
    return(.$class)
  },
  ### Get dataframe name
  getDataName = function(.) {
    return(.$data)
  },
  ### Get dependent variable name
  getDV = function(.) {
    return(.$dv)
  },
  ### Get dependent variable
  getY = function(.) {
    return(.$Rmodel@resp$y)
  },
  ### Get independent variable names
  getIV = function(.) {
    return(.$iv)
  },
  ### Get *fixed* independent variable names
  getFIV = function(.) {
    return(.$iv[ !(.$iv %in% .$getRandomFactorList())])
  },
  getModelVars = function(.) {
    all.vars(terms(.$Rmodel))
  },
  ### Get the names of fixed factor predictors in the model
  getFactorList = function(.) {
    is.fact = sapply(.$Rmodel@frame,is.factor)
    is.fact = names(is.fact)[is.fact]
    is.fact[ !(is.fact %in% .$getRandomFactorList()) && !(is.fact %in% .$dv) ]
  },
  # Get the names of the *random* factors
  getRandomFactorList = function(.) {
    names(getME(.$Rmodel,"flist"))
  },
  ### Get the names of numeric predictors
  getNumVarList = function(.) {
    is.num = sapply(.$Rmodel@frame,is.numeric)
    is.num = names(is.num)[is.num]
    is.num[ !(is.num %in% .$dv) ]
  },
  ### Get the names of indicator (numeric binary) variables
  getIndicList = function(.) {
    if(!length(.$iv)) return(NULL)
    .$iv[sapply(.$getModelData()[,.$iv], function(var) length(unique(var))==2)]
  },
  ### Get model formula
  getFormula = function(.) {
    paste(.$dvField,.$ivField,sep="~")
  },
  ### Get constrained formula
  getConstrainedFormula = function(.) {
  
    Formula = .$ivField
    factor.list = union(.$getFactorList(),.$getIndicList())
    constrF = .$getConstrFactor()
    
    if(constrF == .$translate("Constant")) return("1")
    
    # This replaces all factor occurences by the constraint factor
    for(v in factor.list) Formula = gsub(v,constrF,Formula)
    Formula
  },
  ### Get distribution family in R
  getFamily = function(.) {
    names(.$RLinkList)[.$family]
  },
  getFamilyAsIndex = function(.) {
    return(.$family)
  },
  ### Get link function in R
  getLink = function(.) {
    family = .$getFamily()
    .$RLinkList[[family]][.$link]
  },
  getLinkAsIndex = function(.) {
    return(.$link)
  },
  ### Get the list of available plots for this model
  getAvailablePlots = function(.) {
    return(.$translate(.$graphTypes))
  },
  ### Get the name of the constraint factor
  getConstrFactor = function(.) {
    return(.$constrFactor)
  },
  ### Get the data subtable used in model estimation
  getModelData = function(.) {

    formula = paste(c(.$dv,.$iv),collapse="+")
    
    # If a (non constant) constraint factor is defined, it must be added to the model frame
    if( !(.$constrFactor %in% .$translate(c("No factor","Constant")))) formula = paste(c(formula,.$getConstrFactor()),collapse="+")

    # With model frames, subsetting one variable as a data.frame is possible
    eval(parse(text=paste("model.frame(~",formula,",data=",.$data,",subset=",.$subset,")",sep="")),envir=.GlobalEnv)
  },
  ### Get the full data table
  getFullData = function(.) {
  
    # Original data set is not stored in the model object so we go fetch it in the global env.
    eval(parse(text=.$Rmodel@call[[3]]),env=.GlobalEnv)
  },
  ### Get individual fitted values
  getPrediction = function(.) {
    fitted(.$Rmodel)
  },
  ### Get fixed (group) predicted values
  getFixedPrediction = function(.) {
    .$invLink(.$getFixedLinearTerm())
  },
  ### Get the eta linear term (including random effects)
  getLinearTerm = function(.) {
    .$Rmodel@eta
  },
  ### Get the eta linear term (including random effects)
  getFixedLinearTerm = function(.) {
    getME(.$Rmodel,"X") %*% fixef(.$Rmodel)
  },
  ### Get the linear term on numeric variables only
  getNumLinearTerm = function(.) {
  
    # Factors are binary recoded in design matrix
    X = getME(.$Rmodel,"X")
    whichNum = apply(X,2,function(x) length(unique(x))>2)
    
    # Numeric prediction
    numCoefs = cbind(coef(summary(.$Rmodel))[whichNum,1])
    cbind(X)[,whichNum] %*% numCoefs
  },
  ### Get the prior weights
  getPriorWeights = function(.) {

    w = weights(.$Rmodel)
    if(length(w)) return(w)
    else          return(rep(1,length(w)))
  },
  ### Boolean: Has the model an intercept
  hasIntercept = function(.) {
    attr(terms(.$Rmodel),"intercept")
  },
  ### Evaluate the subset command
  getSubset = function(.) {
    eval(parse(text=paste("with(",.$data,",",.$subset,")",sep="")),envir=.GlobalEnv)
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                             VARIABLE EDITION TOOLS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Get response function (inverse of link function
  invLink = function(.,x) {
  
    link = .$getLink()
    
         if(link == "identity") return(x)
    else if(link == "log")      return(exp(x))
    else if(link == "inverse")  return(1/x)
    else if(link == "logit")    return(plogis(x))
    else if(link == "probit")   return(pnorm(x))
    else if(link == "cauchit")  return(0.5 + atan(x)/pi)
    else if(link == "cloglog")  return(1-exp(-exp(x)))
    else if(link == "sqrt")     return(x**2)
    else if(link == "1/mu^2")   return(sqrt(1/x))
  },
  ### Make a single factor from all random variables
  poolRandomFactors = function(.) {
     factor(do.call("paste",.$Rmodel@flist))
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                 EDITION TOOLS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Remove cbind() from a string
  strip.cbind = function(l) {
    l = sapply(l,function(x) sub("cbind","",x))
    l = sapply(l,function(x) sub("\\(","",x))
    sapply(l,function(x) sub("\\)","",x))
  },
  ### Gettext utility for translating messages
  translate = function(.,...) {
    gettext(..., domain="R-R2STATS")
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                              ADDITIONAL TESTS AND STATISTICS
  #
  #------------------------------------------------------------------------------------------------------------------------
  ### Homoskedasticity test
  varTest = function(.,y,group) {
  
    nl = nlevels(group)
    if(nl < 2) return()
    
    # Two groups
    if(nl == 2) {
      vt = var.test(y~group)
      table = c(vt$statistic,round(vt$parameter),p=vt$p.value)
      names(table) = .$translate(c("F","Num.Df","Den.Df","Prob."))
    }
    
    # More than two groups
    else {
      meds = tapply(y, group, median, na.rm = TRUE)
      resp = abs(y - meds[group])
      table = anova(lm(resp ~ group))[, c(1, 4, 5)]
      colnames(table) = .$translate(c("Df","Levene's F","Pr(>F)"))
      rownames(table) = .$translate(c("Groups","Error"))
    }
    
    table
  },
  ### Test of normality
  normalityTest = function(.,y) {
  
    # Shapiro.test accepts 5000 obs. max.
    if(length(y) > 5000) normtest = shapiro.test(sample(y,5000))
    else                 normtest = shapiro.test(y)

    normtest = data.frame(W=round(normtest$statistic,3),Prob=round(normtest$p.value,3))
    row.names(normtest)="Shapiro-Wilk"
    
    normtest
  },
  ### AIC
  Aic = function(.) { 
    AIC(.$LogLik())
  },
  ### BIC function
  Bic = function(.) {
    BIC(.$LogLik())  
  },  
  ### Loglikelihood
  LogLik = function(.,REML=FALSE) { 
    logLik(.$Rmodel,REML=REML)
  },
  ### Degrees of freedom
  df = function(.) { 
    attr(logLik(.$Rmodel),"df")
  },
  ### Deviance
  Deviance = function(.) { 
    deviance(.$Rmodel)
  },
  ### Percentage of deviance explained
  devExplained = function(.) { 
  
    # Null deviance
    #eval(parse(text="M0 <- glm(",.$dvField,"~1,family=",family,"(link=",link,"),data=",.$data,",subset=",.$subset,",weights=",.$weights,")",sep=""))
    #nullDeviance = deviance()

    #1 - (deviance(.$Rmodel)/nullDeviance)
  },
  ### Get model residuals
  Residuals = function(.,type="raw") {
  
    if(type=="raw")           residuals(.$Rmodel)
    else if(type=="standard") residuals(.$Rmodel)/.$getDispersion()
  },
  getRandomVar = function(.) {
    summary(.$Rmodel)@sigma  
  },
  getDispersion = function(.) {
    attr(VarCorr(.$Rmodel),"sc")
  },
  #------------------------------------------------------------------------------------------------------------------------
  #
  #                                                    SLOTS
  #
  #------------------------------------------------------------------------------------------------------------------------
  # SLOT                    INITIAL VALUE                                               DESCRIPTION
  #------------------------------------------------------------------------------------------------------------------------
  Rmodel        = NULL,                                                              # R model object
  name          = "",                                                                # Model name
  class         = "mer",                                                             # Model type: GLM, GLMM, MANOVA, VGLM...
  func          = "glmer",                                                           # Estimation function
  data          = NULL,                                                              # Data set name
  dvField       = "",                                                                # Content of the DV field as a string
  dv            = NULL,                                                              # Vector of dependent variables names
  ivField       = "",                                                                # Content of the IV field as a string
  iv            = NULL,                                                              # Vector of independent variables names (including rand. vars)
  family        = 1,                                                                 # Index of the family function (see list below)
  link          = 1,                                                                 # Index of the link function (see list below)
  RLinkList     = list(gaussian = c("identity","log","inverse"),                     # List of vectors of link function by family
                      binomial = c("logit","probit","cauchit","log","cloglog"),
                      poisson  = c("log","sqrt","identity"),
                      Gamma    = c("inverse","log","identity"),
                      inverse.gaussian = c("inverse","log","identity","1/mu^2")),      # List of link functions
  weights       = "NULL",                                                            # Name of the weighting variable ("NULL" if none)
  constrFactor  = NULL,                                                              # Name of the constraint factor
  designFactors = NULL,                                                              # Names of the original factors (before constraints)
  subset        = "NULL",                                                            # Character vector containing the subsetting command ("NULL" if none)
  estimated     = FALSE,                                                             # Flag: Model is estimated (TRUE) or not (FALSE)
  groupFullColors   = 1,                                                             # Colors to be used in plots
  groupPastelColors = 1,                                                             # Colors to be used in plots
  groupLabels   = NULL,                                                              # Labels for the groups
  graphTypes    = c("Regression plot","Response distribution",                       # Vector of possible plots for this type of model
                   "Fitted and observed values","Quantile-quantile plot",
                   "Residuals distribution","Fitted values and residuals")
)
