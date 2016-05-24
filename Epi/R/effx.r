## Program to calculate effects
## Michael Hills
## Improved by Bendix Carstensen and Martyn Plummer
## Post Tartu 2007 version June 2007
## Addition allowing a TRUE/FALSE as binary outcome
## and possibility or relative risk for binary outcomes and rate
## difference for failure outcomes, BxC, Ocitober 2012

effx<-function(response,
type="metric",
fup=NULL,
exposure,
strata=NULL,
control=NULL,
weights=NULL,
eff=NULL,
alpha=0.05,
base=1,
digits=3,
data=NULL)
{
  ## stores the variable names for response, etc.

    rname<-deparse(substitute(response))
    ename<-deparse(substitute(exposure))
    if (!missing(strata)) {
      sname<-deparse(substitute(strata))
    }

  ## The control argument is more complex, as it may be a name or
  ## list of names

  if(!missing(control)) {
    control.arg <- substitute(control)
    if (length(control.arg) > 1) {
      control.names <- sapply(control.arg, deparse)[-1]
    }
    else {
      control.names <- deparse(control.arg)
    }
  }

  ## Match the type argument
  type <- match.arg(type, c("metric", "failure", "count", "binary"))

  ## Check for missing arguments
  if (missing(response))
    stop("Must specify the response","\n")
  if (missing(exposure))
    stop("Must specify the exposure","\n")
  if (type == "failure" && missing(fup)) {
    stop("Must specify a follow-up variable when type is failure")
  }

  ## performs a few other checks
  if(rname==ename)stop("Same variable specified as response and exposure")
  if (!missing(strata)) {
    if(rname==sname)stop("Same variable specified as response and strata")
    if(sname==ename)stop("Same variable specified as strata and exposure")
  }


  ## If data argument is supplied, evaluate the arguments in that
  ## data frame.

  if (!missing(data)) {
    exposure <- eval(substitute(exposure), data)
    response <- eval(substitute(response), data)
    if (!missing(strata))
      strata <- eval(substitute(strata), data)
    if (!missing(control))
      control <- eval(substitute(control), data)
    if (!missing(fup))
      fup <- eval(substitute(fup), data)
    if (!missing(weights)) {
      weights <- eval(substitute(weights), data)
    }
  }

  ## Now check validity of evaluated arguments

  if(is.logical(response)) response <- as.numeric(response)
  if(!is.numeric(response))
    stop("Response must be numeric, not a factor")

  if (!missing(weights) && type != "binary") {
    stop("weights only allowed for a binary response")
  }

  if (!missing(strata) && !is.factor(strata))
    stop("Stratifying variable must be a factor")

  if(type=="binary") {
    if( is.null(eff) ) eff<-"OR"
    if( !(eff %in% c("OR","RR") ) ) stop( "Only RR and OR allowed for binary response" )
    response <- as.numeric(response)
    tmp<-(response==0 | response==1)
    if(all(tmp,na.rm=TRUE)==FALSE)
      stop("Binary response must be logical or coded 0,1 or NA")
  }

  # If a count is given we are actually using the fup
  if(type=="count") fup<-is.na(response)*1

  if(type=="failure") {
    if( is.null(eff) ) eff<-"RR"
    response <- as.numeric(response)
    tmp<-(response==0 | response==1)
    if(all(tmp,na.rm=TRUE)==FALSE)
      stop("Failure response must be logical or coded 0,1 or NA")
  }

  ## If exposure is an ordered factor, drops the order.

  if(class(exposure)[1]=="ordered") {
    exposure<-factor(exposure, ordered=F)
  }

  ## Fix up the control argument as a named list
  if (!missing(control)) {
    if (is.list(control)) {
      names(control) <- control.names
    }
    else {
      control <- list(control)
      names(control) <- control.names
    }
  }

  ## prints out some information about variables
  cat("---------------------------------------------------------------------------","\n")
  cat("response      : ", rname, "\n")
  cat("type          : ", type, "\n")
  cat("exposure      : ", ename, "\n")
  if(!missing(control))cat("control vars  : ",names(control),"\n")
  if(!missing(strata)) {
    cat("stratified by : ",sname,"\n")
  }
  cat("\n")
  if(is.factor(exposure)) {
    cat(ename,"is a factor with levels: ")
    cat(paste(levels(exposure),collapse=" / "),"\n")
    exposure <- Relevel( exposure, base )
    cat( "baseline is ", levels( exposure )[1] ,"\n")
  }
  else {
    cat(ename,"is numeric","\n")
  }
  if(!missing(strata)) {
    cat(sname,"is a factor with levels: ")
    cat(paste(levels(strata),collapse="/"),"\n")
  }
  if(type=="metric")cat("effects are measured as differences in means","\n")
  if(type=="binary")
    {
    if( eff=="OR" | is.null(eff))cat("effects are measured as odds ratios","\n")
    if( eff=="RR" )cat("effects are measured as relative risk","\n")
    }
  if(type=="failure")
    {
    if( eff=="RR" | is.null(eff))cat("effects are measured as rate ratios","\n")
    if( eff=="RD" )cat("effects are measured as rate differences","\n")
    }
  cat("---------------------------------------------------------------------------","\n")
  cat("\n")

  ## translates type of response and choice of eff into family
  if ( type=="metric")                  family<-gaussian
  if ( type=="binary")                  family<-binomial(link=logit)
  if ( type=="failure" | type=="count") family<-poisson(link=log)
  if ( !is.null(eff) )
     {
  if (type=="binary"  & eff=="RR")   family<-binomial(link=log)
  if (type=="failure" & eff=="RD")   family<-poisson(link=identity)
     }
  ## gets number of levels for exposure if a factor

  if(is.factor(exposure)) {
    nlevE<-length(levels(exposure))
  }

  ## labels the output

  if(is.factor(exposure)) {
    cat("effect of",ename,"on",rname,"\n")
  }
  else {
    cat("effect of an increase of 1 unit in",ename,"on",rname,"\n")
  }
  if(!missing(control)) {
    cat("controlled for",names(control),"\n\n")
  }
  if(!missing(strata)) {
    cat("stratified by",sname,"\n\n")
  }


  ## no stratifying variable

  if(missing(strata)) {
    if(type=="metric") {
      if(missing(control)) {
        m<-glm(response~exposure,family=family)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response~1,family=family,subset=!is.na(exposure))
      }
      else  {
        m<-glm(response~.+exposure,family=family,
               subset=!is.na(exposure),data=control)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response~.,family=family,
                subset=!is.na(exposure),data=control)
      }
      res<-ci.lin(m,subset=c("Intercept","exposure"),alpha=alpha)
      res<-res[,c(1,5,6)]
    }
    if(type=="binary") {
      if(missing(control)) {
        m<-glm(response~exposure,family=family,weights=weights)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response~1,family=family,subset=!is.na(exposure),weights=weights)
      }
      else  {
        m<-glm(response~.+exposure,family=family,
               subset=!is.na(exposure),data=control,weights=weights)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response~.,family=family,
                subset=!is.na(exposure),data=control,weights=weights)
      }
      res<-ci.lin(m,subset=c("Intercept","exposure"),Exp=TRUE,alpha=alpha)
      res<-res[,c(5,6,7)]
    }
    if (type=="failure" | type=="count") {
      if (missing(control)) {
        m<-glm(response/fup~exposure,weights=fup,family=family)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response/fup~1,weights=fup,family=family,
                subset=!is.na(exposure))
      }
      else  {
        m<-glm(response/fup~.+exposure,weights=fup,family=family, data=control)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response/fup~.,weights=fup,family=family,
                subset=!is.na(exposure),data=control)
      }
      res<-ci.exp(m,subset=c("Intercept","exposure"),alpha=alpha,Exp=(eff=="RR"))
    }
    res<-signif(res,digits)
    colnames(res)[1]<-c("Effect")
    if(is.factor(exposure)) {
      ln <- levels(exposure)
      rownames(res)[2:nlevE]<-paste(ln[2:nlevE],"vs",ln[1])
    }
    aov <- anova(mm,m,test="Chisq")
    print( res[-1,] )
    cat("\nTest for no effects of exposure on",
        aov[2,3],"df:",
        "p-value=",format.pval(aov[2,5],digits=3),"\n")
    invisible(list(res,paste("Test for no effects of exposure on",
                             aov[2,3],"df:","p-value=",format.pval(aov[2,5],digits=3))))

  }
  ## stratifying variable
  if(!missing(strata)) {
    sn <- levels(strata)
    nlevS<-length(levels(strata))
    if(type=="metric") {
      if(missing(control)) {
        m<-glm(response~strata/exposure,family=family)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response~strata+exposure,family=family)
      }
      else {
        m <-glm(response~strata/exposure + .,family=family,
                data=control)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm <-glm(response~strata+exposure + .,family=family,
                 data=control)
      }
      res<-ci.lin(m,subset=c("strata"),alpha=alpha)[c(-1:-(nlevS-1)),c(1,5,6)]
    }
    if(type=="binary") {
      if(missing(control)) {
        m<-glm(response~strata/exposure,family=family,weights=weights)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response~strata+exposure,family=family,weights=weights)
      }
      else {
        m <-glm(response~strata/exposure + .,family=family,
                data=control,weights=weights)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm <-glm(response~strata+exposure + .,family=family,
                 data=control,weights=weights)
      }
      res<-ci.exp(m,subset=c("strata"),alpha=alpha)[c(-1:-(nlevS-1)),]
    }
    if (type=="failure" | type=="count") {
      if(missing(control)) {
        m<-glm(response/fup~strata/exposure,weights=fup,family=family)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response~strata+exposure,weights=fup,family=family)
      }
      else {
        m <-glm(response/fup~.+strata/exposure,weights=fup,family=family,data=control)
        cat("number of observations ",m$df.null+1,"\n\n")
        mm<-glm(response/fup~.+strata+exposure,weights=fup,family=family,data=control)
      }
      res<-ci.exp(m,subset=c("strata"),alpha=alpha)[c(-1:-(nlevS-1)),]
    }

    res<-signif(res,digits)
    colnames(res)[1]<-c("Effect")
    if(is.factor(exposure)) {
      ln<-levels(exposure)
      newrownames<-NULL
      for(i in c(1:(nlevE-1))) {
        newrownames<-c(newrownames,
                       paste("strata",sn[1:nlevS],"level",ln[i+1],"vs",ln[1]))
      }
    }
    else {
      newrownames<-paste("strata",sn[1:nlevS])
    }
    rownames(res)<-newrownames
    aov<-anova(mm,m,test="Chisq")
    print( res )
    cat("\nTest for effect modification on",
        aov[2,3],"df:","p-value=",format.pval(aov[2,5],digits=3),"\n")
    invisible(list(res,paste("Test for effect modification on",
                             aov[2,3],"df:","p-value=",format.pval(aov[2,5],digits=3))))

  }
}
