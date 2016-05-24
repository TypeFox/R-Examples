## Program to calculate effects for matched case-control studies
## Michael Hills
## Improved by BxC and MP
## Post Tartu 2007 version June 2007


effx.match<-function(response,
exposure,
match,
strata=NULL,
control=NULL,
base=1,
digits=3,
alpha=0.05,
data=NULL)
{

  ## stores the variable names for response, etc.

    rname<-deparse(substitute(response))
    ename<-deparse(substitute(exposure))
    if (!missing(strata))sname<-deparse(substitute(strata))

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
  ## If data argument is supplied, evaluate the arguments in that
  ## data frame.

  if (!missing(data)) {
    exposure <- eval(substitute(exposure), data)
    response <- eval(substitute(response), data)
    match <- eval(substitute(match),data)
    if (!missing(strata)) {
      strata <- eval(substitute(strata), data)
    }
    if (!missing(control))
      control <- eval(substitute(control), data)
  }

  ## performs a few other checks

  if(rname==ename)stop("Same variable specified as response and exposure")
  if (!missing(strata)) {
    if(rname==sname)stop("Same variable specified as response and strata")
    if(sname==ename)stop("Same variable specified as strata and exposure")
  }

  if(!is.numeric(response))stop("Response must be numeric, not a factor")
  
  if(!missing(strata)&!is.factor(strata))stop("Stratifying
    variable must be a factor")

  tmp<-(response==0 | response==1)
  if(all(tmp,na.rm=TRUE)==FALSE)
  stop("Binary response must be coded 0,1 or NA")

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
  cat("exposure      : ", ename, "\n")
  if(!missing(control))cat("control vars  : ",names(control),"\n")
  if(!missing(strata)) cat("stratified by : ",sname,"\n")
  cat("\n")
  if(is.factor(exposure)) {
    cat(ename,"is a factor with levels: ")
    cat(paste(levels(exposure),collapse=" / "),"\n")
    cat( "baseline is ", levels( exposure )[base] ,"\n")   
    exposure <- Relevel( exposure, base )
  }
  else {
    cat(ename,"is numeric","\n")
  }
  if(!missing(strata)) {
    cat(sname,"is a factor with levels: ")
    cat(paste(levels(strata),collapse="/"),"\n")
  }

  cat("effects are measured as odds ratios","\n")
  cat("---------------------------------------------------------------------------","\n")
  cat("\n")

  ## gets number of levels for exposure if a factor

  if(is.factor(exposure)) {
    nlevE<-length(levels(exposure))
  }
  else {
    nlevE<-1
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
          if(missing(control)) {            
            m<-clogit(response~exposure+strata(match))
            cat("number of observations ",m$n,"\n\n")
          }
          else  {
              m<-clogit(response~.+exposure+strata(match),
                     subset=!is.na(exposure),data=control)
              cat("number of observations ",m$n,"\n\n")
              mm<-clogit(response~.+strata(match),
                      subset=!is.na(exposure),data=control)
          }
          res<-ci.lin(m,subset=c("exposure"),Exp=TRUE,alpha=alpha)[,c(5,6,7)]
          res<-signif(res,digits)
          if(nlevE<3) {
              names(res)[1]<-c("Effect")
          }
          else {
              colnames(res)[1]<-c("Effect")
              if(is.factor(exposure)) {
                ln <- levels(exposure)
                rownames(res)[1:nlevE-1]<-paste(ln[2:nlevE],"vs",ln[1])
              }
          }
          print(res)
          if(missing(control)) {
              chisq<-round(summary(m)$logtest[1],2)
              df<-round(summary(m)$logtest[2])
              p<-round(summary(m)$logtest[3],3)
              cat("\n")
              cat("Test for no effects of exposure:  ","\n")
              cat("chisq=",chisq, " df=",df, " p-value=",format.pval(p,digits=3),"\n")
              invisible(list(res,paste("Test for no effects of exposure on",
                 df,"df:","p=",format.pval(p,digits=3))))

          }
          else  {
              aov <- anova(mm,m,test="Chisq")
              cat("\nTest for no effects of exposure on",
              aov[2,3],"df:",
              "p-value=",format.pval(aov[2,5],digits=3),"\n")
              invisible(list(res,paste("Test for no effects of exposure on",
                 aov[2,3],"df:","p=",format.pval(aov[2,5],digits=3))))
          }
  }      
   ## stratifying variable

  if(!missing(strata)) {
      sn <- levels(strata)
      nlevS<-length(levels(strata))
          if(missing(control)) {
            m<-clogit(response~strata/exposure+strata(match))
            cat("number of observations ",m$n,"\n\n")
            mm<-clogit(response~strata+exposure+strata(match))
          }
          else {
            m <-clogit(response~strata/exposure + . +strata(match),
            data=control)
            cat("number of observations ",m$n,"\n\n")
            mm <-clogit(response~strata+exposure + . +strata(match),
            data=control)
          }
          res<-ci.lin(m,Exp=TRUE,alpha=alpha,subset="strata")[c(-1:-(nlevS-1)),c(5,6,7)]
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
