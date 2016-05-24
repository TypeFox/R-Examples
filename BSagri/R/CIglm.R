

CIGLM<-function(x, conf.level=0.95, method=c("Raw", "Adj", "Bonf"))
{

method<-match.arg(method)

switch(method,
"Raw"={
CI<-confint(x,level=conf.level, calpha=univariate_calpha())
adjnam<-"Unadjusted"
},

"Adj"={
CI<-confint(x,level=conf.level, calpha=adjusted_calpha())
adjnam<-"Adjusted"
},

"Bonf"={
ncomp<-nrow(x$linfct)
CI<-confint(x,level=1-(1-conf.level)/ncomp, calpha=univariate_calpha())
adjnam<-"Bonferroni-adjusted"
})

return(CI)
}



minus2slash<-function(x)
{
sl <- strsplit(x,"-")
out <- unlist(lapply(X=sl,FUN=function(x){paste(x,collapse="/")}))
return(out)
}


UnlogCI<-function(x){UseMethod("UnlogCI")}


# # # #

UnlogCI.glht<-function(x)
{
if(all(class(x)!="glht"))
 {stop("x must be of class glht")}


# # # fit by glm

if(class(x$model)[1] %in% c("glm", "geeglm"))
{

fam<-x$model$family$family

if( !fam %in% c("poisson", "quasipoisson", "binomial", "quasibinomial") )
 {stop("Family of the fitted glm should be one of poisson, quasipoisson, binomial or quasibinomial")}

link<-x$model$family$link

if(!link %in% c("log", "logit"))
 {stop("The link function of the fitted glm should be log or logit")}

if( fam %in% c("poisson", "quasipoisson") )
 {
 CI<-x$confint
 CIout<-exp(CI)
 cnam <- dimnames(CI)[[1]]
 dimnames(CIout)[[1]]<-minus2slash(cnam)
 para<-"Confidence intervals for the ratios of abundance"

 out<-x
 out$conf.int<-CIout
 out$parameter <- para

 }

if(fam %in% c("binomial", "quasibinomial"))
 {
 CI<-x$confint
 CIout<-exp(CI)
 cnam <- dimnames(CI)[[1]]
 dimnames(CIout)[[1]]<-minus2slash(cnam)
 sf<-dimnames(x$model$model[[1]])[[2]]

 para<-"Confidence intervals for odds ratios"

 if(!is.null(sf[1]))
  {
  ORdef<-paste(c("with the odds defined as p(",sf[1],")/(1-p(",sf[1],"))"), collapse="" )
  }
 else
  {ORdef<-""}

out<-x
out$conf.int<-CIout
out$parameter <- paste(para, ORdef, sep="\n")

 }
}

# # # fit by glm.nb 

if(class(x$model)[1]=="negbin")
{

link<-x$model$family$link

if(!link=="log")
 {stop("The link function of the fitted glm should be log")}

 CI<-x$confint
 CIout<-exp(CI)
 cnam <- dimnames(CI)[[1]]
 dimnames(CIout)[[1]]<-minus2slash(cnam)
 para<-"Confidence intervals for the ratios of abundance"

 out<-x
 out$conf.int<-CIout
 out$parameter <- para

}

# # # fit by lme4 (older versions)

if(class(x$model)[1]=="glmer")
{

fam<-x$model@family$family

if( !fam %in% c("poisson", "quasipoisson", "binomial", "quasibinomial") )
 {stop("Family of the fitted glmer should be one of poisson, quasipoisson, binomial or quasibinomial")}

link<-x$model@family$link

if(!link %in% c("log", "logit"))
 {stop("The link function of the fitted glm should be log or logit")}

if( fam %in% c("poisson", "quasipoisson") )
 {
 CI<-x$confint
 CIout<-exp(CI)
 cnam <- dimnames(CI)[[1]]
 dimnames(CIout)[[1]]<-minus2slash(cnam)
 para<-"Confidence intervals for the ratios of abundance"

 out<-x
 out$conf.int<-CIout
 out$parameter <- para

 }

if(fam %in% c("binomial", "quasibinomial"))
 {
 CI<-x$confint
 CIout<-exp(CI)
 cnam <- dimnames(CI)[[1]]
 dimnames(CIout)[[1]]<-minus2slash(cnam)
 sf<-dimnames(x@frame[[1]])[[2]]


 para<-"Confidence intervals for odds ratios"

 if(!is.null(sf[1]))
  {
  ORdef<-paste(c("with the odds defined as p(",sf[1],")/(1-p(",sf[1],"))"), collapse="" )
  }
 else
  {ORdef<-""}

out<-x
out$conf.int<-CIout
out$parameter <- paste(para, ORdef, sep="\n")
 }


}


# # # fit by lmer, later versions

if(class(x$model)[1]=="mer")
{

callchar<-paste(deparse(x$model@call), collapse="")
bin <- grep("binomial", callchar, value=TRUE)
pois <- grep("poisson", callchar, value=TRUE)

if( length(bin)==0 & length(pois)==0 )
 {stop("Family of the fitted model should be one of poisson, quasipoisson, binomial or quasibinomial")}

link<-grep("link", callchar, value=TRUE)

if(length(link)!=0)
 {warning("Note: The application of this function most probably makes sense only when applied to models with log or logit link!")}

 CI<-x$confint
 CIout<-exp(CI)
 cnam <- dimnames(CI)[[1]]
 dimnames(CIout)[[1]]<-minus2slash(cnam)
 para<-NULL

 out<-x
 out$conf.int<-CIout
 out$parameter <- para

}


# # # fit by gamlss

if(class(x$model)[1]=="gamlss")
{

fam<-x$model$family[1]
paras<-x$model$parameters
link<-x$model$mu.link

if( fam=="NBI" & link=="log")
 {
  CI<-x$confint
  CIout<-exp(CI)
  cnam <- dimnames(CI)[[1]]
  dimnames(CIout)[[1]]<-minus2slash(cnam)
  para<-"Confidence intervals for the ratios of abundance"

  out<-x
  out$conf.int<-CIout
  out$parameter <- para
 }

else{
 if( fam=="BB" & link=="logit")
  {
 CI<-x$confint
 CIout<-exp(CI)
 cnam <- dimnames(CI)[[1]]
 dimnames(CIout)[[1]]<-minus2slash(cnam)
 sf<-dimnames(x$model$model[[1]])[[2]]

 para<-"Confidence intervals for odds ratios"

out<-x
out$conf.int<-CIout
out$parameter <- paste(para, "\n")


  }
 else{
  stop("No methods except for families 'NBI' and 'BB' implemented for objects of class 'gamlss'")  
 } 
}

}


class(out)<-c("UnlogCI", class(out))
return(out)
}




print.UnlogCI<-function(x,...)
{
args<-list(...)

if(is.null(args$digits))
 {args$digits<-3}

call<-x$model$call

print(call)

out<-as.matrix(x$conf.int)
attr(out,"error")<-NULL
attr(out,"calpha")<-NULL
attr(out,"conf.level")<-NULL

args$x<-as.table(out)

cat(x$parameter,"\n")

do.call("print", args)

cat("Estimated quantile = ", round(attr(x$confint, which="calpha"),4), "\n")

}


