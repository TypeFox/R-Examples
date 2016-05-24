###{{{ twinlm

##' Fits a classical twin model for quantitative traits.
##'
##' @title Classic twin model for quantitative traits
##' @return   Returns an object of class \code{twinlm}.
##' @author Klaus K. Holst
##' @seealso \code{\link{bptwin}}, \code{\link{twinlm.time}}, \code{\link{twinlm.strata}}, \code{\link{twinsim}}
##' @aliases twinlm twinlm.strata
##' @export
##' @examples
##' ## Simulate data
##' set.seed(1)
##' d <- twinsim(1000,b1=c(1,-1),b2=c(),acde=c(1,1,0,1))
##' ## E(y|z1,z2) = z1 - z2. var(A) = var(C) = var(E) = 1
##' 
##' ## E.g to fit the data to an ACE-model without any confounders we simply write
##' ace <- twinlm(y ~ 1, data=d, DZ="DZ", zyg="zyg", id="id")
##' ace
##' ## An AE-model could be fitted as
##' ae <- twinlm(y ~ 1, data=d, DZ="DZ", zyg="zyg", id="id", type="ae")
##' ## LRT:
##' lava::compare(ae,ace)
##' ## AIC
##' AIC(ae)-AIC(ace)
##' ## To adjust for the covariates we simply alter the formula statement
##' ace2 <- twinlm(y ~ x1+x2, data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
##' ## Summary/GOF
##' summary(ace2)
##' \donttest{ ## Reduce Ex.Timings
##' ## An interaction could be analyzed as:
##' ace3 <- twinlm(y ~ x1+x2 + x1:I(x2<0), data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
##' ace3
##' ## Categorical variables are also supported##' 
##' d2 <- transform(d,x2cat=cut(x2,3,labels=c("Low","Med","High")))
##' ace4 <- twinlm(y ~ x1+x2cat, data=d2, DZ="DZ", zyg="zyg", id="id", type="ace")
##' }
##' @keywords models
##' @keywords regression
##' @param formula Formula specifying effects of covariates on the response
##' @param data \code{data.frame} with one observation pr row. In
##'     addition a column with the zygosity (DZ or MZ given as a factor) of
##'     each individual much be
##'     specified as well as a twin id variable giving a unique pair of
##'     numbers/factors to each twin pair
##' @param id The name of the column in the dataset containing the twin-id variable.
##' @param zyg The name of the column in the dataset containing the
##'     zygosity variable
##' @param DZ Character defining the level in the zyg variable
##'     corresponding to the dyzogitic twins. If this argument is missing,
##'     the reference level (i.e. the first level) will be interpreted as
##'     the dyzogitic twins
##' @param group Optional. Variable name defining group for interaction analysis (e.g., gender)
##' @param group.equal If TRUE marginals of groups are asummed to be the same
##' @param strata Strata variable name
##' @param weight Weight matrix if needed by the chosen estimator. For use
##'     with Inverse Probability Weights
##' @param type Character defining the type of analysis to be
##'     performed. Should be a subset of "aced" (additive genetic factors, common
##'     environmental factors, unique environmental factors, dominant
##'     genetic factors).
##' @param twinnum The name of the column in the dataset numbering the
##'     twins (1,2). If it does not exist in \code{data} it will
##'     automatically be created.
##' @param binary If \code{TRUE} a liability model is fitted. Note that if the right-hand-side of the formula is a factor, character vector, og logical variable, then the liability model is automatically chosen (wrapper of the \code{bptwin} function).
##' @param keep Vector of variables from \code{data} that are not
##'     specified in \code{formula}, to be added to data.frame of the SEM
##' @param estimator Choice of estimator/model
##' @param constrain Development argument
##' @param control Control argument parsed on to the optimization routine
##' @param messages Control amount of messages shown 
##' @param ... Additional arguments parsed on to lower-level functions
twinlm <- function(formula, data, id, zyg, DZ, group=NULL,
                   group.equal=FALSE, strata=NULL, weight=NULL, type=c("ace"),
                   twinnum="twinnum",
                   binary=FALSE,keep=weight,estimator="gaussian",
                   constrain=TRUE,control=list(),messages=1,...)
{
    
  cl <- match.call(expand.dots=TRUE)
  opt <- options(na.action="na.pass")
  mf <- model.frame(formula,data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  ## formula <- update(formula, ~ . + 1)
  yvar <- getoutcome(formula)
  if (missing(zyg)) stop("Zygosity variable not specified")
  if (!(zyg%in%colnames(data))) stop("'zyg' not found in data")
  if (!(id%in%colnames(data))) stop("'id' not found in data")
  if (missing(id)) stop("Twin-pair variable not specified")

  if (binary | is.factor(data[,yvar]) | is.character(data[,yvar]) | is.logical(data[,yvar])) {
    args <- as.list(cl)
    args[[1]] <- NULL
    return(do.call("bptwin",args,envir=parent.frame()))
  }
  
  formulaId <- unlist(Specials(formula,"cluster"))
  formulaStrata <- unlist(Specials(formula,"strata"))
  formulaSt <- paste("~.-cluster(",formulaId,")-strata(",paste(formulaStrata,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    id <- formulaId
    cl$id <- id
  }
  if (!is.null(formulaStrata)) strata <- formulaStrata
  cl$formula <- formula
 
  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
      fit <- lapply(seq(length(dd)),function(i) {
        if (messages>0) message("Strata '",names(dd)[i],"'")
        cl$data <- dd[[i]]
        eval(cl)
      })
      res <- list(model=fit)
      res$strata <- names(res$model) <- names(dd)
      class(res) <- c("twinlm.strata","twinlm")
      res$coef <- unlist(lapply(res$model,coef))
      res$vcov <- blockdiag(lapply(res$model,vcov))
      res$N <- length(dd)
      res$idx <- seq(length(coef(res$model[[1]])))
      rownames(res$vcov) <- colnames(res$vcov) <- names(res$coef)
      return(res)
    }
  }
  
  type <- tolower(type)
  ## if ("u" %in% type) type <- c("ue")
  
  varnames <- all.vars(formula)
  latentnames <- c("a1","a2","c1","c2","d1","d2","e1","e2")
  if (any(latentnames%in%varnames))
    stop(paste(paste(latentnames,collapse=",")," reserved for names of latent variables.",sep=""))
  
  M <- model.matrix(formula,mf)
  options(opt)  
  covars <- colnames(M)
  hasIntercept <- FALSE
  if (attr(terms(formula),"intercept")==1) {
      hasIntercept <- TRUE
      covars <- covars[-1]
  }
  if(length(covars)<1) covars <- NULL
  
  zygstat <- data[,zyg]
  if(!is.factor(zygstat)) {
    zygstat <- as.factor(zygstat)
  }
  zyglev <- levels(zygstat)
  if (length(zyglev)>2) stop("More than two zygosity levels found. For opposite sex (OS) analysis use the 'group' argument (and regroup OS group as DZ, e.g. DZ=c('OS','DZ'))")

  if (tolower(type)=="cor") type <- "u"
  if (!is.null(group) && type%in%c("u","flex","sat")) stop("Only polygenic models are allowed with 'group' ('type' subset of 'acde'). See also the 'strata' argument.")      
  ## To wide format:
  num <- NULL; if (twinnum%in%colnames(data)) num <- twinnum
  if (!is.null(group)) data[,group] <- as.factor(data[,group])
  data <- cbind(data[,c(yvar,keep,num,zyg,id,group)],M)
  ddd <- fast.reshape(data,id=c(id),varying=c(yvar,keep,covars,group),keep=zyg,num=num,sep=".",labelnum=TRUE)
  groups <- paste(group,".",1:2,sep="")
  outcomes <- paste(yvar,".",1:2,sep="")  

  if (missing(DZ)) {
    warning("Using first level, `",zyglev[1],"', in status variable as indicator for 'dizygotic'", sep="")
    DZ <- zyglev[1]    
  }
  MZ <- setdiff(zyglev,DZ)
  wide1 <- ddd[which(ddd[,zyg]==MZ),,drop=FALSE]
  wide2 <- ddd[which(ddd[,zyg]%in%DZ),,drop=FALSE]

  mm <- nn <- c()
  dd <- list()
  levgrp <- NULL
  if (!is.null(group)) {
      levgrp <- levels(data[,group])
      for (i1 in levgrp) {
          for (i2 in levgrp) {
              idxMZ <- which(wide1[,groups[1]]==i1 & wide1[,groups[2]]==i2)
              dMZ <- wide1[idxMZ,,drop=FALSE]
              idxDZ <- which(wide2[,groups[1]]==i1 & wide2[,groups[2]]==i2)
              dDZ <- wide2[idxDZ,,drop=FALSE]
              m0 <- twinsem1(outcomes,c(i1,i2),
                             levels=levgrp,covars=covars,type=type,
                             data=list(dMZ,dDZ),constrain=constrain,
                             equal.marg=group.equal,intercept=hasIntercept)$model
              if (length(idxMZ)>0) {
                  nn <- c(nn,paste("MZ:",i1,sep=""))
                  dd <- c(dd,list(dMZ))
                  mm <- c(mm,list(m0[[1]]))
              }
              if (length(idxDZ)>0) {
                  nn <- c(nn,paste("DZ:",i1," ",i2,sep=""))
                  dd <- c(dd,list(dDZ))
                  mm <- c(mm,list(m0[[2]]))
              }
          }
      }; names(mm) <- nn; names(dd) <- nn      
  } else {
      mm <- twinsem1(outcomes,NULL,
                     levels=NULL,covars=covars,type=type,
                     data=list(wide1,wide2),constrain=constrain,
                     intercept=hasIntercept)$model
      dd <- list(MZ=wide1,DZ=wide2)
  }
  
  newkeep <- unlist(sapply(keep, function(x) paste(x, 1:2, sep = ".")))
  if (is.null(estimator)) return(multigroup(mm, dd, missing=TRUE,fix=FALSE,keep=newkeep,type=2))
  optim <- list(method="nlminb2",refit=FALSE,gamma=1,start=rep(0.1,length(coef(mm[[1]]))*length(mm)))

    optim$start <- twinlmStart(formula,na.omit(mf),type,hasIntercept,surv=inherits(data[,yvar],"Surv"),model=mm, group=levgrp, group.equal=group.equal)
  if (length(control)>0) {
    optim[names(control)] <- control
  }

  if (inherits(data[,yvar],"Surv")) {
      if (!requireNamespace("lava.tobit",quietly=TRUE)) stop("lava.tobit required")
      if (is.null(optim$method))
          optim$method <- "nlminb1"
      suppressWarnings(e <- estimate(mm,dd,control=optim,...))
  } else {
      suppressWarnings(e <- estimate(mm,dd,weight=weight,estimator=estimator,fix=FALSE,control=optim,...))
  }

  if (!is.null(optim$refit) && optim$refit) {
    optim$method <- "NR"
    optim$start <- pars(e)
    if (inherits(data[,yvar],"Surv")) {
        suppressWarnings(e <- estimate(mm,dd,estimator=estimator,fix=FALSE,control=optim,...))
    } else {
        suppressWarnings(e <- estimate(mm,dd,weight=weight,estimator=estimator,fix=FALSE,control=optim,...))
    }
  }

  e$vcov <- Inverse(information(e,type="hessian"))
  
  counts <- function(dd) {
    dd0 <- apply(dd,2,function(x) !is.na(x))
    pairs <- sum(dd0[,1]*dd0[,2])
    singletons <- sum((!dd0[,1])*dd0[,2] + (!dd0[,2])*dd0[,1])
    return(c(pairs,singletons))
  }
  counts <- lapply(dd, function(x) counts(x[,outcomes]))
  ## mz  <- counts(object$data.mz[,object$outcomes])
  ## dz  <- counts(object$data.dz[,object$outcomes])
  if (!e$model$missing) {
      zygtab <- c("MZ-pairs"=counts[[1]][1],"DZ-pairs"=counts[[2]][1])
  } else {
      zygtab <- c(paste(counts[[1]],collapse="/"),paste(counts[[2]],collapse="/"))
      names(zygtab) <- c("MZ-pairs/singletons","DZ-pairs/singletons")
  }

  res <- list(coefficients=e$opt$estimate, vcov=e$vcov,
              estimate=e, model=mm, call=cl, data=data, zyg=zyg,
              id=id, twinnum=twinnum, type=type,  group=group,
              constrain=constrain, outcomes=outcomes, zygtab=zygtab,
              nam=nn, groups=levgrp, group.equal=group.equal,
              counts=counts)
  class(res) <- "twinlm"
  return(res)
}

###}}} twinlm

###{{{ twinsem1 (create lava model)

##outcomes <- c("y1","y2"); groups <- c("M","F"); covars <- NULL; type <- "ace"
##twinsem1(c("y1","y2"),c("M","F"))
twinsem1 <- function(outcomes,groups=NULL,levels=NULL,covars=NULL,type="ace",data,constrain=TRUE,equal.marg=FALSE,intercept=TRUE,...) {
    isA <- length(grep("a",type))>0 & type!="sat"
    isC <- length(grep("c",type))>0
    isD <- length(grep("d",type))>0
    isE <- length(grep("e",type))>0 | type=="sat" | type=="u"
    lambdas <- c("lambda[a]","lambda[c]","lambda[d]","lambda[e]")
    varidx <- which(c(isA,isC,isD,isE))
    vMZ1 <- c("a1","c1","d1","e1")
    vMZ2 <- c("a1","c1","d1","e2")
    vDZ2 <- c("a2","c1","d2","e2")
    rhoA <- rhoD <- zA <- zD <- NULL        
    if (is.list(outcomes)) {
        if (!is.null(groups) & is.null(levels)) stop("missing levels")
        if (is.null(groups)) groups <- c("","")
        grp <- paste(sort(groups),collapse=" ")
        sameGroup <- groups[1]==groups[2]
        model1 <- outcomes[[1]]
        model2 <- outcomes[[2]]
        if (!is.null(levels)) {
            pars <- c()
            for (i in seq(length(levels)-1)) for (j in seq(i+1,length(levels)))
                pars <- c(pars,paste(sort(levels)[c(i,j)],collapse=" "))
            if (isA) {
                parameter(model1) <- paste("z(A):",pars,sep="")
                parameter(model2) <- paste("z(A):",pars,sep="")
            }
            if (isD) {
                parameter(model1) <- paste("z(D):",pars,sep="")
                parameter(model2) <- paste("z(D):",pars,sep="")
            }
        }
        outcomes <- endogenous(model1)
        if (!(type%in%c("u","flex","sat"))) {
            if (equal.marg) {
                regression(model1,to=outcomes[1],from=vMZ1[varidx]) <-
                    lambdas[varidx]
                regression(model1,to=outcomes[2],from=vMZ2[varidx]) <-
                    lambdas[varidx]            
                regression(model2,to=outcomes[1],from=vMZ1[varidx]) <-
                    lambdas[varidx]
                regression(model2,to=outcomes[2],from=vDZ2[varidx]) <-
                    lambdas[varidx]
            } else {
                regression(model1,to=outcomes[1],from=vMZ1[varidx]) <-
                    paste(lambdas,groups[1],sep="")[varidx]
                regression(model1,to=outcomes[2],from=vMZ2[varidx]) <-
                    paste(lambdas,groups[2],sep="")[varidx]            
                regression(model2,to=outcomes[1],from=vMZ1[varidx]) <-
                    paste(lambdas,groups[1],sep="")[varidx]
                regression(model2,to=outcomes[2],from=vDZ2[varidx]) <-
                    paste(lambdas,groups[2],sep="")[varidx]
            }
        }
        if (sameGroup) {
            if (isA) covariance(model2,a1~a2)  <- 0.5
            if (isD) covariance(model2,d1~d2)  <- 0.25
        } else {            
            if (isA) {
                rhoA <- paste("Kinship[A]:",grp,sep="")
                zA <- paste("z(A):",grp,sep="")
                covariance(model2, a1~a2) <- rhoA
                constrain(model2, rhoA,zA) <- function(x) tanh(x)
            }
            if (isD) {                
                rhoD <- paste("Kinship[D]:",grp,sep="")
                zD <- paste("z(D):",grp,sep="")
                covariance(model2, d1~d2) <- rhoD
                constrain(model2, rhoD,zD) <- function(x) tanh(x)
            }
        }
        return(list(model=list(MZ=model1,DZ=model2),
                    zA=zA, rhoA=rhoA, zD=zD, rhoD=rhoD))
    }

    ### Build model from scratch....
    model1<- lvm()    
    if (!(type%in%c("u","flex","sat"))) {
        model1 <- regression(model1,to=outcomes[1],from=vMZ1[varidx])
        model1 <- regression(model1,to=outcomes[2],from=vMZ2[varidx])
    } else {
        addvar(model1) <- outcomes
        covariance(model1) <- outcomes        
    }    
    latent(model1) <- union(vMZ1[varidx],vMZ2[varidx])
    intercept(model1,latent(model1)) <- 0

        if (!is.null(covars))
            for (i in 1:length(covars)) {
            regression(model1, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta[",i,"]",sep="")
            regression(model1, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta[",i,"]",sep="")
        }
        covariance(model1,outcomes) <- 0
        covariance(model1, latent(model1)) <- 1

        if (!type%in%c("sat","flex")) {    
            intercept(model1,outcomes) <- "mu"
        }
        if (type%in%c("u","flex","sat")) {
            kill(model1) <- ~e1+e2
            covariance(model1,outcomes) <- "v1"
        }
        model2 <- model1
        if (!(type%in%c("u","flex","sat"))) {
            model2 <- cancel(model2,c(outcomes[2],vMZ2[varidx]))
            model2 <- regression(model2,to=outcomes[2],from=vDZ2[varidx])
            latent(model2) <- vDZ2[varidx]
            intercept(model2, latent(model2)) <- 0
            covariance(model2, latent(model2)) <- 1
        }
    if (type=="flex") {
        varMZ <- paste("var(MZ)",groups,sep="")
        varDZ <- paste("var(DZ)",groups,sep="")
        intercept(model1,outcomes) <- "mu1"
        intercept(model2,outcomes) <- "mu2"
        covariance(model1,outcomes) <- varMZ
        covariance(model2,outcomes) <- varDZ
    }
    if (type=="sat") {
        varMZ <- c(paste("var(MZ)1",groups[1],""),
                   paste("var(MZ)2",groups[2],""))
        varDZ <- c(paste("var(DZ)1",groups[1],""),
                   paste("var(DZ)2",groups[2],""))        
        covariance(model1,outcomes) <- varMZ
        covariance(model2,outcomes) <- varDZ
    }
    if (type%in%c("u","flex","sat")) {
        if (constrain) {
            if (type=="sat") {
                model1 <- covariance(model1,outcomes,constrain=TRUE,rname="atanh(rhoMZ)",cname="covMZ",lname="log(var(MZ)).1",l2name="log(var(MZ)).2")
                model2 <- covariance(model2,outcomes,constrain=TRUE,rname="atanh(rhoDZ)",cname="covDZ",lname="log(var(DZ)).1",l2name="log(var(DZ)).2")
            } else {
                if (type=="flex") {
                    model1 <- covariance(model1,outcomes,constrain=TRUE,rname="atanh(rhoMZ)",cname="covMZ",lname="log(var(MZ))")
                    model2 <- covariance(model2,outcomes,constrain=TRUE,rname="atanh(rhoDZ)",cname="covDZ",lname="log(var(DZ))")
                }  else {
                    model1 <- covariance(model1,outcomes,constrain=TRUE,rname="atanh(rhoMZ)",cname="covMZ",lname="log(var)")
                    model2 <- covariance(model2,outcomes,constrain=TRUE,rname="atanh(rhoDZ)",cname="covDZ",lname="log(var)")
                }        
            }     
        } else {
            covariance(model1,outcomes[1],outcomes[2]) <- "covMZ"
            covariance(model2,outcomes[1],outcomes[2]) <- "covDZ"
        }
    }

    if (!is.null(covars) & type%in%c("flex","sat")) {
        sta <- ""
        if (type=="sat") sta <- "b"
        for (i in 1:length(covars)) {
            regression(model1, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta1[",i,"]",sep="")         
            regression(model1, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta1",sta,"[",i,"]",sep="")
            regression(model2, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta2[",i,"]",sep="")
            regression(model2, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta2",sta,"[",i,"]",sep="")
        }
    }

    if (!intercept) {
        intercept(model1,outcomes) <- 0
        intercept(model2,outcomes) <- 0
    }
   
    ## Full rank covariate/design matrix?
    if (!missing(data))
    for (i in covars) {
        myvars <- paste(i,c(1,2),sep=".")
        dif <- data[[1]][,myvars[1]]-data[[1]][,myvars[2]]   
        mykeep <- myvars
        if (all(na.omit(dif)==00)) {
            mykeep <- mykeep[-2]
        }   
        trash <- setdiff(myvars,mykeep)
        if (length(mykeep)==1) {
            regression(model1, to=outcomes[2], from=mykeep) <- lava::regfix(model1)$label[trash,outcomes[2]]
            kill(model1) <- trash
        }
        dif <- data[[2]][,myvars[1]]-data[[2]][,myvars[2]]   
        mykeep <- myvars
        if (all(na.omit(dif)==00)) {
            mykeep <- mykeep[-2]
        }  
        trash <- setdiff(myvars,mykeep)
        if (length(mykeep)==1) {
            regression(model2, to=outcomes[2], from=mykeep) <- lava::regfix(model2)$label[trash,outcomes[2]]
            kill(model2) <- trash
        }
    }
    twinsem1(list(MZ=model1,DZ=model2),groups=groups,type=type,levels=levels,constrain=constrain,equal.marg=equal.marg)
}

###}}} twinsem1

###{{{ twinlmStart (starting values)

twinlmStart <- function(formula,mf,type,hasIntercept,surv=FALSE,model,group=NULL,group.equal,...)  {
    if (surv) {
        l <- survival::survreg(formula,mf,dist="gaussian")
        beta <- coef(l)
        sigma <- l$scale
    } else {
        l <- lm(formula,mf)
        beta <- coef(l)
        sigma <- summary(l)$sigma
    }    
    start <- rep(sigma/sqrt(nchar(type)),nchar(type))
    if (hasIntercept) {
        start <- c(beta[1],start)
        start <- c(start,beta[-1])
    } else start <- c(start,beta)
    if (type=="sat") {
        varp <- c(rep(log(sigma^2),2),0.5)
        start <- c()
        if (hasIntercept) {
            start <- rep(beta[1],4)
            beta <- beta[-1]
        }
        start <- c(start,rep(c(rep(beta,2),varp),2))
    }
    if (type=="flex") {
        varp <- c(log(sigma^2),0.5)
        start <- c()
        if (hasIntercept) {
            start <- c(beta[1],beta[1])
            beta <- beta[-1]
        }
        start <- c(start,rep(c(beta,varp),2))
    }
    if (type=="u") {
        start <- c(log(sigma^2),0.5,0.5)
        start <- c(beta,start)
    }
    names(start) <- NULL

    if (!is.null(group)) {
        cc <- coef(model[[1]])        
        iA <- grep(c("z(A):"),cc,fixed=TRUE)
        iD <- grep(c("z(D):"),cc,fixed=TRUE)
        nplus <- length(iA) + length(iD) + max(length(iA),length(iD))
        start <- c(start,rep(.2,nplus))
        ii <- sort(c(iA,iD))
        start <- c(start[seq(ii[1]-1)],rep(0.3,length(ii)),start[seq(ii[1]+1,length(start))])
    }
    
    return(start)
}

###}}} twinlmStart
