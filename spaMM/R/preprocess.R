checkRandLink <- function(rand.family) {
  lcrandfamfam <- tolower(rand.family$family) ## tolower once and for all
  oklink <- F
  ## cases where g(u)=th(u)
  if (lcrandfamfam=="gaussian" && rand.family$link=="identity") oklink <- T          
  if (lcrandfamfam=="gamma" && rand.family$link=="log") oklink <- T          
  if (lcrandfamfam=="inverse.gamma" && rand.family$link=="-1/mu") oklink <- T
  if (lcrandfamfam=="beta" && rand.family$link=="logit") oklink <- T
  ## cases where g(u)!=th(u)
  if (lcrandfamfam=="inverse.gamma" && rand.family$link=="log") oklink <- T 
  if (lcrandfamfam=="gamma" && rand.family$link=="identity") oklink <- T ## gamma(identity)
  if ( ! oklink) {
    allowed <- switch(lcrandfamfam,
                      gaussian= "is 'identity'",
                      gamma= "is 'log'", ## gamma(identity) not yet working
                      beta= "is 'logit'",
                      "inverse.gamma" = "are '-1/mu' and 'log'"
    )
    mess <- paste("(!) rand.family/link combination not handled;\nallowed link(s) for rand.family '",rand.family$family,"' ",allowed,sep="")
    stop(mess)
  }
  lcrandfamfam
}

checkRandLinkS <- function(rand.families) {
  rand.families <- lapply(rand.families, function(rf) {
    if (is.character(rf)) {
      rf <-switch(tolower(rf),
                  gaussian = gaussian(),
                  gamma = Gamma(link=log), ## NOT the default link
                  beta = Beta(), 
                  "inverse.gamma" = inverse.Gamma(),
                  stop("rand.family argument not valid"))
    }
    return(rf) 
  })
  lcrandfamfam <- unlist(lapply(rand.families,checkRandLink)) ## a _vector_ of lcrandfamfam := tolower(rand.family$family)
  return(list(rand.families=rand.families,lcrandfamfam=lcrandfamfam))
}

checkRespFam <- function(family) {
  ## four lines from glm()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (family$family=="Gamma") family <- spaMM_Gamma(family$link)
  family
}

getProcessed <- function(object,element,from=NULL) {
  if ( ! is.null(attr(object,"multiple"))) {
    if (is.null(from)) from <- seq_len(length(object))
    if (length(from)>1L) {
      resu <- lapply(from,function(id) {eval(parse(text=paste("object[[",id,"]]$",element,sep="")))})
      names(resu) <- names(object[from])
      return(resu) ## a list, e.g a data list
    } else {
      object <- object[[from]]
    } 
  } 
  return(eval(parse(text=paste("object$",element,sep=""))))
}


setProcessed <- function(object,element,value=1) {
  if ( ! is.null(attr(object,"multiple"))) {
    #    for (it in seq_len(length((object))) {eval(parse(text=paste("object[[",it,"]]$",element," <- ",value,sep="")))} 
    ## fails bc nam loses its enclosing "s :
    #for (nam in names(object)) {eval(parse(text=paste("object[[",nam,"]]$",element," <- ",value,sep="")))} 
    for (nam in names(object)) {eval(parse(text=paste("object[[\"",nam,"\"]]$",element," <- ",value,sep="")))} 
  } else eval(parse(text=paste("object$",element," <- ",value,sep="")))
  return(object)
}

generateInitPhi <- function(formula,data,family,weights=NULL) {
  ## look for replicates to estimate phi
  form <- subbarsMM(asStandardFormula(formula)) ## removes all random effect decorum (but retains its variables)
  # lhs
  mf <- model.frame(form,data=data)
  Y <- model.response(mf) ## evaluated rhs (e.g. log(y) rather than y...)
  # rhs
  rhs_terms <- delete.response(terms(form))
  mf <- model.frame(rhs_terms, data)
  # builds a model which indexes responses with identical predictor [ULI(mf)] 
  # and retains data that are replicates for each level of this index [uli]
  # but (g)lm complains about a model where uli has a single level [though this is meaningful]
  uli <- ULI(mf) 
  # selection of data for replicates
  mf <- data.frame(y=Y,uli=as.factor(uli)) ## what's needed for both sides
  tuli <- table(uli)
  phiinfo <- which(tuli>1L) ## which levels of uli are informative
  whichRows <- uli %in% phiinfo
  mf <- mf[whichRows,] ## keep lines with replicates of all predictor variables
  if (!is.null(weights)) weights <- weights[whichRows]
  if (length(phiinfo)>1L) {
    locform <- y ~ uli 
  } else locform <- y ~ 1 ## only one level of uli had replicates
  # estimation of residual var
  if (NROW(mf)>1L) {
    # formula must retain any operation on the lhs
    locglm <- glm(formula=locform,data=mf,family=family,weights=weights)  
    phi_est <- as.numeric(deviance(locglm)/locglm$df.residual)
  } else phi_est <- NULL
  return(phi_est)
}


preprocess <- function(control.HLfit, ranFix=NULL, HLmethod, 
                       predictor, resid.model,
                       REMLformula, data, family,
                       BinomialDen, rand.families, etaFix, prior.weights,
                       control.glm ) {
  callargs <- match.call() ## to make easy to change these arguments in a later fit
  if ( ! is.list(resid.model)) resid.model <- list(formula=resid.model,family=control.HLfit$resid.family)

  ################ handling list of data #######################
  if (inherits(data,"list")) {
    #     data <- lapply(data,function(dt) {
    #       validdata <- validData(formula=predictor,resid.formula=resid.model$formula,
    #                              data=dt) ## will remove rows with NA's in required variables
    #       dt[rownames(validdata),,drop=FALSE] ## ## before Predictor is called and an LMatrix is added, etc.
    #     })
    locargs <- as.list(callargs)
    famfam <- family$family
    processed <- lapply(data,function(dd) {
      locargs$data <- dd
      if ( ! is.null(famfam) && famfam=="multi") locargs$family <- family$binfamily  
      eval(as.call(locargs)) ## call("preprocess",...) on each data set
    })
    attr(processed,"multiple") <- TRUE ## but NULL otherwise hence test it as not null  
    return(processed)
  }
  ###############################################################
  
  # remove rows with NA's in required variables:
  validdata <- validData(formula=predictor,resid.formula=resid.model$formula,
                         data=data) 
  if (!inherits(data,"environment")) {
    data <- data[rownames(validdata),,drop=FALSE]  
  } else data <- validdata
  # add easily testable family name
  famfam <- family$family
  processed <- list(data=data,family=family)
  #
  stop.on.error <- control.HLfit$stop.on.error ##  
  if (is.null(stop.on.error)) stop.on.error <- FALSE
  processed$stop.on.error <- stop.on.error ##  
  break_conv_logL <- control.HLfit$break_conv_logL ##whether to stop if logL (p_v) appears to have converged  
  if (is.null(break_conv_logL)) break_conv_logL <- FALSE
  processed$break_conv_logL <- break_conv_logL ##  
  AIC <-control.HLfit$AIC ##  
  if (is.null(AIC)) AIC <- FALSE
  processed$AIC <- AIC ##  
  essai <-control.HLfit$essai ##  
  if (is.null(essai)) essai <- FALSE
  processed$essai <- essai ##  
  ## numerical control parameters 
  conv.threshold <-control.HLfit$conv.threshold ## 
  if (is.null(conv.threshold)) conv.threshold <- 1e-05
  processed$conv.threshold <- conv.threshold  
  #
  HL2.threshold <- control.HLfit$HL2.threshold 
  if (is.null(HL2.threshold)) HL2.threshold <- 3
  #    processed$HL2.threshold <- HL2.threshold  
  #
  iter.mean.dispFix <-control.HLfit$iter.mean.dispFix ## private control
  if (is.null(iter.mean.dispFix)) iter.mean.dispFix <- control.HLfit$max.iter.mean ## public control
  if (is.null(iter.mean.dispFix)) iter.mean.dispFix <- 200 ## control of inner loop when no disp param is estimated ## was 40, 06/2014
  processed$iter.mean.dispFix <- iter.mean.dispFix  
  #
  iter.mean.dispVar <-control.HLfit$iter.mean.dispVar ## private control
  if (is.null(iter.mean.dispVar)) iter.mean.dispVar <- control.HLfit$max.iter.mean ## public control ## control of inner loop when some disp param is estimated
  if (is.null(iter.mean.dispVar)) iter.mean.dispVar <- 50 ## control of inner loop when some disp param is estimated  ## was 20, 06/2014
  processed$iter.mean.dispVar <- iter.mean.dispVar  
  #
  max.iter <- control.HLfit$max.iter  ## control of outer loop 
  if (is.null(max.iter)) max.iter <-200
  processed$max.iter <- max.iter  
  #
  ## only check, no modif of data -> not returned
  if ( ! (inherits(data,"data.frame") || inherits(data,"environment"))) {
    mess <- pastefrom("'data' is not a data.frame not an environment.",prefix="(!) From ")
    stop(mess)
  }
  #
  resid.family <- resid.model$family
  if (is.null(resid.family)) {
    resid.family <- spaMM_Gamma(log)
  } else {
    if (resid.family$family!= "Gamma") stop("resid.family must be Gamma.")
    resid.family <- spaMM_Gamma(resid.family$link) ## we will need the returned link, not the promise 
  }
  processed$resid.family <- resid.family ## resid.predictor will also be returned
  #
  if (! inherits(predictor,"predictor")) predictor <- Predictor(predictor) 
  MeanFrames <- HLframes(formula=predictor,data=data) ## design matrix X, Y...
  #
  y <- MeanFrames$Y
  if ( family$family == "binomial" && NCOL(y)==2L) y <- y[,1L,drop=FALSE] ## that is, we have the cbind syntax up to this fn
  ## binomial denom determined later, not using y => FR->FR not clear why we need two column format up to this point.
  processed$y <- y
  nobs <- NROW(y)
  #
  X.pv <- MeanFrames$X
  ###   OFFSET
  off <- model.offset(MeanFrames$mf) ## look for offset from (ori)Formula 
  if ( is.null(off) ) { ## no offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
    offsetObj <- attr(predictor,"offsetObj")
    off <- offsetObj$offsetArg 
  } else {
    offsetObj <- list()
    keepOriForm <- attr(predictor,"oriFormula")
    predictor <- noOffset(predictor)
    attr(predictor,"oriFormula") <- keepOriForm
  }
  if (!is.null(off)) {
    off <- pmax(log(.Machine$double.xmin),off) ## handling log(0) ## but if input off were NULL, output off would be is numeric(0) where it should remain NULL
  }
  colnames(X.pv) <- colnames(MeanFrames$X)
  ## reimplementation of etaFix$beta (2015/03)
  X.Re <- X.pv ## can be overwritten according to etaFix and REMLformula
  namesOri <- colnames(X.pv)
  if ( length(betaFix <- etaFix$beta)>0 ) {
    namesbetafix <- names(betaFix)
    if (is.null(namesbetafix)) {
      message("(!) The elements of etaFix$beta should be named and the names should match the column names of the design matrix.")
      stop("    I exit.")
    }
    if (length(setdiff(namesbetafix,namesOri))==0L) { ## if no incorrect name
      offFromEtaFix <- X.pv[,namesbetafix,drop=FALSE] %*% betaFix
      namesbetavar <- setdiff(colnames(X.pv),namesbetafix)
      X.pv <- structure(X.pv[,namesbetavar,drop=FALSE],namesbetavar=namesbetavar)
      offsetObj$betaFix <- betaFix
      if (is.null(off)) {
        off <- offFromEtaFix
      } else off <- off + offFromEtaFix
      ## TRUE by default:
      if ( is.null( keepInREML <- attr(betaFix,"keepInREML") ) ||  ( ! keepInREML) ) X.Re <- X.pv ## can be overwritten according to REMLformula
    } else {
      message("(!) The names of elements of etaFix$beta should all match column names of the design matrix.")
      stop("    I exit.")
    }
  } 
  attr(X.pv,"namesOri") <- namesOri
  if (is.null(off)) { ## model.frame.default(formula = locform, offset = off,...) expects a vector....
    offsetObj$total <- rep(0,nobs)
    offsetObj$nonZeroInfo <- FALSE
  } else {
    offsetObj$total <- off
    offsetObj$nonZeroInfo <- TRUE
  }
  #
  attr(predictor,"offsetObj") <- offsetObj ## has $total, $nonZeroInfo, possibly $offsetArg (argument of Predictor()) and $betaFix. As nothing from the formula itself 
  processed$predictor <- predictor  
  processed$X.pv <- X.pv ## further modified in HLfit if non-estimable parameters
  # handling offset in *resid.predictor*
  resid.predictor <- resid.model$formula
  if (! inherits(resid.predictor,"predictor")) resid.predictor <- Predictor(resid.predictor)
  roff <- model.offset(model.frame(resid.predictor,data=data)) ## look for offset from (ori)Formula 
  if ( is.null(roff) ) { ## ## no offset (ori)Formula term. Check attribute (Predictor ensures offset is not both in formula and attr)
    roff <- attr(resid.predictor,"offsetObj")$total 
  } else {
    keepOriForm <- attr(resid.predictor,"oriFormula")
    resid.predictor <- noOffset(resid.predictor)
    attr(resid.predictor,"oriFormula") <- keepOriForm
  }
  if (is.null(roff)) { ## model.frame.default(formula = locform, offset = off,...) expects a vector....
    attr(resid.predictor,"offsetObj") <- list(total=rep(0,nobs),nonZeroInfo=FALSE)
  } else {
    attr(resid.predictor,"offsetObj") <- list(total=roff,nonZeroInfo=TRUE) ## subtly, will be of length 1 if   original offset was a constant...
  }
  processed$resid.predictor <- resid.predictor  
  #
  # comme l'EQL est dans un monde quasi gaussien, il se ramene aux LMM et utilise les leverage standard,
  # Pour les GLMM non LMM, Dans la mesure ou HL(0,.) utilise h et non q+, ce n'est pas l'EQL pour les fixed params
  # Dans la mesure ou on utilise les leverages standard pour les dispersion param, c'est de l'EQL
  ## glm() syntax:  
  ## With binomial data the response can be either a vector or a matrix with two columns.
  ## If the response is a vector, it is treated as a binary factor with the first level representing "success" and all others representing "failure". 
  ## In this case R generates a vector of ones to represent the binomial denominators.
  ## Alternatively, the response can be a matrix where the first column shows the number of "successes" and the second column shows the number of "failures". 
  ## In this case R adds the two columns together to produce the correct binomial denominator. 
  if (family$family=="binomial") {
    BinDenForm <- attr(predictor,"BinDenForm")
    if ( ! is.null(BinDenForm)) { ## the cbind syntax was used in the formula 
      BinomialDen <- eval(parse(text=BinDenForm),envir=data)
      ## la suite ducode suppose que pas cbind => ie essentially obsolete syntax 
    } else if (missing(BinomialDen) || is.null(BinomialDen)) { ## then this should be a binary response
      checkResp <- eval(parse(text=as.character(predictor[2])),envir=data) ## 'positives'
      if (any(checkResp>1)) {
        mess <- pastefrom("binomial, non-binary response. Please use the",prefix="(!) From ")
        message(mess)
        message("    standard glm() syntax with _cbind_: 'cbind(<successes>, <failures>) ~ <etc>'")
        stop()
      } else BinomialDen <- rep(1,nobs) ## response appears to be binary...
    }
    no.info <- (BinomialDen == 0)
    if (any(no.info)) {
      mess <- pastefrom("please remove missing data (i.e. for which binomial sample size is 0).",prefix="(!) From ")
      stop(mess)
    }
    ## It's not really possible to remove data at this stage as this may not match the dimension of the distance matrices
    ## moreover one cannot simply remove rows of a matrix "root"...
    ## it _could_ be useful to be able to hand BinomilaDen=0 by the general code but...
    if (var(y)==0 && var(BinomialDen)==0) {
      warning("var(response) = 0.")
    }  
  } else {
    BinomialDen <- rep(1,nobs)
    if (var(y)==0) {
      warning("var(response) = 0.")
    }  
  }
  processed$BinomialDen <- BinomialDen  
  ## conversion from user-friendly format to standard 'XX(...)' format
  ## first index is for (0) h / (1) p_v(h) / (2) p^s_v(h) ie whether h lik or marginal lik is used for fixed effect estimation
  ## Other components determine three options wrt to leverages, some for ReML correction, some for notEQL.
  ## ML() vs HL() determines whether the hatval leverages are computed ie whether some basic ReML correction in applied
  ## second index is for further ReML correction: no (0) or yes (1) D log Det / d log lambda correction (2) further p^s_bv correction
  ## third index is for use of (0) EQL deviance residuals (this affects the leverages) or not (1) (this is not an ReML correction... but impatcs only dispersion estimation..)
  ## thus overall we have <ReML/not>( <h/l> , <more ReML/not> , <not/EQL> )
  ## NohL07 table 1 has interesting terminology and further tables show even more cases
  terms_ranefs <- parseBars(predictor) ## a vector of char strings with attribute(s), 
  ##    otherwise similar to bars <- spMMexpandSlash(findbarsMM(formula[[length(formula)]])): FR->FR maybe simplification of code possible here?
  models <- list(eta="",lambda="",phi="")
  if ( length(terms_ranefs) > 0L ) {
    processed$lambdaFamily <- Gamma(link="log")
    models[["eta"]] <- "etaHGLM" 
    FL <- spMMFactorList(predictor, MeanFrames$mf, 0L, drop=TRUE) ## this uses the spatial information in the formula, even if an explicit distMatrix was used elsewhere
    ZAlist <- FL$Design ## : is a list of design matrices (temporarily only Z)
    attr(ZAlist,"ranefs") <- terms_ranefs
    attr(ZAlist,"namesTerms") <- FL$namesTerms
    AMatrix <- attr(predictor,"AMatrix")
    if (!is.null(AMatrix)) {
      ## logic is Z[nresp,nUniqueRespLoc].A[nUniqueRespLoc,nHiddenv].L[nHiddenv,nHiddenv]
      for (iMat in seq(length(ZAlist))) {
        ZAlist[[iMat]] <- ZAlist[[iMat]] %*% AMatrix[[iMat]]
      }
    }
  } else {
    models[["eta"]] <- "etaGLM" 
    ZAlist <- structure(list(),anyRandomSlope=FALSE)
  }
  processed$ZAlist <- ZAlist ## attributes will be added later !
  nrand <- length(ZAlist)
  #
  if (inherits(rand.families,"family")) rand.families <- list(rand.families) ## I should pass rand.families to preprocess
  if (nrand != 1L && length(rand.families)==1L) rand.families <- rep(rand.families,nrand) 
  rfblob <- checkRandLinkS(rand.families)  
  rand.families <- rfblob$rand.families
  lcrandfamfam <- rfblob$lcrandfamfam
  if (HLmethod=="ML") {
    HLmethod <- "ML(1,1,1)" ## here there could be a special hack for (family$family=="binomial" && mean(BinomialDen)<HL2.threshold) 
  } else if (HLmethod=="SEM") {
    HLmethod <- "ML('SEM',NA,NA)" 
    if ( ! (family$family=="binomial" && family$link=="probit")) {
      stop("SEM is applicable only to binomial(probit) models.")
    }
    if (sum(BinomialDen) != nobs) stop("(!) SEM procedure: the data do not seem binary; non-binary data are not handled.")
    SEMseed <- control.HLfit$SEMseed ##  
    if (is.null(SEMseed)) SEMseed <- 123 ## OK pour SEM *unique* mais remplace par NULL dans optimthroughSmooth
    SEMargs <- list(SEMseed=SEMseed)
    SEMargs$nMCint <- control.HLfit$nMCint ##  as is => SEM procedure must handle null value
    SEMargs$control_pmvnorm$maxpts <- control.HLfit$pmvnorm_maxpts ##  as is => SEM procedure must handle null value
    SEMlogL <- control.HLfit$SEMlogL
    if (is.null(SEMlogL)) {
      SEMlogL <- "pmvnorm"
    }
    SEMargs$SEMlogL <- SEMlogL
    nSEMiter <- control.HLfit$nSEMiter ##  
    if ( (! is.null(nSEMiter)) && nSEMiter < 10) {
      stop(" 'nSEMiter' should be >9")
    } else SEMargs$nSEMiter <- nSEMiter
    SEMargs$ngibbs <-control.HLfit$ngibbs ##  
    SEMargs$SEMsample <- control.HLfit$SEMsample ## stays NULL if NULL
    processed$SEMargs <- SEMargs
  } else if (HLmethod=="REML") {
    HLmethod <- "HL(1,1,1)" ## here there could be a special hack for (family$family=="binomial" && mean(BinomialDen)<HL2.threshold) 
  } else if (HLmethod=="REPQL" || HLmethod=="PQL") {
    if (any(lcrandfamfam!="gaussian"))  stop("PQL is not defined for HGLMs in general. Do you mean 'EQL-'?") 
    HLmethod <- "HL(0,0,1)" ## (=REPQL, equivalent to HL(0,0,0) ie EQL- for GLMMs )
  } else if (HLmethod=="PQL/L") { ## again no D log Det / d log lambda correction
    if (any(lcrandfamfam!="gaussian"))  stop("PQL is not defined for HGLMs in general. Do you mean 'EQL-'?") 
    HLmethod <- "ML(0,0,1)" ## (equivalent to ML(0,0,0) for GLMMs)
  } else if (HLmethod=="EQL-") { ## version LeeNP06 p.212 incomplete 1st order ## probably hglm package
    ## thus overall we have <ReML->HL >( <h->0> , <not more ReML->0> , <EQL -> 0> )
    HLmethod <- "HL(0,0,0)" ## 
  } else if (HLmethod=="EQL+") { ## version LeeN01 complete 1st order
    ## thus overall we have <ReML->HL >( <h->0> , <more ReML->1> , <EQL -> 0> )
    HLmethod <- "HL(0,1,0)" ## (0,...): gaussianise everything, hence no a(1) correction ## there is no HL(1) a(1) correction in GLM.MME
  }
  HL <- eval(parse(text=paste("c",substr(HLmethod,3,100),sep=""))) ## extracts the (...) part into a vector
  if (length(HL)==2) HL <- c(HL,1)
  ## HLmethod is not a member of the return object
  processed$HL <- HL  
  if (substr(HLmethod,0,2)=="ML") { # && HL[1]!="SEM") { ## FR->FR c'est bizarre d'exclure le SEM là... p_bv est il vraiment utilisé ?
    if ( ! is.null(REMLformula)) {
      message("Confusing combination of arguments: 'HLmethod=ML(...)' with non-null 'REMLformula'.")
      stop("  Make sure what you mean and simplify the arguments.")
    }
    lhs <- paste(predictor)[[2]] ## extract response, either cbind or not
    if ( ! is.null(terms_ranefs)) {
      ## build formula with only random effects
      REMLformula <- as.formula(paste(lhs,"~",paste(terms_ranefs,collapse="+")))
    } else REMLformula <- as.formula(paste(lhs,"~ 0")) 
  }
  processed$REMLformula <- REMLformula  
  #
  processed$loglfn.fix <- selectLoglfn(family)
  if ( family$family %in% c("binomial","poisson","COMPoisson")) {
    ## the response variable should always be Counts
    if (max(abs(y-as.integer(y)))>1e-05) {
      mess <- pastefrom("response variable should be integral values.",prefix="(!) From ")
      stop(mess)
    }
    if (DEPARSE(resid.predictor) != "~1") {
      warning(paste("resid.model is ignored in ",family$family,"-response models",sep=""))
    }
  }
  ## code derived from the glm() function in the safeBinaryRegression package
  if(family$family == "binomial" && length(unique(y)) == 2L && ncol(X.pv)>0L) {
    separation <- separator(X.pv, as.numeric(y), purpose = "test")$separation
    if(separation) {
      message("Separation exists among the sample points.\n\tThis model cannot be fit by maximum likelihood.")
      message("The following terms are causing separation among the sample points:")
      separation <- separator(X.pv, as.numeric(y), purpose = "find")$beta
      separating.terms <- dimnames(X.pv)[[2]][abs(separation) > 1e-09]
      if(length(separating.terms)) message(paste(separating.terms, collapse = ", "))
      stop()
    }
  }
  #
  if ( ! is.null(REMLformula) ) { ## differences affects only REML estimation of dispersion params, ie which p_bv is computed
    REMLFrames <- HLframes(formula=REMLformula,data=data) ## design matrix X, Y...
    X.Re <- REMLFrames$X  
  } ## else keep previously computed X.Re 
  processed$X.Re <- X.Re
  #
  canonicalLink <- FALSE
  if (family$family=="gaussian" && family$link=="identity") {
    canonicalLink <- TRUE
  } else if (family$family=="poisson" && family$link=="log") {
    canonicalLink <- TRUE
  } else if (family$family=="binomial" && family$link=="logit") {
    canonicalLink <- TRUE
  } else if (family$family=="Gamma" && family$link=="inverse") {
    canonicalLink <- TRUE
  } else if (family$family=="COMPoisson" && family$link=="loglambda") {
    canonicalLink <- TRUE
  }
  processed$canonicalLink <- canonicalLink  
  #
  GLMMbool <- (nrand>0 && all(lcrandfamfam=="gaussian") )
  if (GLMMbool) {
    LMMbool <- (family$family=="gaussian" && canonicalLink)
  } else LMMbool <- FALSE
  processed$LMMbool <- LMMbool
  processed$GLMMbool <- GLMMbool
  #
  if (is.null(prior.weights)) prior.weights <- rep(1,nobs)
  attr(prior.weights,"unique") <- (length(unique(prior.weights))==1L) 
  processed$prior.weights <- prior.weights
  ## algorithms (control of defaults remains in the HLfit code)
  betaFirst <-control.HLfit$betaFirst ##  
  if (is.null(betaFirst)) {
    betaFirst <- FALSE
  } else if (betaFirst && HL[1]=="SEM") {
    message("betaFirst && HLmethod= SEM: betaFirst turned to FALSE")
    betaFirst <- FALSE
  }
  processed$betaFirst <- betaFirst ##
  LevenbergM <- control.HLfit$LevenbergM ## 
  if (is.null(LevenbergM)) { 
    if (HL[1]=="SEM") {
      LevenbergM <- FALSE 
    } else if (betaFirst) {
      LevenbergM <- FALSE 
    } else {
      if (LMMbool) {  
        LevenbergM <- FALSE ## because no reweighting when beta_eta changes => no IRWLS necess   
      } else LevenbergM <- .spaMM.data$options$LevenbergM
    }
  }
  processed$LevenbergM <- LevenbergM
  #
  if (nrand>0) {   
    processed$lambda.Fix <- getPar(ranFix,"lambda")
    nrand_lambda <- 0
    models[["lambda"]] <- FL$termsModels
    ################################################################################
    # for a random slope term, ie v= v_1+x v_2 , the x went into the general ZAL matrix 
    # (construction of ZAlist by spMMFactorList), and
    # we are still estimating the lambda's using a X_lamres with 0/1's only
    # unless there is a non-trivial model for the lambdas
    ################################################################################
    if (all(models[["lambda"]]=="lamScal")) { ## all mixed models handled in 06/2015 (random slope, adjacency...) 
      Xi_cols <- unlist(lapply(FL$namesTerms,length))
      if (any(Xi_cols>1 & !lcrandfamfam=="gaussian")) {
        stop("(!) random slope models with correlated non-gaussian random effects are not fitted.")
      }
      cum_Xi_cols <- cumsum(c(0, Xi_cols)) ## if two ranef,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
      n_u_h <- rep(0, sum(Xi_cols))
      for (i in 1:nrand) n_u_h[(cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]] <-  nlevels(FL$Subject[[i]]) 
      cum_h_u_h <- cumsum(c(0, n_u_h)) ## if two "Intercept" ranefs,  with n_u_h=(3,3), this is 0,3,6. cum_h_u_h[nrand+1] is then 6, the total # of realizations
        ## if (1+X|...) +(1|...),  with n_u_h=(3,4), this is 0,3,6,10. cum_h_u_h[sum(Xi_cols)+1] is then 10, the total # of realizations
      X_lamres <- matrix(0,cum_h_u_h[sum(Xi_cols)+1L],sum(Xi_cols))
      colnames(X_lamres) <- unlist(FL$namesTerms)
      for (i in seq(nrand)) {
        for (j in (cum_Xi_cols[i]+1L):cum_Xi_cols[i+1L]) {
          X_lamres[(cum_h_u_h[j]+1L):cum_h_u_h[j+1L],j] <- 1L ## this maps the deviance residuals to the lambda's to be estimated from them. None of the random-slope columns is a constant full column because each dev res is used for estimating only one lambda. Nevertheless, the first col will be called "(Intercept)", and this makes a valid output.
        }
      } 
    } else {  ## linear predictor for variance of random effects (lambda) (lamGLM or lamHGLM) 
      if (any(models[["lambda"]]=="lamHGLM")) { ##need distinct calls... to fit each lambda model  
        if (length(formulaLambda)==2) formulaLambda <- as.formula(paste('"lambda"',paste(formulaLambda,collapse=" ")))
        lambda_ranefs <- findbarsMM(formulaLambda)
        if (!is.null(lambda_ranefs)) {  ## lamHGLM
          nrand_lambda <- length(lambda_ranefs)
          models[["lambda"]] <- list("lamHGLM")
        } else models[["lambda"]] <- list("lamGLM")  
        colnames(X_lambda) <- colnames(fr_lambda$X) ## but code not effective, fr_lambda not computed
      } else { ## can use a single design matrix for all random effects, which is convenient.
        mess <- pastefrom("LIKELY missing code to handle linear predictor for lambda.")
        stop(mess)
        # la suite c'est dexu residus de code a assembler: il faut une liste de HLframes du type
        fr_lambda <- HLframes(formula=formulaLambda,data=data) ## but the "data" should probably be distinct data here, with nrow=number of reals of ranefs 
        # (pobablement calculee en amont pour determiner lamScal aussi...) ensuite extraire les design matrices
        #X_lamres ? Xi_cols ?
      }
    } 
    processed$X_lamres <- X_lamres ## for glm for lambda
    attr(processed$ZAlist,"Xi_cols") <- Xi_cols ## used to handle random slope...
    attr(processed$ZAlist,"anyRandomSlope") <- any(Xi_cols>1L) ## used to handle random slope...
  } 
  #
  phi.Fix <- getPar(ranFix,"phi")
  if (is.null(phi.Fix)) {
    if (family$family %in% c("poisson","binomial","COMPoisson")) phi.Fix <- 1 
  } else if (any(phi.Fix==0)) stop("phi cannot be fixed to 0.")
  processed$phi.Fix <- phi.Fix
  if ( is.null(phi.Fix)) {
    formulaDisp <- resid.predictor
    fr_disp <- HLframes(formula=formulaDisp,data=data) 
    X_disp <- fr_disp$X
    ## if formula= ~1 and data is an environment, there is no info about nobs, => X_disp has zero rows, which is a problem later 
    if(nrow(X_disp)==0L) X_disp=matrix(1,nrow=nobs)
    namesX_disp <- colnames(fr_disp$X)
    colnames(X_disp) <- namesX_disp
    random_dispersion<-findbarsMM(formulaDisp) ## random effect in mean predictor of dispersion phi
    if (!is.null(random_dispersion)) {
      FL_disp <- spMMFactorList(formulaDisp, fr_disp$mf, 0L, drop=TRUE)
      Z_disp <- FL_disp$Design  ## now Matrix...
      namesRE_disp <- attr(Z_disp,"Groupings")
      models[["phi"]] <- "phiHGLM"
      mess <- pastefrom("LIKELY missing code to handle random effect for linear predictor for phi.")
      stop(mess)
    } else {
      if (length(namesX_disp)==1 && namesX_disp[1]=="(Intercept)") {
        models[["phi"]] <- "phiScal"
        processed$init_phi <- generateInitPhi(formula=predictor,data=data,family=family,weights=prior.weights)
      } else { 
        models[["phi"]] <- "phiGLM"
      }
    } 
    processed$X_disp <- X_disp
  } 
  #
  if (LMMbool) {
    ## identifiability checks cf modular.R -> checkNlevels() in lmer:
    ## FR->FR one could also add an isNested check as in https://github.com/lme4/lme4/blob/master/R/utilities.R, but presumably not here
    LMatrix <- attr(predictor,"LMatrix")
    for (iMat in seq(length(ZAlist))) {
      nc <- ncol(ZAlist[[iMat]])
      if (nc < 2 && is.null(phi.Fix)) {
        mess <- paste("Only ",nc," level for random effect ",terms_ranefs[iMat],
                      ";\n   this model cannot be fit unless phi is fixed.",sep="")
        warning(mess)
      }
      if (is.null(LMatrix) && substr(terms_ranefs[iMat], 1, 1)=="(" && nc == nobs && is.null(phi.Fix)) {
        mess <- paste("Number of levels = number of observations \n   for random effect ",
                      terms_ranefs[iMat],
                      ";\n   this model cannot be fit unless phi is fixed\n   or a correlation matrix is given.",sep="")
        stop(mess)
      }
    }
  } 
  #
  processed$models <- models
  attr(rand.families,"lcrandfamfam") <- lcrandfamfam
  attr(rand.families,"unique.psi_M") <- sapply(lcrandfamfam, function(v) {
    switch(v, 
           gaussian = 0,
           gamma = 1, 
           beta = 1/2, 
           "inverse.gamma" = 1
    )
  })
  processed$rand.families <- rand.families
  processed$callargs <- callargs
  processed$fixef_terms <- MeanFrames$fixef_terms ## added 2015/12/09 for predict
  processed$fixef_levels <- MeanFrames$fixef_levels ## added 2015/12/09 for predict
  processed$control.glm <- do.call("glm.control", control.glm) ## added 04/2016
  class(processed) <- c("arglist","list")
  return(processed)
}

eval.update.call <- function(mc,...) {
  mc <- as.list(mc)
  dotlist <- list(...)
  mc[names(dotlist)] <- dotlist ## a un moment j'ai mis cette ligne en commentaire, ce qui rend la fonction ineffective !
  eval(as.call(mc))  
}
