##
## This file contains MODIFIED COPIES of the lm, summary.lm, confint.lm and anova.lm methods from the base package stats (2014-10-31),
## and Anova.lm from the package car (2014-10-10).
##

# The lm function replaces stats::lm, organizes fixed and random effects, removes r() from formula and parses to lm or lmer.

anova.lmm <- function(object, ...){
  if(!is.null(object$random)){
    stop("Only type III error Anova() supported in mixed/random models")
  } else {
    class(object) <- 'lm'
    return(anova(object,...))
  }
}


Anova.lmm <- function(mod, ...){
  mf <- match.call()
  if(!is.null(mod$random)){
    if(mf$type=="II" || mf$type=="2"){
      warning("Using (default) type III error supported for mixed/random models")
    }
    return(AnovaMix(mod))
  } else {
    class(mod) <- 'lm'
    return(Anova(mod, ...))
  }
}


## FIXME: Her antas modellen å kun inneholde faktorer, ingen andre effekter. Kan gi rare resultater!
# Mixed model ANOVA
AnovaMix <- function(object){
  formula         <- formula(object)
  formula.text    <- as.character(formula)
  all.effects     <- object$random$all							  # All model effects and interactions
  fixed.effects   <- object$random$fixed							  # All fixed effects
  random.effects  <- object$random$random						  # All random effects
  main.rands.only.inter <- object$random$main.rands.only.inter     # Random effects only present in interactions
  restrictedModel <- !object$random$unrestricted
  data    <- object$model
  n.effects    <- length(all.effects)
  main.effects <- fparse(formula)							  # All main effects (even though only included in interactions)
  n.levels     <- numeric(length(main.effects))
  for(i in 1:length(main.effects)){
    n.levels[i] <- length(levels(data[,main.effects[i]])) # Number of levels per main effect
  }
  names(n.levels) <- main.effects
  N <- dim(data)[1]
  
  ind.randoms <- numeric()
  ind.randoms <- match(random.effects,all.effects) # Placement of random effects in "all.effects"
  ind.fixed   <- match(fixed.effects,all.effects)  # Placement of fixed effects in "all.effects"
  ind.fixed <- setdiff(1:n.effects,ind.randoms)										
  n.randoms    <- length(ind.randoms)
  
  # Estimate fixed effect Anova
  opt <- options("contrasts")
  options(contrasts=c('contr.sum','contr.poly'))
  noRandom <- update(object)
  noRandom$random <- NULL
  class(noRandom) <- "lm"
  fixed.model <- as.data.frame(Anova(noRandom, type='III', singular.ok=TRUE))
  options(contrasts=opt$contrasts)
  # 	fixed.model <- fixed.model[-1,] # Remove intercept
  fixed.model <- fixed.model[c(all.effects,"Residuals"),] # Sort according to all.effects
  if(!any("Mean Sq"%in%colnames(fixed.model))){
    fixed.model <- cbind(fixed.model[,"Sum Sq"]/fixed.model[,"Df"], fixed.model)
    colnames(fixed.model)[1] <- "Mean Sq"
  }
  
  # Check which effects should use interactions as denominators instead of error
  approved.interactions <- list()
  approved.interactions.fixed <- list()
  for(i in 1:n.effects){
    this.effect <- strsplit(all.effects[i],":")[[1]]
    which.contains <- numeric()
    for(j in 1:n.effects){ # Find all other effects containing this.effect
      effect.names <- is.element(strsplit(all.effects[j],":")[[1]],this.effect)
      if(i!=j && sum(effect.names)==length(this.effect) && length(effect.names)>length(this.effect)){
        which.contains <- union(which.contains,j)}
    }
    which.contains <- sort(which.contains)
    if(length(which.contains)>0){
      approved.interaction <- numeric(length(which.contains))
      approved.interaction.fixed <- numeric(length(which.contains))
      for(j in 1:length(which.contains)){
        if(restrictedModel){
          approved.interaction[j] <- prod(is.element(setdiff(strsplit(all.effects[which.contains],":")[[j]],strsplit(all.effects[i],":")[[1]]),c(random.effects,main.rands.only.inter)))}
        else{
          if(any(is.element(ind.fixed,i))){
            approved.interaction.fixed[j] <- prod(is.element(setdiff(strsplit(all.effects[which.contains],":")[[j]],strsplit(all.effects[i],":")[[1]]),fixed.effects))}
          approved.interaction[j] <- 1-prod(!is.element(strsplit(all.effects[which.contains],":")[[j]],c(random.effects,main.rands.only.inter)))}
      }
      if(length(which(approved.interaction==1))>0){
        approved.interactions[[i]] <- which.contains[which(approved.interaction==1)]}
      else{
        approved.interactions[[i]] <- FALSE}
      if(length(which(approved.interaction.fixed==1))>0){
        approved.interactions.fixed[[i]] <- which.contains[which(approved.interaction.fixed==1)]}
      else{
        approved.interactions.fixed[[i]] <- FALSE}
    }
    else{
      approved.interactions[[i]] <- FALSE
      approved.interactions.fixed[[i]] <- FALSE}
  }
  
  # Find variance components (except MSerror), 
  # and find linear combinations needed to produce denominators of F-statistics
  mix.model.attr <- list()
  denom.df <- numeric(n.effects+1)
  exp.mean.sq <- rep(paste("(",n.effects+1,")", sep=""), n.effects+1)
  var.comps <- numeric(n.effects+1)*NA
  var.comps[n.effects+1] <- fixed.model[n.effects+1,"Mean Sq"]
  errors <- numeric(n.effects)
  for(i in 1:n.effects) {
    if(!is.logical(approved.interactions[[i]])){
      # Set up matrix A and vector b to find linear combinations of effects to use as denominators in F statistics
      ## This is probably where unbalancedness should be included !!!!!!!!
      lap <- length(approved.interactions[[i]])
      A <- matrix(0,lap+1,n.effects+1)
      b <- rep(1,lap+1)
      for(j in 1:lap){
        A[j,approved.interactions[[approved.interactions[[i]][j]]]] <- 1
        A[j,approved.interactions[[i]][j]] <- 1
        k <- length(approved.interactions[[i]])+1-j
        exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[approved.interactions[[i]][k]],":")[[1]]]), " (",which(all.effects==all.effects[approved.interactions[[i]][k]]),")", sep="")
      }
      A[, n.effects+1] <- 1
      A <- A[,apply(A,2,sum)>0]
      denominator <- solve(t(A),b)
      denominator.id <- c(approved.interactions[[i]],n.effects+1)
      denominator.id <- denominator.id[denominator!=0]
      mix.model.attr[[i]] <- denominator <- denominator[denominator!=0]
      names(mix.model.attr[[i]]) <- denominator.id
      if(length(denominator)==1){ # Original df
        denom.df[i] <- fixed.model[denominator.id,"Df"]}
      else{ # Satterthwaite's df correction
        denom.df[i] <- sum(fixed.model[denominator.id,"Mean Sq"]*denominator)^2/sum((fixed.model[denominator.id,"Mean Sq"]*denominator)^2/fixed.model[denominator.id,"Df"])} 
    } else{
      denominator.id <- n.effects+1
      mix.model.attr[[i]] <- 1
      names(mix.model.attr[[i]]) <- denominator.id
      denom.df[i] <- fixed.model[denominator.id,"Df"]
      denominator <- 1
    }
    if(sum(ind.randoms==i)>0){
      exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]), " (",i,")", sep="")
      var.comps[i] <- (fixed.model[i,"Mean Sq"]-fixed.model[denominator.id,"Mean Sq"]%*%denominator)/(N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]))}
    else{
      if(!is.logical(approved.interactions.fixed[[i]])){
        ex.ind <- paste(",", paste(approved.interactions.fixed[[i]], sep="", collapse=","),sep="")}
      else{
        ex.ind <- ""}
      exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]), " Q[",i,ex.ind,"]", sep="")
    }
    errors[i] <- fixed.model[denominator.id,"Mean Sq"]%*%denominator
    fixed.model[i,"F value"] <- fixed.model[i,"Mean Sq"]/(fixed.model[denominator.id,"Mean Sq"]%*%denominator)
    if(is.na(fixed.model[i,"F value"]) || fixed.model[i,"F value"]<0){
      fixed.model[i,"F value"] <- NA
    }
    fixed.model[i,"Pr(>F)"] <- 1-pf(fixed.model[i,"F value"],fixed.model[i,"Df"],denom.df[i])
  }
  names(denom.df) <- rownames(fixed.model)
  object <- list(lm=object, anova=fixed.model, err.terms=c(mix.model.attr,NA), denom.df=denom.df, restricted=restrictedModel,
                 exp.mean.sq=exp.mean.sq, var.comps=var.comps, random.effects=random.effects, ind.randoms=ind.randoms, formula.text=formula.text, errors=errors)
  class(object) <- "AnovaMix"
  object
}


## Print method for object from AnovaMix
print.AnovaMix <- function(x,...){
  object <- x
  N <- length(object$err.terms)
  output1 <- object$anova
  Fs <- PrF <- character(N)
  PrF[!is.na(output1$"Pr(>F)")] <- format(round(output1$"Pr(>F)"[!is.na(output1$"Pr(>F)")],4), digits=1, scientific=FALSE, nsmall=4)
  PrF[is.na(output1$"Pr(>F)")] <- "-"
  output1$"Pr(>F)" <- PrF
  Fs[!is.na(output1$"F value")] <- format(output1$"F value"[!is.na(output1$"F value")], digits=1, scientific=FALSE, nsmall=2)
  Fs[is.na(output1$"F value")] <- "-"
  output1$"F value" <- Fs
  output1$"Sum Sq" <- format(output1$"Sum Sq", digits=1, scientific=FALSE, nsmall=2)
  output1$"Mean Sq" <- format(output1$"Mean Sq", digits=1, scientific=FALSE, nsmall=2)
  
  err.terms <- character(length(object$err.terms))
  for(i in 1:N){
    if(length(object$err.terms[[i]])==1 && is.na(object$err.terms[[i]])){
      err.terms[i] <- "-"
    }
    else{
      err.terms[i] <- paste(ifelse(object$err.terms[[i]][1]>1,paste(object$err.terms[[i]][1],"*",sep=""),""),"(",names(object$err.terms[[i]][1]),")",sep="")
      if(length(object$err.terms[[i]])>1){
        for(j in 2:length(object$err.terms[[i]])){
          if(object$err.terms[[i]][j]<0){
            err.terms[i] <- paste(err.terms[i], paste(ifelse(object$err.terms[[i]][j]<(-1),paste(abs(object$err.terms[[i]][j]),"*",sep=""),""), "(", names(object$err.terms[[i]][j]), ")", sep=""), sep=" - ")
          } else {
            err.terms[i] <- paste(err.terms[i], paste(ifelse(object$err.terms[[i]][j]>1,paste(object$err.terms[[i]][j],"*",sep=""),""), "(", names(object$err.terms[[i]][j]), ")", sep=""), sep=" + ")}
        }
      }
    }
  }
  var.comps <- format(object$var.comps, digits=3)
  var.comps[setdiff(1:(N-1), object$ind.randoms)] <- "fixed"
  
  denom.df <- character(N)
  denom.df[!is.na(object$denom.df)&round(object$denom.df)==object$denom.df] <- format(object$denom.df[!is.na(object$denom.df)&round(object$denom.df)==object$denom.df], digits=3)
  denom.df[!is.na(object$denom.df)&round(object$denom.df)!=object$denom.df] <- format(object$denom.df[!is.na(object$denom.df)&round(object$denom.df)!=object$denom.df], digits=3)
  denom.df[object$denom.df==0] <- "-"
  output2 <- data.frame("Err.terms"=err.terms, "Denom.df"=denom.df, "VC(SS)"=var.comps)
  colnames(output2) <- c("Err.term(s)", "Err.df", "VC(SS)")
  output3 <- data.frame("E(MS)"=format(object$exp.mean.sq))
  colnames(output3) <- "Expected mean squares"
  rownames(output2) <- paste(1:N," ",rownames(object$anova), sep="")
  rownames(output3) <- rownames(object$anova)
  if(!object$restricted){
    un <- "un"}
  else{
    un <- ""}
  cat("Analysis of variance (", un, "restricted model)\n", sep="")
  cat("Response: ", object$formula.text[2], "\n", sep="")
  print(format(output1, digits=3))
  cat("\n")
  print(output2)
  cat("(VC = variance component)\n\n")
  print(output3)
  if(!is.balanced(object$lm)){
    cat("\nWARNING: Unbalanced data may lead to poor estimates\n")
  }
}

# lm from stats, edited to use treatment names in sum contrasts
# and enable limited classical least squares mixed models
lmer <- lme4::lmer
lm <- function (formula, data, subset, weights, na.action,
                method = "qr", model = TRUE, x = FALSE, y = FALSE,
                qr = TRUE, singular.ok = TRUE, contrasts = NULL,
                offset, unrestricted = TRUE, REML = NULL, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  ## Edited by KHL
  mfd <- match(c("formula","data"), names(mf), 0L)
  if(length(mfd)==2){ # Has formula and data
    is.random <- TRUE
    if( any(grepl("r(",formula,fixed=TRUE)) ){
      rw <- random.worker(formula, data, REML)
    } else {
      rw <- list(0)
    }
    if(length(rw) == 1){
      is.random <- FALSE
    } else { # Removed r() from formula
      formula <- rw$formula
      mf$formula <- rw$formula
      rw$unrestricted <- unrestricted
      if(is.logical(REML)){ # Perform 
        cl[[1]] <- as.name("lmer")
        cl[["formula"]] <- rw$reml.formula
        object <- eval(cl,parent.frame())
        object@call <- cl
        return(object)
      }
    }
  } else {
    is.random <- FALSE
  }
  ## End of edit
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
            domain = NA)
  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  w <- as.vector(model.weights(mf))
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y))
      matrix(,0,3) else numeric(), residuals = y,
      fitted.values = 0 * y, weights = w, rank = 0L,
      df.residual = if(!is.null(w)) sum(w != 0) else
        if (is.matrix(y)) nrow(y) else length(y))
    if(!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    ## Edited by KHL
    if(is.null(contrasts) && (options("contrasts")[[1]][1]!="contr.treatment" || options("contrasts")[[1]][1]!="contr.poly") && !missing(data)){
      col.names   <- effect.labels(mt,mf) # mt er "terms" fra formula, x er model.matrix
      if(length(col.names)==length(colnames(x))){
        colnames(x) <- effect.labels(mt,mf)
        effect.sources <- effect.source(mt,data)
      }
    }
    ## End edit
    z <- if(is.null(w)) lm.fit(x, y, offset = offset,
                               singular.ok=singular.ok, ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
  }
  if(is.matrix(y)){
    class(z) <- c( "mlm","lm")
  } else {
    class(z) <- c( "lmm", "lm")
  }
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr) z$qr <- NULL
  ## Edited by KHL
  if( is.random ){
    z$random <- rw
    if(!all(grepl("factor",attr(mt,"dataClasses")[-1])|grepl("ordered",attr(mt,"dataClasses")[-1]))){
      stop("Mixed models containing continuous effects not supported")
    }
  }
  if(exists("effect.sources") && !is.null(effect.sources))
    z$effect.sources <- effect.sources
  ## End edit
  z
}

## Collect and extract randomness
random.worker <- function(formula, data, REML = NULL){
  formula <- formula(formula)
  terms <- terms(formula)
  effsr <- attr(terms,"term.labels")
  effs  <- attr(terms(rparse(formula)),"term.labels")
  if(length(effs)==0){
    return( list(0) )
  }
  
  has.intercept <- attr(terms,"intercept")==1
  rands <- sort(unique(c(grep("[:]r[(]",effsr),   # Match random in interaction
                         grep("^r[(]",  effsr),   # Match random in the beginning
                         grep("[(]r[(]",effsr)))) # Match random inside function
  
  # which.rands <- match(rands,effsr)
  eff.splits <- list()
  for(i in 1:length(effs)){ # Split effect to look for hidden random interactions
    eff.splits[[i]] <- fparse(formula(paste("1~", effs[i],sep="")))
  }
  eff.lengths <- lapply(eff.splits,length)
  main.effs   <- effs[eff.lengths==1]
  main.rands  <- main.effs[main.effs%in%effs[rands]]
  main.rands.only.inter <- character(0)
  for(i in rands){
    main.rands.only.inter <- c(main.rands.only.inter, setdiff(eff.splits[[i]],main.effs)) # Random main effects only present in interactions
  }
  inter.rands <- which(unlist(lapply(eff.splits,function(i) any(main.rands%in%i))))
  # Check if any interactions containing random effects are not labeled as random
  if(any(is.na(match(inter.rands,rands)))){
    extra.randoms <- inter.rands[which(is.na(match(inter.rands,rands)))]
    warning(paste(paste(effs[extra.randoms],sep="",collapse=", "), " included as random interaction",ifelse(length(extra.randoms)==1,"","s"),sep=""))
    rands <- cbind(rands,extra.randoms)
    effs  <- effs[!(extra.randoms%in%effs)]
  }
  if(length(rands)==0){
    return( list(0) ) 
  } else {
    if(is.logical(REML)){
      remleffs     <- c(effs[setdiff(1:length(effs),rands)],paste("(1|",effs[rands],")",sep=""))
      reml.formula <- formula(paste(formula[[2]],"~",paste(remleffs,collapse="+"),ifelse(has.intercept,"","-1"),sep=""))
      
      return( list(formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")), random = effs[rands], main.rands.only.inter = main.rands.only.inter, fixed = effs[setdiff(1:length(effsr),rands)], all = effs, has.intercept = has.intercept, remleffs = remleffs, reml.formula = reml.formula))
    } else {
      return( list(formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")), random = effs[rands], main.rands.only.inter = main.rands.only.inter, fixed = effs[setdiff(1:length(effsr),rands)], all = effs, has.intercept = has.intercept))
    }
  }
}


###########################################
# summary.lm from stats, edited to enable limited classical least squares mixed models
summary.lmm <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate", 
                                               "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms)) 
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(object, "lm")) 
    warning("calling summary.lm(<fake-lm-object>) ...")
  Qr <- qr.lmm(object)
  n <- NROW(Qr$qr)
  if (is.na(z$df.residual) || n - p != z$df.residual) 
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept")) 
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  p1 <- 1L:p # Moved up for the sake of "random"
  if(is.null(object$random)){
    resvar <- rss/rdf
  } else {
    An <- Anova(object,type=3)
    effect.names <- rownames(An$anova)
    if(is.null(object$effect.sources))
      stop("Effect sources not found")
    errors <- An$errors
    err.df <- An$denom.df
    inds <- match(object$effect.sources[-1],effect.names)
    #  resvar <- c(rss/rdf,errors[inds]/err.df[inds])
    resvar <- c(rss/rdf,errors[inds])[Qr$pivot[p1]]
  }
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) * 
      1e-30) 
    warning("essentially perfect fit: summary may be unreliable")
  R    <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se   <- sqrt(diag(R) * resvar)
  est  <- z$coefficients[Qr$pivot[p1]]
  tval <- est/se
  ans  <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
                                                  rdf, lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans$aliased <- is.na(coef(object))
  ans$sigma <- sqrt(resvar[1])
  ans$df <- c(p, rdf, NCOL(Qr$qr))
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
                                                       df.int)/rdf[1])
    ans$fstatistic <- c(value = (mss/(p - df.int))/resvar[1], 
                        numdf = p - df.int, dendf = rdf[1])
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
                                                             1)]
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se)
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action)) 
    ans$na.action <- z$na.action
  if(!is.null(object$random) && !is.balanced(object)){
    cat("\nWARNING: Unbalanced data may lead to poor estimates\n")
  }
  class(ans) <- "summary.lmm"
  ans
}


###########################################
# confint.lm from stats, edited to enable limited classical least squares mixed models
confint.lmm <- function (object, parm, level = 0.95, ...) 
{
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qt(a, object$df.residual)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             pct))
  if(is.null(object$random)){
    ses <- sqrt(diag(vcov(object)))[parm]
  } else {
    p <- object$rank
    p1 <- 1L:p # Moved up for the sake of "random"
    Qr <- qr.lmm(object)
    rdf <- object$df.residual
    r <- object$residuals
    w <- object$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    } else {
      rss <- sum(w * r^2)
    }
    An <- Anova(object,type=3)
    effect.names <- rownames(An$anova)
    if(is.null(object$effect.sources))
      stop("Effect sources not found")
    errors <- An$errors
    err.df <- An$denom.df
    inds <- match(object$effect.sources[-1],effect.names)
    resvar <- c(rss/rdf,errors[inds])[Qr$pivot[p1]]
    R    <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    ses   <- sqrt(diag(R) * resvar)[match(parm,colnames(vcov(object)))]
  }
  
  ci[] <- cf[parm] + ses %o% fac
  ci
}
format.perc <- function (probs, digits) 
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), 
        "%")
