sandcov <- function(model, id) {
  model$weights <- pmax(1e-8, model$weights)
  l <- length(coef(model))
  vv <- summary(model)$cov.unscaled
  mm <- model.matrix(model)
  ii <- match(colnames(vv), colnames(mm))
  infl <- (residuals(model, "response")*model$prior.weights*mm[,ii])%*%vv
  jj <- match(colnames(mm), colnames(vv))
  infl <- infl[,jj]
  colnames(infl) <- colnames(mm)
  infl <- rowsum(infl,id)
  t(infl)%*%infl
}

haplo <- function(..., mode=c("additive", "dominant", "recessive"), group.rare=TRUE, rare.freq=0.02) {

  if(group.rare==FALSE) rare.freq <- 0

  rval <- as.matrix(cbind(...))

  mode <- match.arg(mode)

  attr(rval, "inheritance.mode") <- mode
  attr(rval, "group.rare") <- group.rare
  attr(rval, "rare.freq") <- rare.freq

  rval

}

count.haps <- function(h1, h2) {
  if(length(h1)!=length(h2)) stop("Error in haplotype length!")
  hh <- sort(union(h1,h2))
  table(1:length(h1), factor(h1, levels=hh))+table(1:length(h2), factor(h2, levels=hh))
}

return.haps <- function(mat) {
  mat <- as.matrix(mat)
  lis <- NULL
  for(i in 1:nrow(mat)){lis <- c(lis, as.numeric(paste(mat[i,], collapse="")))}
  return(lis)
}

one <- function(..., mode) {
  args <- list(...)
  rep(1, NROW(args[[1]]))
}

haplo.ccs <- function(formula, data=NULL, ...) {

  cl <- mf <- match.call(expand.dots=FALSE)

  if(missing(data))
    data <- environment(formula)

  mf$... <- NULL
  mf[[1]] <- as.name("model.frame")
  mf.new <- mf
  
  formula.new <- do.call("substitute", list(formula, list(haplo=as.name("one"))))
  formula.new <- eval(formula.new, parent.frame())
  mf.new$formula <- formula.new

  model.terms <- terms(formula, specials=c("haplo"))
  model.terms.new <- terms(formula.new, specials=c("one"))

  mf <- eval(mf, parent.frame())
  mf.new <- eval(mf.new, parent.frame())
  y <- model.response(mf, "numeric")

  one.main <- untangle.specials(model.terms.new, "one", order=1)$terms
  one.int <- untangle.specials(model.terms.new, "one", order=2:9)$terms

  mm <- model.matrix(model.terms.new, mf.new)
  assgn <- attr(mm, "assign")

  if(sum(is.na(match(colnames(mm), colnames(mf)))!=1)==0 & sum(one.int)==0)
    x <- NULL

  if(sum(is.na(match(colnames(mm), colnames(mf)))!=1)!=0 & sum(one.int)==0) {
    x <- mm[, assgn!=one.main]
    x <- x[, -1]
    names.x <- colnames(mm)[assgn!=one.main]
    names.x <- names.x[-1]
  }

  if(sum(is.na(match(colnames(mm), colnames(mf)))!=1)!=0 & sum(one.int)!=0) {
    test <- matrix(nrow=length(one.int), ncol=length(colnames(mm)))
    for(i in 1:length(one.int)) { test[i,] <- as.numeric(assgn==one.int[i]) }
    test <- rowSums(t(test))
    x <- mm[, (assgn!=one.main) & (test!=1)]
    x <- x[, -1]
    names.x <- colnames(mm)[(assgn!=one.main) & (test!=1)]
    names.x <- names.x[-1]
    int <- mm[, test==1]
    names.int <- gsub(":$","",gsub("^:","",gsub("(:|^)one\\(.*\\)(:|$)",":",colnames(mm)[test==1])))
  }

  inherit.mode <- attr(mf$haplo, "inheritance.mode")
  group.rare <- attr(mf$haplo, "group.rare")
  rare.freq <- attr(mf$haplo, "rare.freq")

  if(is.null(x)==1 & sum(one.int)==0)
    model <- haplo.ccs.fit(y=y, x=NULL, int=NULL, geno=mf$haplo, inherit.mode=inherit.mode, group.rare=group.rare, rare.freq=rare.freq, names.x=NULL, names.int=NULL, ...)

  if(is.null(x)==0 & sum(one.int)==0)
    model <- haplo.ccs.fit(y=y, x=data.frame(x), int=NULL, geno=mf$haplo, inherit.mode=inherit.mode, group.rare=group.rare, rare.freq=rare.freq, names.x=names.x, names.int=NULL, ...)

  if(is.null(x)==0 & sum(one.int)!=0)
    model <- haplo.ccs.fit(y=y, x=data.frame(x), int=data.frame(int), geno=mf$haplo, inherit.mode=inherit.mode, group.rare=group.rare, rare.freq=rare.freq, names.x=names.x, names.int=names.int, ...)

  class(model) <- c("haplo.ccs")

  model$formula <- formula
  model$call <- cl

  return(model)

}

haplo.ccs.fit <- function(y, x, int, geno, inherit.mode, group.rare, rare.freq, referent=return.haps(em$haplotype)[em$hap.prob==max(em$hap.prob)], names.x, names.int, ...) {

  em <- haplo.em(geno)

  y <- y[em$indx.sub]
  x <- x[em$indx.sub,]
  int <- int[em$indx.sub,]
  geno <- geno[em$indx.sub,]
  prob <- em$post
  id <- em$indx.sub

  rares <- (1:length(em$hap.prob))[em$hap.prob<=rare.freq]

  if((group.rare==TRUE)&(length(rares)!=0)) {
    rare.haps <- return.haps(em$haplotype[rares,])
    commons <- (1:length(em$hap.prob))[em$hap.prob>rare.freq]
    common.haps <- return.haps(em$haplotype[commons,])
    hap1test <- em$hap1code %in% rares
    hap2test <- em$hap2code %in% rares
    hap1code <- ifelse(hap1test==TRUE, 1+max(em$hap1code, em$hap2code), em$hap1code)
    hap2code <- ifelse(hap2test==TRUE, 1+max(em$hap1code, em$hap2code), em$hap2code)
    haplo.mat <- count.haps(hap1code, hap2code)
    names <- c(return.haps(em$haplotype)[commons], paste("Rare (<", rare.freq, ")", sep=""))
    colnames(haplo.mat) <- names
    haplo.mat <- haplo.mat[,c(colnames(haplo.mat)!=referent)==1]
  } else {
    haplo.mat <- count.haps(em$hap1code, em$hap2code)
    colnames(haplo.mat) <- paste(return.haps(em$haplotype))
    haplo.mat <- haplo.mat[,c(colnames(haplo.mat)!=referent)==1]
  }

  if(inherit.mode=="dominant")
    haplo.mat <- ifelse(haplo.mat==2, 1, haplo.mat)

  if(inherit.mode=="recessive")
    haplo.mat <- ifelse(haplo.mat==1, 0, ifelse(haplo.mat==2, 1, 0))
  
  if(is.null(x)==1 & is.null(int)==1)
    fit <- glm(y ~ haplo.mat, family=quasibinomial(link=logit), weights=prob, ...)

  if(is.null(x)==0 & is.null(int)==1) {
    x <- as.matrix(x)
    fit <- glm(y ~ haplo.mat + x, family=quasibinomial(link=logit), weights=prob, ...)
  }

  if(is.null(x)==0 & is.null(int)==0) {
    x <- as.matrix(x)
    int <- as.matrix(int)
    haplo.int <- NULL
    for(i in 1:dim(int)[2]){haplo.int <- cbind(haplo.int, apply(haplo.mat, 2, function(x){int[,i]*x}))}
    fit <- glm(y ~ haplo.mat + x + haplo.int, family=quasibinomial(link=logit), weights=prob, ...)
  }

  if((group.rare==TRUE)&(length(rares)!=0)) {
    fit$hap.probs <- as.matrix(c(em$hap.prob[c(return.haps(em$haplotype)==referent)==1],
                                 em$hap.prob[c((return.haps(em$haplotype)!=referent)&(return.haps(em$haplotype)%in%common.haps))==1],
                                 sum(em$hap.prob[c((return.haps(em$haplotype)!=referent)&(return.haps(em$haplotype)%in%rare.haps))==1])))
    fit$hap.names <- c(paste(referent, "(Ref)"), colnames(haplo.mat))
    rownames(fit$hap.probs) <- c(fit$hap.names)
    colnames(fit$hap.probs) <- c("Frequency")
  } else {
    fit$hap.probs <- as.matrix(c(em$hap.prob[c(return.haps(em$haplotype)==referent)==1],
                                 em$hap.prob[c(return.haps(em$haplotype)!=referent)==1]))
    fit$hap.names <- c(paste(referent, "(Ref)"), colnames(haplo.mat))
    rownames(fit$hap.probs) <- c(fit$hap.names)
    colnames(fit$hap.probs) <- c("Frequency")
  }

  fit$coefficients <- coef(fit)
  fit$covariance <- sandcov(fit, id)

  if(is.null(x)==1 & is.null(int)==1) {
    names(fit$coefficients) <- colnames(fit$covariance) <- rownames(fit$covariance)<- fit$hap.names
  }

  if(is.null(x)==0 & is.null(int)==1) {
    names(fit$coefficients) <- colnames(fit$covariance) <- rownames(fit$covariance)<- c(fit$hap.names, names.x)
  }

  if(is.null(x)==0 & is.null(int)==0) {
    int.names <- NULL
    for(i in 1:length(names.int)){int.names <- c(int.names, paste(fit$hap.names[-1], ":", names.int[i], sep=""))}
    names(fit$coefficients) <- colnames(fit$covariance) <- rownames(fit$covariance)<- c(fit$hap.names, names.x, int.names)
  }

  names(fit$residuals) <- names(fit$fitted.values) <- names(fit$linear.predictors) <- names(fit$weights) <- names(y) <- id

  if(is.null(fit$offset)!=1)
    names(fit$offset) <- id

  fit$df <- fit$df.residual

  fit$inheritance.mode <- inherit.mode
  fit$rare.freq <- rare.freq

  fit$em.lnlike <- em$lnlike
  fit$em.lr <- em$lr
  fit$em.df.lr <- em$df.lr
  fit$em.converged <- as.logical(em$converge)
  fit$prior.weights <- em$post
  fit$hap1 <- return.haps(em$haplotype)[em$hap1code]
  fit$hap2 <- return.haps(em$haplotype)[em$hap2code]
  
  fit$em.nreps <- em$nreps
  fit$em.max.pairs <- em$max.pairs
  fit$em.control <- em$control

  names(fit$hap1) <- names(fit$hap2) <- names(fit$prior.weights) <- id
  fit$id <- id

  c(fit[c("coefficients", "covariance", "residuals", "fitted.values", "linear.predictors", "df", 
          "rank", "family", "iter", "weights", "prior.weights", "y", "id", "converged", "boundary", 
          "model", "terms", "offset", "contrasts", "xlevels", "inheritance.mode", "rare.freq",
          "em.lnlike", "em.lr", "em.df.lr", "hap1", "hap2", "hap.names", "hap.probs", 
          "em.converged", "em.nreps", "em.max.pairs", "em.control")])

}

print.haplo.ccs <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nFormula:", deparse(x$formula), "\n\n")
  cat("Relative Risks:", "\n") 
  print(round(exp(x$coefficients), digits))
  cat("\n")
  cat("Haplotypes:", "\n")
  print(round(x$hap.probs, digits))
  cat("\n")
}

summary.haplo.ccs <- function(object, ...) {
  ans <- c(object[c("call", "formula", "coefficients", "covariance", "terms", "contrasts", 
                   "iter", "df", "hap.names", "hap.probs")])

  class(ans) <- c("summary.haplo.ccs")
  return(ans)
}

print.summary.haplo.ccs <- function(x, digits=max(3, getOption("digits")-3), ...) {
  coef <- x$coefficients
  se <- sqrt(diag(x$covariance))
  tstat <- coef/se
  p <- 2 * pt(-abs(tstat), x$df)

  summary <- cbind(exp(coef), se, tstat, p)
  colnames(summary) <- c("Relative Risk", "Robust SE", "t Value", "P(T>|t|)")
  rownames(summary) <- names(coef)

  cat("\nFormula:", deparse(x$formula), "\n\n")
  cat("Estimates:", "\n") 
  print(round(summary, digits))
  cat("\n")
  cat("Haplotypes:", "\n")
  print(round(x$hap.probs, digits))
  cat("\nNumber of Fisher Scoring Iterations:", x$iter, "\n\n")
}

coef.haplo.ccs <- function(object, ...) {
  return(object$coefficients)
}

vcov.haplo.ccs <- function(object, ...) {
  return(object$covariance)
}

residuals.haplo.ccs <- function(object, ...) {
  return(object$residuals)
}

fitted.haplo.ccs <- function(object, ...) {
  return(object$fitted.values)
}

haplo.freq <- function(object, ...) {
  return(object$hap.probs)
}

anova.haplo.ccs <- function(object, ...) {
  stop("haplo.ccs does not fit by maximum likelihood.")
}

AIC.haplo.ccs <- function(object, ..., k) {
  stop("haplo.ccs does not fit by maximum likelihood.")
}

logLik.haplo.ccs <- function(object, ...) {
  stop("haplo.ccs does not fit by maximum likelihood.")
}

haplo.hpp <- function(model, prob=0.95) {

  prob <- as.vector(prob)

  hpp <- vector(length=length(unique(model$id)))
  for(i in unique(model$id)) {
    hpp[i] <- max(model$prior.weights[model$id==i])
  }

  out <- matrix(nrow=length(hpp), ncol=length(prob))
  for(i in 1:length(prob)) {
    for(j in 1:length(hpp)) {
      out[j,i] <- hpp[j]>=prob[i]
    }
  }

  num <- apply(out, 2, sum)
  total <- length(unique(model$id))
  summary <- t(as.matrix(num/total))
  colnames(summary) <- prob
  rownames(summary) <- "Proportion"

  return(round(summary, max(3, getOption("digits")-3)))

}
