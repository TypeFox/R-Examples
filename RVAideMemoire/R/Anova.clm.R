# Fonctions ou méthodes redéfinies car non exportées par car :
#   - df.terms.clm_m
#   - relatives
#   - responseName.clm_m

df.terms.clm_m <- function(model,term,...){
  if (!missing(term) && 1==length(term)) {
    assign <- attr(model.matrix(model)$X,"assign")
    which.term <- which(term==labels(terms(model)))
    if (0==length(which.term)) {stop(paste(term,"is not in the model."))}
    sum(assign==which.term)
  } else {
    terms <- if (missing(term)) {
	labels(terms(model))
    } else {
	term
    }
    result <- numeric(0)
    for (term in terms) {
	result <- c(result,Recall(model,term))
    }
    names(result) <- terms
    result
  }
}

relatives <- function(term,names,factors){
  is.relative <- function(term1,term2) {
    all(!(factors[,term1]&(!factors[,term2])))
  }
  if(length(names)==1) {return(NULL)}
  which.term <- which(term==names)
  (1:length(names))[-which.term][sapply(names[-which.term],
    function(term2) is.relative(term,term2))]
}

responseName.clm_m <- function (model,...) {
  deparse(attr(terms(model),"variables")[[2]])
}

Anova.clm <- function(mod,type=c("II","III",2,3),...) {
  type <- as.character(type)
  type <- match.arg(type)
  if (is.null(mod$beta) && (type=="2" || type=="II")) {
    type <- "III"
    warning("the model contains only intercepts: Type III test substituted")
  }
  switch(type,II=Anova.II.clm(mod,...),III=Anova.III.clm(mod,...),
    `2`=Anova.II.clm(mod,...),`3`=Anova.III.clm(mod,...))
}

Anova.II.clm <- function(mod,...) {
  fac <- attr(mod$terms,"factors")
  names <- labels(terms(mod))
  n.terms <- length(names)
  p <- LR <- rep(0,n.terms)
  df <- df.terms.clm_m(mod)
  for (term in 1:n.terms) {
    rels <- names[relatives(names[term],names,fac)]
    exclude.1 <- if (length(rels)==0) {
	names[term]
    } else {
	paste(names[term],paste(rels,collapse="+"),sep="+")
    }
    mod.1 <- update(mod,as.formula(paste0(".~.-(",exclude.1,")")))
    dev.1 <- -2*logLik(mod.1)
    mod.2 <- if (length(rels)==0) { 
	mod
    } else {
	exclude.2 <- paste(rels,collapse="+")
	update(mod,as.formula(paste0(".~.-(",exclude.2,")")))
    }
    dev.2 <- -2*logLik(mod.2)
    LR[term] <- dev.1-dev.2
    p[term] <- pchisq(LR[term],df[term],lower.tail=FALSE)
  }
  result <- data.frame(LR,df,p)
  row.names(result) <- names
  names(result) <- c("LR Chisq","Df","Pr(>Chisq)")
  class(result) <- c("anova","data.frame")
  attr(result,"heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
    paste("Response:",responseName.clm_m(mod)))
  result
}

Anova.III.clm <- function(mod,...) {
  names <- labels(terms(mod))
  n.terms <- length(names)
  p <- LR <- rep(0,n.terms)
  df <- df.terms.clm_m(mod)
  deviance <- -2*logLik(mod)
  for (term in 1:n.terms) {
    mod.1 <- update(mod,as.formula(paste0(".~.-(",names[term],")")))
    LR[term] <- -2*logLik(mod.1)-deviance
    p[term] <- pchisq(LR[term],df[term],lower.tail=FALSE)
  }
  result <- data.frame(LR,df,p)
  row.names(result) <- names
  names(result) <- c("LR Chisq","Df","Pr(>Chisq)")
  class(result) <- c("anova","data.frame")
  attr(result,"heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
    paste("Response:",responseName.clm_m(mod)))
  result
}
