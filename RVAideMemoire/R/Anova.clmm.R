# La méthode model.matrix.clmm est identique à model.matrix.clm
#   (qui est dans ordinal), mais il n'y a pas de méthode prévue
#   pour les clmm
model.matrix.clmm <- function(object,type=c("design","B"),...) {
  type <- match.arg(type)
  if (type == "design") {
    mf <- update(object,method="model.frame")
    keep <- c("X","NOM","S")
    select <- match(keep,names(mf),nomatch=0)
    return(mf[select])
  } else {
    env <- update(object,doFit=FALSE)
    ans <- list(B1=env$B1,B2=env$B2)
    ans$S <- env$S
    return(ans)
  }
}

Anova.clmm <- function(mod,type=c("II","III",2,3),...) {
  type <- as.character(type)
  type <- match.arg(type)
  if (is.null(mod$beta) && (type=="2" || type=="II")) {
    type <- "III"
    warning("the model contains only intercepts: Type III test substituted")
  }
  switch(type,II=Anova.II.clmm(mod,...),III=Anova.III.clmm(mod,...),
    `2`=Anova.II.clmm(mod,...),`3`=Anova.III.clmm(mod,...))
}

Anova.II.clmm <- function(mod,...) {
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
    mod.1 <- if (length(names)-(length(rels)+1)>0) {
	update(mod,as.formula(paste0(".~.-(",exclude.1,")")))
    } else {
	null.form <- update.formula(formula(mod),as.formula(paste0(".~.-(",exclude.1,")")))
	null.form2 <- paste0(as.character(null.form)[2],"~1+",as.character(null.form)[3])
	call <- mod$call
	call$formula <- as.formula(null.form2)
	eval(call)
    }
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

Anova.III.clmm <- function(mod,...) {
  names <- labels(terms(mod))
  n.terms <- length(names)
  p <- LR <- rep(0,n.terms)
  df <- df.terms.clm_m(mod)
  deviance <- -2*logLik(mod)
  for (term in 1:n.terms) {
    mod.1 <- if (n.terms>1) {
	update(mod,as.formula(paste0(".~.-(",names[term],")")))
    } else {
	null.form <- update.formula(formula(mod),as.formula(paste0(".~.-(",names[term],")")))
	null.form2 <- paste0(as.character(null.form)[2],"~1+",as.character(null.form)[3])
	call <- mod$call
	call$formula <- as.formula(null.form2)
	eval(call)
    }
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

