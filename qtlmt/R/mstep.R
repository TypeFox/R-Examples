
####################################
#  functions for model selection:  #
# multivariate multiple regression #
################################################################################

addscope<- function(object,scope,...){
# scope: upper model or terms to add
  terms1<- attr(terms(object),"term.labels")
  if(!is.character(scope)){
    terms2<- attr(terms(update.formula(object, scope,...)),"term.labels")
    if(any(is.na(match(terms1,terms2)),na.rm = FALSE))
      stop("upper scope does not include model...")
    setdiff(terms2,terms1)
  }else{
    setdiff(scope,terms1)
  }
}

dropscope<- function(object,scope,...){
# scope: lower model or terms to drop
  terms1<- attr(terms(object),"term.labels")
  if(!is.character(scope)){
    terms2<- attr(terms(update.formula(object, scope,...)),"term.labels")
    if(any(is.na(match(terms2,terms1)),na.rm = FALSE))
      stop("lower scope does not include model...")
    fdrop<- setdiff(terms1,terms2)
  }else{
    drop.rhs <- paste(scope, collapse = "+")
    drop.rhs <- eval(parse(text = paste("~ . +", drop.rhs)))
    terms2 <- attr(terms(update.formula(object, drop.rhs,...)),"term.labels")
    if(any(is.na(match(scope,terms1)),na.rm = FALSE))
      stop("scope is not in model...")
    fdrop<- scope
  }
  fdrop
}

fscope<- function(scope){
# scope: vector of predictors
  if(length(scope)){
    formula(paste("~",paste(scope,collapse="+",sep="")))
  }else{
    formula(paste("~",1))
  }
}

mlogLik<- function(object){
# object: regression object
  rs<- as.matrix(object$residuals)
  n<- dim(rs)[1]
  p<- dim(rs)[2]
  E<- crossprod(rs)/n
  detE<- det(E)
  loglik<- -n*p/2*(log(2*pi)+1) - n/2*log(detE)
  
  loglik
}

#################################
###  multivariate regression  ###
### add single terms to model ###
#################################
mAdd1 <-
   function(object,
            scope,
            test = c("none", "Chisq", "F"),
            k = 0,
            ...)
{
   UseMethod("mAdd1")
}
mAdd1.default <- function(object, scope, test=c("none", "Chisq", "F"), k=0, ...)
{#scope: upper model that contains the base or terms to add
  Fstat <- function(WL,p,vH,vE) { #WL: Wilks' Lamda
    w<- vE + vH - (p+vH+1)/2
    t<- rep(1,length(w))
    choice<- p^2+vH^2-5 > 0 & vH > 0
    choice<- !is.na(choice) & choice
    if(sum(choice)>0){
      t[choice]<- sqrt((p^2*vH^2-4)[choice]/(p^2+vH^2-5)[choice])
    }
    df1<- p*vH
    df2<- w*t-(df1-2)/2
    Fs<- (1-WL^(1/t))/(WL^(1/t))*df2/df1
    Fs[df1 < .Machine$double.xmin] <- NA
    P <- Fs
    nnas <- !is.na(Fs)
    P[nnas] <- pf(Fs[nnas], df1[nnas], df2[nnas], lower.tail=FALSE)
    list(Fs=Fs, P=P)
  }
  Chistat <- function(WL,p,vH,vE) { #WL: Wilks' Lamda
    f<- vE-(p-vH+1)/2
    Chi2<- -f*log(WL)
    df<- p*vH
    Chi2[df < .Machine$double.xmin] <- NA
    P <- Chi2
    nnas <- !is.na(Chi2)
    P[nnas] <- pchisq(Chi2[nnas], df[nnas], lower.tail=FALSE)
    list(Chi2=Chi2, P=P)
  }
  Edet<- function(object){
    rs<- object$residuals
    nobs<- object$rank+object$df.residual
    E<- crossprod(rs)/nobs
    det(E)
  }

  if(missing(scope) || is.null(scope)) stop("no terms in scope...")
  scope<- addscope(object,scope)
  if(!length(scope))
    stop("no terms in scope for adding to object...")
  oTerms <- attr(object$terms, "term.labels")
  int <- attr(object$terms, "intercept")
  ns <- length(scope)
  y <- object$residuals + predict(object)
    y<- as.matrix(y)
  n<- dim(y)[1]
  p<- dim(y)[2]
  vH<- numeric(ns+1)
  Es <- numeric(ns+1)
  names(vH) <- names(Es) <- c("<none>", scope)
  vE<- vH
  AIC<- vH
  vH[1] <- object$rank
  vE[1]<- object$df.residual
  Es[1] <- Edet(object)
  const<- -n*p/2*(log(2*pi)+1)
  if(Es[1] > 0){
     lik<- const - n/2*log(Es[1])
     AIC[1]<- -2*lik + k*(vH[1]*p + p*(p+1)/2)
  }else AIC[1]<- Inf
  for(tt in scope) {
    fit<- NULL
    fit <- update(object, paste("~ . + ", tt), evaluate = TRUE,...)
    vH[tt] <- fit$rank
    vE[tt]<- fit$df.residual
    Es[tt] <- Edet(fit)
    if(Es[tt] > 0){
       lik<- const - n/2*log(Es[tt])
       AIC[tt]<- -2*lik + k*(vH[tt]*p + p*(p+1)/2)
    }else AIC[tt]<- Inf
  }
  vH[-1] <- vH[-1] - vH[1]; vH[1]<- NA
  WLs<- c(NA, Es[-1]/Es[1])
  aod <- data.frame(Df = vH, "Wilks' Lambda" = WLs, AIC=AIC,
      row.names = names(vH), check.names = FALSE)
  test <- match.arg(test)
  if(test == "Chisq") {
    aod[, c("Modified Chisq", "Pr(Chi)")] <- Chistat(WLs, p, vH, vE)
  } else if(test == "F") {
    aod[, c("F value", "Pr(F)")] <- Fstat(WLs, p, vH, vE)
  }
  head <- c("Single term additions", "\nModel:",
    deparse(as.vector(formula(object))))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

####################################
###  multivariate regression   ###
### drop single terms from model ###
####################################
mDrop1 <-
   function(object,
            scope,
            test = c("none", "Chisq", "F"),
            k = 0,
            ...)
{
   UseMethod("mDrop1")
}

mDrop1.default <- function(object, scope, test=c("none", "Chisq", "F"), k=0, ...)
{#scope: lower model that is contained in the base or terms to drop
  Fstat <- function(WL,p,vH,vE) { #WL: Wilks' Lamda
    w<- vE + vH - (p+vH+1)/2
    t<- rep(1,length(w))
    choice<- p^2+vH^2-5 > 0 & vH > 0
    choice<- !is.na(choice) & choice
    if(sum(choice)>0){
      t[choice]<- sqrt((p^2*vH^2-4)[choice]/(p^2+vH^2-5)[choice])
    }
    df1<- p*vH
    df2<- w*t-(df1-2)/2
    Fs<- (1-WL^(1/t))/(WL^(1/t))*df2/df1
    Fs[df1 < .Machine$double.xmin] <- NA
    P <- Fs
    nnas <- !is.na(Fs)
    P[nnas] <- pf(Fs[nnas], df1[nnas], df2[nnas], lower.tail=FALSE)
    list(Fs=Fs, P=P)
  }
  Chistat <- function(WL,p,vH,vE) { #WL: Wilks' Lamda
    f<- vE-(p-vH+1)/2
    Chi2<- -f*log(WL)
    df<- p*vH
    Chi2[df < .Machine$double.xmin] <- NA
    P <- Chi2
    nnas <- !is.na(Chi2)
    P[nnas] <- pchisq(Chi2[nnas], df[nnas], lower.tail=FALSE)
    list(Chi2=Chi2, P=P)
  }
  Edet<- function(object){
    rs<- object$residuals
    nobs<- object$rank+object$df.residual
    E<- crossprod(rs)/nobs
    det(E)
  }

  x <- model.matrix(object)
  iswt <- !is.null(wt <- object$weights)
  n <- nrow(x)
  asgn <- attr(x, "assign")
  tl <- attr(object$terms, "term.labels")
  if(missing(scope)){
    scope <- attr(terms(object),"term.labels")
  }else {
    scope<- dropscope(object,scope)
    if(!all(match(scope, tl, FALSE)))
      stop("scope is not a subset of term labels")
  }
  ndrop <- match(scope, tl)
  ns <- length(scope)
  y <- object$residuals + predict(object)
    y<- as.matrix(y)
  p<- dim(y)[2]
  vH<- numeric(ns)
  Es <- numeric(ns)
  AIC<- vH

  const<- -n*p/2*(log(2*pi)+1)
  for(i in 1:ns) {
    ii <- seq(along=asgn)[asgn == ndrop[i]]
    jj <- setdiff(seq(ncol(x)), ii)
    z <- if(iswt) lm.wfit(x[, jj, drop = FALSE], y, wt)
    else lm.fit(x[, jj, drop = FALSE], y)
    vH[i] <- z$rank
    Es[i] <- Edet(z)
    if(Es[i] > 0){
       lik<- const - n/2*log(Es[i])
       AIC[i]<- -2*lik + k*(vH[i]*p + p*(p+1)/2)
    }else AIC[i]<- Inf
  }
  vH<- object$rank - vH; vH<- c(NA,vH)
  vE<- object$df.residual
  detObj<- Edet(object)
  WLs<- detObj/Es; WLs<- c(NA,WLs)
  if(detObj > 0){
     lik<- const - n/2*log(detObj)
     AIC<- c(-2*lik + k*(object$rank*p + p*(p+1)/2),AIC)
  }else AIC<- c(Inf, AIC)
  aod <- data.frame(Df = vH, "Wilks' Lambda" = WLs, AIC=AIC,
        row.names =c("<none>" ,scope), check.names = FALSE)
  test <- match.arg(test)
  if(test == "Chisq") {
    aod[, c("Modified Chisq", "Pr(Chi)")] <- Chistat(WLs, p, vH, vE)
  } else if(test == "F") {
    aod[, c("F value", "Pr(F)")] <- Fstat(WLs, p, vH, vE)
  }
  head <- c("Single term deletions", "\nModel:",
      deparse(as.vector(formula(object))))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

###################################
###   multivariate regression   ###
###  stepwise model selection   ###
###################################
mStep<-
   function(object,
            scope,
            direction = c("both", "backward", "forward"),
            trace = FALSE,
            keep = TRUE,
            steps = 1000,
            k = 2,
            ...)
{
   UseMethod("mStep")
}

mStep.default<- function (object, scope, direction = c("both", "backward", "forward"),
  trace = FALSE, keep = TRUE, steps = 1000, k = 2, ...){
  y <- object$residuals + predict(object)
    y<- as.matrix(y)
  n<- dim(y)[1]
  p<- dim(y)[2]
  Edet<- function(object){
    rs<- object$residuals
    nobs<- object$rank+object$df.residual
    E<- crossprod(rs)/nobs
    det(E)
  }
  myaic<- function(object,k,p){
    const<- -n*p/2*(log(2*pi)+1)
    lik<- const - n/2*log(Edet(object))
    -2*lik + k*(object$rank*p + p*(p+1)/2)
  }
  cut.string <- function(string) {
    if (length(string) > 1) 
      string[-1] <- paste("\n", string[-1], sep = "")
    string
  }
  Terms <- terms(object)
  object$call$formula <- object$formula <- Terms
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric(0)
    fadd <- attr(Terms, "factors")
    if (md) 
      forward <- FALSE
  }else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) 
         attr(terms(update.formula(object, fdrop)), "factors")
      else numeric(0)
      fadd <- if (!is.null(fadd <- scope$upper)) 
         attr(terms(update.formula(object, fadd)), "factors")
    }else {
      fadd <- if (!is.null(fadd <- scope)) 
        attr(terms(update.formula(object, scope)), "factors")
      fdrop <- numeric(0)
    }
  }
  fit <- object
  fit0<- fit
  bAIC <- myaic(fit,k,p)
  edf <- fit$rank
  nm <- 1
  Terms <- fit$terms
  if (trace) 
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n", 
    cut.string(deparse(as.vector(formula(fit)))), "\n\n")
  out<- NULL
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && length(scope$drop)) {
      aod <- mDrop1(fit, scope$drop, test="none", k = k, ...)
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1], paste("-", rn[-1], sep = " "))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        change <- rev(rownames(aod)[zdf])[1]
      }
    }
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- mAdd1(fit, scope$add, test="none", k = k, ...)
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1], paste("+", rn[-1], sep = " "))
        aod <- if (is.null(aod)) 
          aodf
        else rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df)) 
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nc <- match("AIC", names(aod))
      nc <- nc[!is.na(nc)][1]
      o <- order(aod[, nc])
      if (trace) 
        print(aod[o, ])
      if (o[1] == 1 || !is.finite(aod[o[1], nc]))
        break
      change <- rownames(aod)[o[1]]
    }
    fit0 <- update(fit, paste("~ .", change), evaluate = TRUE,...)
#    fit0 <- eval.parent(fit0)
    if (dim(as.matrix(fit0$residuals))[1] != n)
      stop("number of rows in use has changed: remove missing values?")
    bAIC <- myaic(fit0,k,p)
    if (bAIC >= AIC + 1e-7) 
      break
    fit<- fit0
    out<- rbind(out,aod[change,])
    Terms <- terms(fit)
    edf <- fit$rank
    if (trace) 
      cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n", 
        cut.string(deparse(as.vector(formula(fit)))), "\n\n")
    nm <- nm + 1
  }
  if(keep) fit$keep<- out
  fit
}

################################################################################
# the end #
###########

if(FALSE){ # no longer in need

##################################
###  multivariate regression   ###
###  backward model selection  ###
##################################
mbackwd<- function (object, scope, test=c("none","Chisq","F"), 
  k=Inf, keep=TRUE, trace = FALSE, steps = 1000, ...){
# multivariate regression model selection
# scope: lower model that is contained in base
  cut.string <- function(string) {
    if (length(string) > 1) 
      string[-1] <- paste("\n", string[-1], sep = "")
    string
  }
  Edet<- function(object){
    rs<- object$residuals
    nobs<- object$rank+object$df.residual
    E<- crossprod(rs)/nobs
    det(E)
  }
  Terms <- terms(object)
  object$call$formula <- object$formula <- Terms
  test <- match.arg(test)
  mt <- missing(test) || test=="none"
  n <- dim(as.matrix(object$residuals))[1]
  p<- dim(as.matrix(object$residuals))[2]
  fit <- object
  fdrop<- dropscope(fit,scope)
  nm <- 1
  Terms <- fit$terms
  if (trace) 
    cat("Start:", "\n", cut.string(deparse(as.vector(formula(fit)))), "\n\n")

  const<- -n*p/2*(log(2*pi)+1)
  lik<- const - n/2*log(Edet(fit))
  bAIC<- -2*lik + k*(fit$rank*p + p*(p+1)/2)
  out<- NULL
  while (steps > 0) {
    steps <- steps - 1
    aod <- NULL
    change <- NULL
    if (length(fdrop)) {
      aod <- mDrop1(fit, fdrop, test=test, k=k, ...)
    }
    if (is.null(aod) || ncol(aod) == 1) break
    nc <- match(c("Wilks' Lambda"), colnames(aod))
    nc <- nc[!is.na(nc)][1]
    o <- order(aod[, nc],decreasing=T)
    aod<- aod[o,]
    change <- rownames(aod)[1]
    fdrop<- setdiff(fdrop,change)
    
    if(aod$AIC[1]>bAIC+1e-7) break else bAIC<- aod$AIC[1]
    out<- rbind(out,aod[1,])
    if (trace){
      rn <- rownames(aod)
      rownames(aod) <- paste("-", rn, sep = " ")
      print(aod)
    }
    fit <- update(fit, paste("~ . -", change), evaluate = TRUE,...)
#    fit <- eval.parent(fit)
    if (dim(as.matrix(fit$residuals))[1] != n) 
      stop("number of rows in use has changed: remove missing values?")
    if (trace) 
      cat("\nStep", nm, "\n", cut.string(deparse(as.vector(formula(fit)))),"\n\n")
    nm <- nm + 1
  }

  attr(out,"heading")<- NULL
  if(keep) fit$keep<- out
  fit
}

##################################
###  multivariate regression   ###
###  forward model selection   ###
##################################
mforwd<- function (object, scope, test=c("none","Chisq","F"),
  k=0, keep=TRUE, trace = FALSE, steps = 1000, ...){
# multivariate regression forward selection
# scope: upper model that contains base
# maxr: max number of terms to select
  cut.string <- function(string) {
    if (length(string) > 1) 
      string[-1] <- paste("\n", string[-1], sep = "")
    string
  }
  Edet<- function(object){
    rs<- object$residuals
    nobs<- object$rank+object$df.residual
    E<- crossprod(rs)/nobs
    det(E)
  }
  Terms <- terms(object)
  object$call$formula <- object$formula <- Terms
  test <- match.arg(test)
  mt <- missing(test) || test=="none"
  n <- dim(as.matrix(object$residuals))[1]
  p<- dim(as.matrix(object$residuals))[2]
  fit <- object
  fadd <- addscope(fit, scope)
  nm <- 1
  Terms <- fit$terms
  if (trace) 
    cat("Start:", "\n", cut.string(deparse(as.vector(formula(fit)))), "\n\n")

  const<- -n*p/2*(log(2*pi)+1)
  lik<- const - n/2*log(Edet(fit))
  bAIC<- -2*lik + k*(fit$rank*p + p*(p+1)/2)
  out<- NULL
  while (steps > 0) {
    steps <- steps - 1
    aod <- NULL
    change <- NULL
    if (length(fadd)) {
      aod <- mAdd1(fit, fadd, test=test, k=k, ...)
    }
    if (is.null(aod) || ncol(aod) == 1) break
    nc <- match("Wilks' Lambda", names(aod))
    nc <- nc[!is.na(nc)][1]
    o <- order(aod[, nc],decreasing=F)
    aod<- aod[o,]
    change <- rownames(aod)[1]
    fadd<- setdiff(fadd,change)

    if(aod$AIC[1]>bAIC+1e-7) break else bAIC<- aod$AIC[1]
    out<- rbind(out,aod[1,])
    if (trace){
      rn <- rownames(aod)
      rownames(aod) <- paste("+", rn, sep = " ")
      print(aod)
    }
    fit <- update(fit, paste("~ . + ", change), evaluate = TRUE,...)
#    fit <- eval.parent(fit)
    if (dim(fit$residuals)[1] != n) 
      stop("number of rows in use has changed: remove missing values?")
    if (trace) 
      cat("\nStep", nm, "\n", cut.string(deparse(as.vector(formula(fit)))),"\n\n")
    nm <- nm + 1
  }

  attr(out,"heading")<- NULL
  if(keep) fit$keep<- out
  fit
}

##############################################################################
# functions for model selection: multivariate multiple regression -- extensive
##############################################################################

mMove<- function(obj,scope,...){
# scope: all candidate predictors
  ob<- obj
  ins<- attr(terms(ob),"term.labels")
  if(any(is.na(match(ins,scope)),na.rm = FALSE))
    stop("scope does not include model...")
  yes<- TRUE
  nn<- length(ins); if(nn<1) yes<- FALSE
  outs<- setdiff(scope,ins); if(length(outs)<1) yes<- FALSE
  lik<- mlogLik(ob)
  while(yes){
    yes<- FALSE
    for(n in 1:nn){
      ins1<- ins
      for(x in outs){
        ins1[n]<- x; scp1<- sort(ins1)
        scp1<- fscope(ins1)
        ob1<- update(ob, scp1, evaluate = TRUE,...)
        lik1<- mlogLik(ob1)
        if(lik1>lik+1e-8){
          ins<- ins1
          lik<- lik1
          ob<- ob1
          yes<- TRUE
        }
      }
      outs<- setdiff(scope,ins)
    }
  }
  ob
}

#################################
###  multivariate regression  ###
### add single terms to model ###
#################################
mAdd1 <-
   function(obj,
            scope,
            ext = FALSE,
            ...)
{
   UseMethod("mAdd1")
}

mAdd1.default <- function(obj, scope, ext=FALSE, ...){
# scope: upper model that contains the base or terms to add
  if(missing(scope) || is.null(scope)) stop("no terms in scope...")
  ins<- attr(terms(obj),"term.labels")
  if(is.null(scope)||missing(scope)){
    scope<- list()
    scope$lower<- numeric(0)
    scope$upper<- ins
  }else if(is.list(scope)){
    if(length(scope)==1){
      scope0<- scope[[1]]
      scope<- list()
      scope$lower<- ins
      scope$upper<- scope0
    }else if(length(scope)==2){
      flg1<- any(is.na(match(scope[[1]],scope[[2]])),na.rm = FALSE)
      flg2<- any(is.na(match(scope[[2]],scope[[1]])),na.rm = FALSE)
      if(flg1&&flg2){
        stop("scope: wrong...")
      }else if(flg1){
        scope$lower<- scope[[2]]
        scope$upper<- scope[[1]]
      }else{
        scope$lower<- scope[[1]]
        scope$upper<- scope[[2]]
      }
    }else stop("scope: wrong...")
  }else{
    scope0<- scope
    scope<- list()
    scope$lower<- ins
    scope$upper<- scope0
  }

  if(any(is.na(match(ins,scope$upper)),na.rm = FALSE))
    stop("scope$upper does not include model...")

  yes<- TRUE
  outs<- setdiff(scope$upper,ins)
  if(!length(outs)){
#    cat("no terms in scope for adding...")
    yes<- FALSE
  }

  if(ext) ob<- mMove(obj,scope$upper,...) else ob<- obj
  add<- FALSE
  if(yes){
    yes<- FALSE
    lik0<- mlogLik(ob)
    ob0<- ob
    for(tt in outs) {
      ob1<- NULL
      ob1 <- update(ob, paste("~ . + ", tt), evaluate = TRUE,...)
      lik1<- mlogLik(ob1)
      if(lik1>lik0+1e-8){
        ob0<- ob1
        lik0<- lik1
        add<- TRUE
        yes<- TRUE
      }
    }
    if(yes){
      ob<- ob0
      if(ext) ob<- mMove(ob,scope$upper,...)
    }
  }

  list(fit=ob,add=add)
}

####################################
###  multivariate regression   ###
### drop single terms from model ###
####################################
mDrop1 <-
   function(obj,
            scope,
            ext = FALSE,
            ...)
{
   UseMethod("mDrop1")
}

mDrop1.default <- function(obj, scope, ext=FALSE, ...){
# scope: lower model that is contained in the base or terms to drop
  if(missing(scope) || is.null(scope)) stop("no terms in scope...")
  ins<- attr(terms(obj),"term.labels")
  if(is.null(scope)||missing(scope)){
    scope<- list()
    scope$lower<- numeric(0)
    scope$upper<- ins
  }else if(is.list(scope)){
    if(length(scope)==1){
      scope0<- scope[[1]]
      scope$lower<- scope0
      scope$upper<- ins
    }else if(length(scope)==2){
      flg1<- any(is.na(match(scope[[1]],scope[[2]])),na.rm = FALSE)
      flg2<- any(is.na(match(scope[[2]],scope[[1]])),na.rm = FALSE)
      if(flg1&&flg2){
        stop("scope: wrong...")
      }else if(flg1){
        scope$lower<- scope[[2]]
        scope$upper<- scope[[1]]
      }else{
        scope$lower<- scope[[1]]
        scope$upper<- scope[[2]]
      }
    }else stop("scope: wrong...")
  }else{
    scope0<- scope
    scope<- list()
    scope$lower<- scope0
    scope$upper<- ins
  }

  if(any(is.na(match(scope$lower,ins)),na.rm = FALSE))
    stop("scope$lower is not a subset of term labels")

  yes<- TRUE
  outs<- setdiff(ins,scope$lower)
  if(!length(outs)){
#    cat("no terms in scope for dropping...")
    yes<- FALSE
  }

  if(ext) ob<- mMove(obj,scope$upper,...) else ob<- obj
  drop<- FALSE
  if(yes){
    yes<- FALSE
    lik0<- -Inf
    ob0<- ob
    for(n in 1:length(outs)) {
      scp1<- setdiff(ins,outs[n])
      scp1<- fscope(scp1)
      ob1<- update(ob, scp1, evaluate = TRUE,...)
      lik1<- mlogLik(ob1)
      if(lik1>lik0+1e-8){
        ob0<- ob1
        lik0<- lik1
        drop<- TRUE
        yes<- TRUE
      }
    }
    if(yes){
      ob<- ob0
      if(ext) ob<- mMove(ob,scope$upper,...)
    }
  }

  list(fit=ob,drop=drop)
}

###################################
###   multivariate regression   ###
###  stepwise model selection   ###
###################################
mStep<-
   function (obj,
             scope,
             cv,
             direction = c("both", "backward", "forward"),
             steps = 1000,
             ext = FALSE,
             ...)
{
   UseMethod("mStep")
}

mStep.default<- function (obj, scope, cv, direction=c("both","backward","forward"),
  steps = 1000, ext=FALSE, ...){
  direction<- match.arg(direction)
  ins<- attr(terms(obj),"term.labels")
  if(is.null(scope)||missing(scope)){
    scope<- list()
    scope$lower<- numeric(0)
    scope$upper<- ins
  }else if(is.list(scope)){
    if(length(scope)==1){
      scope0<- scope[[1]]
      scope$lower<- numeric(0)
      scope$upper<- scope0
    }else if(length(scope)==2){
      flg1<- any(is.na(match(scope[[1]],scope[[2]])),na.rm = FALSE)
      flg2<- any(is.na(match(scope[[2]],scope[[1]])),na.rm = FALSE)
      if(flg1&&flg2){
        stop("scope: wrong...")
      }else if(flg1){
        scope$lower<- scope[[2]]
        scope$upper<- scope[[1]]
      }else{
        scope$lower<- scope[[1]]
        scope$upper<- scope[[2]]
      }
    }else stop("scope: wrong...")
  }else{
    scope0<- scope
    scope<- list()
    scope$lower<- numeric(0)
    scope$upper<- scope0
  }

  if(ext) o<- mMove(obj,scope$upper,...) else o<- obj
  lik<- mlogLik(o)
  if(direction=="both"){
    yes<- TRUE
    if(missing(cv)){
      stop("\a\n  cv: should be positive but missing...\n\n")
    }
    while(yes){
      yes<- FALSE
      od<- mDrop1(o,scope$lower,ext=ext,...)
      lik1<- mlogLik(od$fit)
      if(od$drop & 2*(lik-lik1)<cv){
        o<- od$fit
        lik<- mlogLik(o)
        yes<- TRUE
      }else{
        oa<- mAdd1(o,scope$upper,ext=ext,...)
        lik1<- mlogLik(oa$fit)
        if(oa$add && 2*(lik1-lik)>cv){
          o<- oa$fit
          lik<- mlogLik(o)
          yes<- TRUE
        }else{
          od<- mDrop1(oa$fit,scope$lower,ext=ext,...)
          lik1<- mlogLik(od$fit)
          if(lik1>lik+1e-8){
             o<- od$fit
             lik<- mlogLik(o)
             yes<- TRUE
          }
        }
      }
    }
  }else if(direction=="backward"){
    if(missing(cv)){
      cv<- Inf
    }
    yes<- TRUE
    while(yes){
      yes<- FALSE
      od<- mDrop1(o,scope$lower,ext=ext,...)
      lik1<- mlogLik(od$fit)
      if(od$drop && 2*(lik-lik1)<cv){
        o<- od$fit
        lik<- mlogLik(o)
        yes<- TRUE
      }
    }
  }else if(direction=="forward"){
    if(missing(cv)){
      cv<- 0
    }
    yes<- TRUE
    while(yes){
      yes<- FALSE
      oa<- mAdd1(o,scope$upper,ext=ext,...)
      lik1<- mlogLik(oa$fit)
      if(oa$add && 2*(lik1-lik)>cv){
        o<- oa$fit
        lik<- mlogLik(o)
        yes<- TRUE
      }
    }
  }

  o
}

}
