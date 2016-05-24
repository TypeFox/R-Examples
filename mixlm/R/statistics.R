# $Id: statistics.R 38 2014-02-18 21:41:43Z khliland $

##################################
# Stepwise forward/backward
stepWise <- function(model, alpha.enter=0.15, alpha.remove=0.15, full=FALSE){
  # NB! The consistency of hiearchy is somewhat unclear in this method
  #     as add() and drop() work differently in this respect.
  cat("Stepwise regression (forward-backward), alpha-to-enter: ", alpha.enter, ", alpha-to-remove: ", alpha.remove, "\n\nFull model: ", sep="")
  full.formula <- as.formula(model$call$formula)        # Formula for full model
  print(full.formula)
  cat("\n")
  current.formula <- update(full.formula, ~ 1)          # Response against intercept formula
  current.model <- update(model, current.formula)       # ---------- || -----------  model
  n.effects <- length(attr(terms(model),"term.labels")) # Number of possible extra regressors
  n.obs <- length(model$residuals)
  S2 <- sum(model$residuals^2)/(n.obs-length(model$coefficients)) # Sum of squares for full model
  the.summary <- 'No iterations performed (re-run using extended output for details)'
  F <- 0
  i <- 1
  while(F<alpha.enter && n.effects>0){ # Add variables while extra regressors contribute
    added <- add1(current.model, full.formula, test="F") # All possible extended models
    F.added <- added[,6]            # p values of F statistic
    F <- min(F.added, na.rm = TRUE) # Minimum of p values
    n <- which(F.added == F)-1      # Number of best extra regressor
    if(full){ # Extended output
      cat("--= Step (forward)", i, "=--", "\n", sep=" ")
      print(added)
      cat("\n\n")}
    if(F<alpha.enter){
      current.formula <- update(current.formula, paste("~.+",rownames(added[n+1,,drop=FALSE]),sep="")) # Add extra regressor to formula
      current.model <- update(current.model, current.formula) # Update model with extra regressor
      n.effects <- n.effects - 1
      if(!full){
        p <- length(current.model$coefficients)
        if(i==1){ # Add current update to the final summary
          the.summary <- cbind(Step=1, InOut=1, added[n+1,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, added[n+1,5:6])
        } else{
          the.summary <- rbind(the.summary,cbind(Step=i, InOut=1, added[n+1,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, added[n+1,5:6]))}}
      i <- i+1
      
      
      # Remove extra regressors
      removed <- drop1(current.model, test="F")
      F.removed <- removed[,6]            # p values of F statistic
      F.r <- max(F.removed, na.rm = TRUE) # Maximum of p values
      n.r <- which(F.removed == F.r)-1    # Number of best removed regressor
      if(F.r>alpha.remove){
        current.formula <- update(current.formula, paste("~.-",rownames(removed[n.r+1,,drop=FALSE]),sep="")) # Add extra regressor to formula
        current.model <- update(current.model, current.formula) # Update model with extra regressor
        F <- 0
        if(!full){ # Add current update to the final summary
          p <- length(current.model$coefficients)
          if(i==1){
            the.summary <- cbind(Step=1, InOut=-1, removed[n.r+1,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, removed[n.r+1,5:6])}
          else{
            the.summary <- rbind(the.summary,cbind(Step=i, InOut=-1, removed[n.r+1,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, removed[n.r+1,5:6]))
          }
        } else {
          cat("--= Step (backward)", i, "=--", "\n", sep=" ")
          print(removed)
          cat("\n\n")
        }
        i <- i+1
      }
    }		
  }
  if(!full){
    if(i == 1){
      print(the.summary)}
    else {
      printCoefmat(the.summary,cs.ind=NULL)}}
  return(current.model)
}


##################################
# Stepwise backward/forward
stepWiseBack <- function(model, alpha.remove=0.15, alpha.enter=0.15, full=FALSE){
  cat("Stepwise regression (backward-forward), alpha-to-remove: ", alpha.remove, ", alpha-to-enter: ", alpha.enter, "\n\nFull model: ", sep="")
  full.formula <- as.formula(model$call$formula)        # Formula for full model
  print(full.formula)
  cat("\n")
  current.formula <- full.formula                       # Response against intercept formula
  current.model <- model						          # ---------- || -----------  model
  n.effects <- length(attr(terms(model),"term.labels")) # Number of possible regressors
  the.summary <- 'No iterations performed (re-run using extended output for details)'
  n.obs <- length(model$residuals)
  S2 <- sum(model$residuals^2)/(n.obs-length(model$coefficients)) # Sum of squares for full model
  F <- 1
  i <- 1
  while(F>alpha.remove && n.effects>0){ # Remove variables while some included regressor does not contribute
    removed <- drop1(current.model, test="F") # All possible reduced models
    F.removed <- removed[,6]          # p values of F statistic
    F <- max(F.removed, na.rm = TRUE) # Minimum of p values
    n <- which(F.removed == F)-1      # Number of best worst regressor
    if(full){ # Extended output
      cat("--= Step (backward)", i, "=--", "\n", sep=" ")
      print(removed)
      cat("\n\n")}
    if(F>alpha.remove){
      current.formula <- update(current.formula, paste("~.-",rownames(removed[n+1,,drop=FALSE]),sep="")) # Remove regressor from formula
      current.model <- update(current.model, current.formula) # Update model with extra regressor
      n.effects <- n.effects - 1
      if(!full){ # Add current update to the final summary
        p <- length(current.model$coefficients)
        if(i==1){
          the.summary <- cbind(Step=1, InOut=-1, removed[n+1,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, removed[n+1,5:6])}
        else{
          the.summary <- rbind(the.summary,cbind(Step=i, InOut=-1, removed[n+1,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, removed[n+1,5:6]))}}
      i <- i+1
      
      # Add extra regressors
      added <- add1(current.model, full.formula, test="F")
      F.added <- added[,6]              # p values of F statistic
      F.a <- min(F.added, na.rm = TRUE) # Maximum of p values
      n.a <- which(F.added == F.a)-1    # Number of best added regressor
      if(F.a<alpha.enter){
        current.formula <- update(current.formula, paste("~.+",rownames(added[n.a+1,,drop=FALSE]),sep="")) # Add extra regressor to formula
        current.model <- update(current.model, current.formula) # Update model with extra regressor
        F <- 0
        if(!full){ # Add current update to the final summary
          p <- length(current.model$coefficients)
          if(i==1){
            the.summary <- cbind(Step=1, InOut=1, added[n.a+1,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, added[n.a+1,5:6])}
          else{
            the.summary <- rbind(the.summary,cbind(Step=i, InOut=1, added[n.a+1,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, added[n.a+1,5:6]))}}
        else {
          cat("--= Step (forward)", i, "=--", "\n", sep=" ")
          print(added)
          cat("\n\n")}
        i <- i+1
      }
    }		
  }
  if(!full){
    if(i == 1){
      print(the.summary)}
    else {
      printCoefmat(the.summary,cs.ind=NULL)}}
  return(current.model)
}


##################################
## Best subsets
best.subsets <- function(model, nbest=5, nvmax, digits, force.in='NULL'){
  if(missing(nvmax)){
    nvmax <- length(attr(terms(model),"term.labels")) # Number of possible regressors
  }
  full.formula <- as.formula(model$call$formula)        # Formula for full model
  current.data <- model$call$data
  subsets <- eval(parse(text = paste("regsubsets(full.formula, data=", current.data, ", nbest=nbest, nvmax=nvmax, force.in=",force.in,")")))
  ss <- summary(subsets)
  if(missing(digits)){
    with(ss,print(data.frame(outmat,RSS=rss,R2=rsq,R2adj=adjr2,Cp=cp)))
  } else {
    with(ss,print(format(data.frame(outmat,RSS=rss,R2=rsq,R2adj=adjr2,Cp=cp), digits=digits)))
  }
}


##################################
# Forward addition
forward <- function(model, alpha=0.2, full=FALSE, force.in=NULL){
  cat("Forward selection, alpha-to-enter: ", alpha, "\n\nFull model: ", sep="")
  full.formula <- as.formula(model$call$formula)       # Formula for full model
  print(full.formula)
  cat("\n")
  if(is.null(force.in)){
    current.formula <- update(full.formula, ~ 1)         # Response against intercept formula
  } else {
    current.formula <- update(full.formula, paste('~',force.in, sep=""))}         # Response against intercept formula
  current.model <- update(model, current.formula)      # ---------- || -----------  model
  possible.effects <- attr(terms(model),"term.labels") # Vector of possible extra regressors
  the.summary <- 'No iterations performed (re-run using extended output for details)'
  n.obs <- length(model$residuals)
  S2 <- sum(model$residuals^2)/(n.obs-length(model$coefficients)) # Sum of squares for full model.
  F <- 0
  i <- 1
  the.summary <- NULL
  while(F<alpha && length(possible.effects)>0){ # Add variables while extra regressors contribute
    added <- add1(current.model, full.formula, test="F") # All possible extended models
    F.added <- added[,6]            # p values of F statistic
    F <- min(F.added, na.rm = TRUE) # Minimum of p values
    n <- which(F.added == F)      # Number of best extra regressor
    if(full){ # Extended output
      cat("--= Step", i, "=--", "\n", sep=" ")
      print(added)
      cat("\n\n")}
    if(F<alpha){
      current.formula <- update(current.formula, paste("~.+",attr(added,"row.names")[n],sep="")) # Add extra regressor to formula
      current.model <- update(current.model, current.formula) # Update model with extra regressor
      possible.effects <- possible.effects[-(n-1)]                # Remove extra regressor from possible regressors
      if(!full){ # Add current update to the final summary
        p <- length(current.model$coefficients)
        if(i==1){
          the.summary <- cbind(Step=1, added[n,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, added[n,5:6])}
        else{
          the.summary <- rbind(the.summary,cbind(Step=i, added[n,3:4], R2pred=R2pred(current.model), Cp=sum(current.model$residuals^2)/S2-n.obs+2*p, added[n,5:6]))}}
      i <- i+1
    }
  }
  if(!full){
    if(i == 1){
      print(the.summary)}
    else {
      printCoefmat(the.summary,cs.ind=NULL)}}
  return(current.model)
}


##################################
# Backward elimination
backward <- function(model, alpha=0.2, full=FALSE, hierarchy=TRUE, force.in=NULL){
  cat("Backward elimination, alpha-to-remove: ", alpha, "\n\nFull model: ", sep="")
  current.dropable <- current.formula <- as.formula(model$call$formula)    # Formula for full model
  print(current.formula)
  if(is.null(force.in)){
    possible.effects <- attr(terms(model),"term.labels") # Vector of possible removed regressors
  } else {
    current.dropable <- update(current.formula, paste("~.-", paste(force.in, collapse="-"), sep=""))
    possible.effects <- setdiff(attr(terms(model),"term.labels"), force.in)} # Vector of possible removed regressors
  cat("\n")
  n.obs <- length(model$residuals)
  S2 <- sum(model$residuals^2)/(n.obs-length(model$coefficients)) # Sum of squares for full model
  the.summary <- 'No iterations performed (re-run using extended output for details)'
  F <- Inf
  i <- 1
  while(F>alpha && length(possible.effects)>0){ # Remove variables until all regressor contributes
    if(hierarchy){
      removed <- drop1(model, setdiff(drop.scope(model),ifelse(is.null(force.in),"",strsplit(force.in,"\\+")[[1]])), test="F")} # All possible extended models
    else{
      removed <- drop1(model, current.dropable, test="F")} # All possible extended models
    F.removed <- removed[,6]          # p values of F statistic
    F <- max(F.removed, na.rm = TRUE) # Maximum of p values
    n <- which(F.removed == F)-1  # Number of best removed regressor
    if(full){ # Extended output
      cat("--= Step", i, "=--", "\n", sep=" ")
      print(removed)
      cat("\n\n")}
    if(F>alpha){
      current.formula  <- update(current.formula, paste("~.-",attr(removed,"row.names")[n+1],sep="")) # Remove regressor from formula
      current.dropable <- update(current.dropable, paste("~.-",attr(removed,"row.names")[n+1],sep="")) # Remove regressor from dropable
      model <- update(model, current.formula)  # Update model without removed regressor
      possible.effects <- possible.effects[-n] # Remove removed regressor from possible regressors
      if(!full){ # Add current update to the final summary
        p <- length(model$coefficients)
        if(i==1){
          the.summary <- cbind(Step=1, removed[n+1,3:4], R2pred=R2pred(model), Cp=sum(model$residuals^2)/S2-n.obs+2*p, removed[n+1,5:6])}
        else{
          the.summary <- rbind(the.summary,cbind(Step=i, removed[n+1,3:4], R2pred=R2pred(model), Cp=sum(model$residuals^2)/S2-n.obs+2*p, removed[n+1,5:6]))}}
      i <- i+1
    }
  }
  if(!full){
    if(i == 1){
      print(the.summary)}
    else {
      printCoefmat(the.summary,cs.ind=NULL)}}
  return(model)
}


#######################
## PRESS statistics
PRESS.res <- function(object=NULL, ncomp=NULL) {
  if(is.null(object) && "package:Rcmdr"%in%search()){
    try(eval(parse(text=paste("object <- ", activeModel(), sep=""))))
  }
  if(class(object)[1]=="lm" || class(object)[1]=="lmm" || class(objects)[1]=="glm"){
    return(residuals(object)/(1-lm.influence(object)$hat))
  }
  if(class(object)[1]=="mvr"){
    if(is.null(object$validation)){
      if(dim(object$model)[1] < 10){
        warning("Refitting with leave-one-out cross-validation as no cross-validation was done")
        object <- update(object,validation="LOO")
      } else {
        warning("Refitting with 10-fold cross-validation as no cross-validation was done")
        object <- update(object,validation="CV")
      }
    }
    if(is.null(ncomp)){
      ncomp <- object$validation$ncomp
    } else {
      ncomp <- min(ncomp, object$validation$ncomp)
    }
    return(model.response(model.frame(object))-object$validation$pred[,,ncomp])
  }
}
PRESS.pred <- function(object=NULL, ncomp=NULL) {
  if(is.null(object)){
    try(eval(parse(text=paste("object <- ", activeModel(), sep=""))))
  }
  if(class(object)[1]=="lm" || class(object)[1]=="lmm" || class(objects)[1]=="glm"){
    return(model.response(model.frame(object)) - residuals(object)/(1-lm.influence(object)$hat))
  }
  if(class(object)[1]=="mvr"){
    if(is.null(object$validation)){
      if(dim(object$model)[1] < 10){
        warning("Refitting with leave-one-out cross-validation as no cross-validation was done")
        object <- update(object,validation="LOO")
      } else {
        warning("Refitting with 10-fold cross-validation as no cross-validation was done")
        object <- update(object,validation="CV")
      }
    }
    if(is.null(ncomp)){
      ncomp <- object$validation$ncomp
    } else {
      ncomp <- min(ncomp, object$validation$ncomp)
    }
    return(object$validation$pred[,,ncomp])
  }
}
PRESS <- function(object=NULL) {
  if(is.null(object) && "package:Rcmdr"%in%search()){
    try(eval(parse(text=paste("object <- ", activeModel(), sep=""))))
  }
  the.PRESS <- NULL
  if(class(object)[1]=="mvr"){ # PCR/PLSR
    temp.model <- object
    hasVal <- !is.null(temp.model$validation)
    hasLOO <- logical(0)
    if(hasVal){ # Check if validated
      hasLOO <- attr(temp.model$validation$segments, 'type')=='leave-one-out'
    }
    if(length(hasLOO)!=1 || hasLOO==FALSE){ # Wrong/no cross-validation
      temp.model <- update(temp.model,validation="LOO")
    }
    n <- dim(temp.model$scores)[1]
    comp0 <- temp.model$validation$PRESS0; names(comp0) <- "(intercept)"
    the.PRESS <- c(comp0,temp.model$validation$PRESS[1,])
  }
  if(class(object)[1]=="lm" || class(object)[1]=="lmm"){ # Linear regression
    the.PRESS <- sum(residuals(object)^2/(1-lm.influence(object)$hat)^2)
  }
  the.PRESS
}
R2pred <- function(object=NULL) {
  if(is.null(object) && "package:Rcmdr"%in%search()){
    try(eval(parse(text=paste("object <- ", activeModel(), sep=""))))
  }
  R2pred <- NULL
  if(class(object)[1]=="mvr"){ # PCR/PLSR
    temp.model <- object
    hasVal <- !is.null(temp.model$validation)
    hasLOO <- logical(0)
    if(hasVal){ # Check if validated
      hasLOO <- attr(temp.model$validation$segments, 'type')=='leave-one-out'
    }
    if(length(hasLOO)!=1 || hasLOO==FALSE){ # Wrong/no cross-validation
      temp.model <- update(temp.model,validation="LOO")
    }
    n <- dim(temp.model$scores)[1]
    comp0 <- temp.model$validation$PRESS0; names(comp0) <- "(intercept)"
    R2pred <- 1-c(comp0,temp.model$validation$PRESS[1,])/((n-1)*var(temp.model$model[,1]))
  }
  if(class(object)[1]=="lm" || class(object)[1]=="lmm"){ # Linear regression
    R2pred <- 1 - sum(residuals(object)^2/(1-lm.influence(object)$hat)^2) /
      (var(object$model[,1])*(length(object$model[,1])-1))
  }
  R2pred
}


####################
## Confusion matrix
confusion <- function(true, predicted){
  n.lev <- length(levels(predicted))
  a <- table(Predicted=predicted,True=true)
  b <- rbind(a,Total=apply(a,2,sum),Correct=diag(a),Proportion=diag(a)/apply(a,2,sum))
  aa <- dimnames(a)
  aa$Predicted <- c(aa$Predicted,"Total","Correct","Proportion")
  dimnames(b) <- aa 
  print(b[-(n.lev+3),])
  cat("\nProportions correct\n")
  print(b[n.lev+3,])
  cat(paste('\nN correct/N total = ', sum(b[n.lev+2,]), '/', sum(b[n.lev+1,]), ' = ', format(sum(b[n.lev+2,])/sum(b[n.lev+1,])), '\n', sep=''))
}


#################################
# ANOVA for regression
anova_reg <- function(lm.object){
  anova.result <- anova(lm.object)
  p <- dim(anova.result)[1]
  new.anova <- anova.result[c(1,p),]
  new.anova[1,1] <- sum(anova.result[1:(p-1),1])
  new.anova[1,2] <- sum(anova.result[1:(p-1),2])
  new.anova[1,3] <- new.anova[1,2]/new.anova[1,1]
  new.anova[1,4] <- new.anova[1,3]/new.anova[2,3]
  new.anova[1,5] <- pf(new.anova[1,4], new.anova[1,1], new.anova[2,1], lower.tail=FALSE)
  if(p>2)
    rownames(new.anova)[1] <- "Regression"
  attributes(new.anova)$heading <- "Analysis of Variance Table"
  new.anova
}


################################
# RMSEP for lm and mvr
RMSEP <- function(object){
  if(is.null(object) && "package:Rcmdr"%in%search()){
    try(eval(parse(text=paste("object <- ", activeModel(), sep=""))))
  }
  the.RMSEP <- NULL
  if(class(object)[1]=="mvr"){ # PCR/PLSR
    the.RMSEP <- pls::RMSEP(object, estimate="all")
  }
  if(class(object)[1]=="lmm" || class(object)[1]=="lm"){ # Linear regression
    the.RMSEP <- sqrt(mean(residuals(object)^2/(1-lm.influence(object)$hat)^2))
  }
  the.RMSEP
}
rmsep <- function(object){
  mixlm::RMSEP(object)
}

################################
# One sample proportion test
prop.test.ordinary <- function (x, n, p = NULL, alternative = c("two.sided", "less", 
                                                                "greater"), conf.level = 0.95, correct = TRUE, pooled=TRUE)
{
  DNAME <- deparse(substitute(x))
  if (is.table(x) && length(dim(x)) == 1L) {
    if (dim(x) != 2L) 
      stop("table 'x' should have 2 entries")
    l <- 1
    n <- sum(x)
    x <- x[1L]
  }
  else if (is.matrix(x)) {
    if (ncol(x) != 2L) 
      stop("'x' must have 2 columns")
    l <- nrow(x)
    n <- rowSums(x)
    x <- x[, 1L]
  }
  else {
    DNAME <- paste(DNAME, "out of", deparse(substitute(n)))
    if ((l <- length(x)) != length(n)) 
      stop("'x' and 'n' must have the same length")
  }
  OK <- complete.cases(x, n)
  x <- x[OK]
  n <- n[OK]
  if ((k <- length(x)) < 1L) 
    stop("not enough data")
  if (any(n <= 0)) 
    stop("elements of 'n' must be positive")
  if (any(x < 0)) 
    stop("elements of 'x' must be nonnegative")
  if (any(x > n)) 
    stop("elements of 'x' must not be greater than those of 'n'")
  if (is.null(p) && (k == 1)) 
    p <- 0.5
  if (!is.null(p)) {
    DNAME <- paste(DNAME, ", null ", ifelse(k == 1, "probability ", 
                                            "probabilities "), deparse(substitute(p)), sep = "")
    if (length(p) != l) 
      stop("'p' must have the same length as 'x' and 'n'")
    p <- p[OK]
    if (any((p <= 0) | (p >= 1))) 
      stop("elements of 'p' must be in (0,1)")
  }
  alternative <- match.arg(alternative)
  if (k > 2 || (k == 2) && !is.null(p)) 
    alternative <- "two.sided"
  if ((length(conf.level) != 1L) || is.na(conf.level) || (conf.level <= 
                                                          0) || (conf.level >= 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  correct <- as.logical(correct)
  ESTIMATE <- x/n
  names(ESTIMATE) <- if (k == 1) 
    "p"
  else paste("prop", 1L:l)[OK]
  NVAL <- p
  CINT <- NULL
  YATES <- ifelse(correct && (k <= 2), 0.5, 0)
  if (k == 1) {
    z <- ifelse(alternative == "two.sided", qnorm((1 + conf.level)/2), 
                qnorm(conf.level))
    YATES <- min(YATES, abs(x - n * p))
    p.c <- ESTIMATE + YATES/n
    p.u <- if (p.c >= 1) 
      1
    else (p.c + z * sqrt(p.c * (1 - p.c)/n))
    p.c <- ESTIMATE - YATES/n
    p.l <- if (p.c <= 0) 
      0
    else (p.c - z * sqrt(p.c * (1 - p.c)/n))
    CINT <- switch(alternative, two.sided = c(max(p.l, 0), 
                                              min(p.u, 1)), greater = c(max(p.l, 0), 1), less = c(0, 
                                                                                                  min(p.u, 1)))
  }
  else if ((k == 2) & is.null(p)) {
    DELTA <- ESTIMATE[1L] - ESTIMATE[2L]
    YATES <- min(YATES, abs(DELTA)/sum(1/n))
    WIDTH <- (switch(alternative, two.sided = qnorm((1 + 
                                                       conf.level)/2), qnorm(conf.level)) * sqrt(sum(ESTIMATE * 
                                                                                                       (1 - ESTIMATE)/n)) + YATES * sum(1/n))
    CINT <- switch(alternative, two.sided = c(max(DELTA - 
                                                    WIDTH, -1), min(DELTA + WIDTH, 1)), greater = c(max(DELTA - 
                                                                                                          WIDTH, -1), 1), less = c(-1, min(DELTA + WIDTH, 1)))
  }
  if (!is.null(CINT)) 
    attr(CINT, "conf.level") <- conf.level
  METHOD <- paste(ifelse(k == 1, "1-sample proportions test", 
                         paste(k, "-sample test for ", ifelse(is.null(p), "equality of", 
                                                              "given"), " proportions", sep = "")), ifelse(YATES, 
                                                                                                           "with", "without"), "continuity correction")
  if (is.null(p)) {
    p <- sum(x)/sum(n)
    PARAMETER <- k - 1
  }
  else {
    PARAMETER <- k
    names(NVAL) <- names(ESTIMATE)
  }
  names(PARAMETER) <- "df"
  if(k == 1 || pooled == TRUE)
    x <- cbind(x, n - x)
  E <- cbind(n * p, n * (1 - p))
  if (any(E < 5)) 
    warning("Approximation may be incorrect")
  STATISTIC <- sum((abs(x - E) - YATES)^2/E)
  if (k == 1) 
    z <- sign(ESTIMATE - p) * sqrt(STATISTIC)
  else { 
    if(!exists("DELTA"))
      DELTA <- 1
    if(pooled)
      z <- sign(DELTA) * sqrt(STATISTIC)
    else {
      p <- x/n
      z <- sign(DELTA)*(abs(DELTA) + YATES * sum(1/n))/sqrt(sum(p*(1-p)/n))
      STATISTIC <- z^2
    }
  }
  if (alternative == "two.sided") {
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  } else {
    PVAL <- pnorm(z, lower.tail = (alternative == "less"))
  }
  STATISTIC <- z
  names(STATISTIC) <- "z"
  PARAMETER <- Inf
  names(PARAMETER) <- "df"
  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
               p.value = as.numeric(PVAL), estimate = ESTIMATE, null.value = NVAL, 
               conf.int = CINT, alternative = alternative, method = METHOD, 
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}



################################
# Standardized Pearson residuals
spearson <- function(object){
  residuals(object, type="pearson")/sqrt(1-hatvalues(object))
}


#####################################
# z and t tests for summarized data
z_test_sum <- function(means, sds, ns, alternative = c("two.sided", "less", "greater"),
                       mu = 0, var.equal = FALSE, conf.level = 0.95, z.test=TRUE, ...){
  t_test_sum(means, sds, ns, alternative,
             mu, var.equal, conf.level, z.test=TRUE, ...)
}
t_test_sum <- function(means, sds, ns, alternative = c("two.sided", "less", "greater"),
                       mu = 0, var.equal = FALSE, conf.level = 0.95, z.test=FALSE, ...)
{
  alternative <- match.arg(alternative)
  
  if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if(!missing(conf.level) &&
     (length(conf.level) != 1 || !is.finite(conf.level) ||
      conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if( length(means)==2 ) {
    dname <- "two samples"#paste(deparse(substitute(x)),"and",
    #  deparse(substitute(y)))
  }
  else {
    dname <- "one sample" #deparse(substitute(x))
  }
  nx <- ns[1]
  mx <- means[1]
  vx <- sds[1]^2
  estimate <- mx
  if(length(means)==1) {
    if(nx < 2) stop("not enough 'x' observations")
    df <- ifelse(z.test,Inf,nx-1)
    stderr <- sqrt(vx/nx)
    if(stderr < 10 *.Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    tstat <- (mx-mu)/stderr
    method <- ifelse(z.test,"One Sample z-test","One Sample t-test")
    names(estimate) <- "mean of x"
  } else {
    ny <- ns[2]
    if(nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if(var.equal && nx+ny < 3) stop("not enough observations")
    my <- means[2]
    vy <- sds[2]^2
    method <- paste(if(!var.equal)"Welch", ifelse(z.test,"Two Sample z-test","Two Sample t-test"))
    estimate <- c(mx,my)
    names(estimate) <- c("mean of x","mean of y")
    if(var.equal) {
      df <- nx+ny-2
      v <- 0
      if(nx > 1) v <- v + (nx-1)*vx
      if(ny > 1) v <- v + (ny-1)*vy
      v <- v/df
      stderr <- sqrt(v*(1/nx+1/ny))
      df <- ifelse(z.test,Inf,nx+ny-2)
    } else {
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- ifelse(z.test,Inf,stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1)))
    }
    if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df) )
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- ifelse(z.test,"z","t")
  names(df) <- "df"
  names(mu) <- if(length(means)==2) "difference in means" else "mean"
  attr(cint,"conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               conf.int=cint, estimate=estimate, null.value = mu,
               alternative=alternative,
               method=method, data.name=dname)
  class(rval) <- "htest"
  return(rval)
}

#####################################
# z tests for ordinary data
z_test <- function(x, ...) UseMethod("z_test")
z_test.formula <-
  function(formula, data, subset, na.action, ...)
  {
    if(missing(formula)
       || (length(formula) != 3L)
       || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
      stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    m[[1L]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
      stop("grouping factor must have exactly 2 levels")
    DATA <- split(mf[[response]], g)
    names(DATA) <- c("x", "y")
    y <- do.call("z_test", c(DATA, list(...)))
    y$data.name <- DNAME
    if(length(y$estimate) == 2L)
      names(y$estimate) <- paste("mean in group", levels(g))
    y
  }

z_test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                           mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95, sds=NULL,
                           ...)
{
  z.test <- TRUE
  alternative <- match.arg(alternative)
  
  if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if(!missing(conf.level) &&
     (length(conf.level) != 1 || !is.finite(conf.level) ||
      conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if( !is.null(y) ) {
    dname <- paste(deparse(substitute(x)),"and",
                   deparse(substitute(y)))
    if(paired)
      xok <- yok <- complete.cases(x,y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    dname <- deparse(substitute(x))
    if( paired ) stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if( paired ) {
    x <- x-y
    y <- NULL
  }
  nx <- length(x)
  mx <- mean(x)
  vx <- ifelse(z.test,sds[1]^2,var(x))
  estimate <- c(mx,sqrt(vx))
  if(is.null(y)) {
    if(nx < 2) stop("not enough 'x' observations")
    df <- ifelse(z.test,Inf,nx-1)
    stderr <- sqrt(vx/nx)
    if(stderr < 10 *.Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    tstat <- (mx-mu)/stderr
    method <- ifelse(paired,ifelse(z.test,"Paired z-test","Paired t-test"),ifelse(z.test,"One Sample z-test","One Sample t-test"))
    names(estimate) <- c(ifelse(paired,"mean of the differences","mean of x"),ifelse(paired," std.dev. of the differences"," std.dev. of x"))
  } else {
    ny <- length(y)
    if(nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if(var.equal && nx+ny < 3) stop("not enough observations")
    my <- mean(y)
    vy <- ifelse(z.test,sds[2]^2,var(y))
    method <- paste(if(!var.equal)"Welch", ifelse(z.test,"Two Sample z-test","Two Sample t-test"))
    if(var.equal) {
      df <- nx+ny-2
      v <- 0
      if(nx > 1) v <- v + (nx-1)*vx
      if(ny > 1) v <- v + (ny-1)*vy
      v <- v/df
      stderr <- sqrt(v*(1/nx+1/ny))
      df <- ifelse(z.test,Inf,nx+ny-2)
      estimate <- c(mx,my,sqrt(v))
      names(estimate) <- c("mean of x"," mean of y"," pooled std.dev.")
    } else {
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- ifelse(z.test,Inf,stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1)))
      estimate <- c(mx,my,sqrt(vx),sqrt(vy))
      names(estimate) <- c("mean of x"," mean of y", " std.dev. of x", " std.dev. of y")
    }
    if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df) )
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- ifelse(z.test,"z","t")
  names(df) <- "df"
  names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
  attr(cint,"conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               conf.int=cint, estimate=estimate, null.value = mu,
               alternative=alternative,
               method=method, data.name=dname)
  class(rval) <- "htest"
  return(rval)
}

#####################################
# t tests for ordinary data (extended output)
t_test <- function(x, ...) UseMethod("t_test")
t_test.default <-
  function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
           mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
           ...)
  {
    alternative <- match.arg(alternative)
    
    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
      stop("'mu' must be a single number")
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
      stop("'conf.level' must be a single number between 0 and 1")
    if( !is.null(y) ) {
      dname <- paste(deparse(substitute(x)),"and",
                     deparse(substitute(y)))
      if(paired)
        xok <- yok <- complete.cases(x,y)
      else {
        yok <- !is.na(y)
        xok <- !is.na(x)
      }
      y <- y[yok]
    }
    else {
      dname <- deparse(substitute(x))
      if( paired ) stop("'y' is missing for paired test")
      xok <- !is.na(x)
      yok <- NULL
    }
    x <- x[xok]
    if( paired ) {
      x <- x-y
      y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    estimate <- c(mx,sqrt(vx))
    if(is.null(y)) {
      if(nx < 2) stop("not enough 'x' observations")
      df <- nx-1
      stderr <- sqrt(vx/nx)
      if(stderr < 10 *.Machine$double.eps * abs(mx))
        stop("data are essentially constant")
      tstat <- (mx-mu)/stderr
      method <- ifelse(paired,"Paired t-test","One Sample t-test")
      names(estimate) <- c(ifelse(paired,"mean of the differences","mean of x"),ifelse(paired," std.dev. of the differences", " std.dev. of x"))
    } else {
      ny <- length(y)
      if(nx < 1 || (!var.equal && nx < 2))
        stop("not enough 'x' observations")
      if(ny < 1 || (!var.equal && ny < 2))
        stop("not enough 'y' observations")
      if(var.equal && nx+ny < 3) stop("not enough observations")
      my <- mean(y)
      vy <- var(y)
      method <- paste(if(!var.equal)"Welch", "Two Sample t-test")
      estimate <- c(mx,my,sqrt(vx),sqrt(vy))
      names(estimate) <- c("mean of x"," mean of y"," std.dev. of x"," std.dev. of y")
      if(var.equal) {
        df <- nx+ny-2
        v <- 0
        if(nx > 1) v <- v + (nx-1)*vx
        if(ny > 1) v <- v + (ny-1)*vy
        v <- v/df
        stderr <- sqrt(v*(1/nx+1/ny))
        estimate <- c(mx,my,sqrt(v))
        names(estimate) <- c("mean of x"," mean of y"," pooled std.dev.")
      } else {
        stderrx <- sqrt(vx/nx)
        stderry <- sqrt(vy/ny)
        stderr <- sqrt(stderrx^2 + stderry^2)
        df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))
      }
      if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
        stop("data are essentially constant")
      tstat <- (mx - my - mu)/stderr
    }
    if (alternative == "less") {
      pval <- pt(tstat, df)
      cint <- c(-Inf, tstat + qt(conf.level, df) )
    }
    else if (alternative == "greater") {
      pval <- pt(tstat, df, lower.tail = FALSE)
      cint <- c(tstat - qt(conf.level, df), Inf)
    }
    else {
      pval <- 2 * pt(-abs(tstat), df)
      alpha <- 1 - conf.level
      cint <- qt(1 - alpha/2, df)
      cint <- tstat + c(-cint, cint)
    }
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
    attr(cint,"conf.level") <- conf.level
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
                 conf.int=cint, estimate=estimate, null.value = mu,
                 alternative=alternative,
                 method=method, data.name=dname)
    class(rval) <- "htest"
    return(rval)
  }

t_test.formula <-
  function(formula, data, subset, na.action, ...)
  {
    if(missing(formula)
       || (length(formula) != 3L)
       || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
      stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    m[[1L]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2L)
      stop("grouping factor must have exactly 2 levels")
    DATA <- split(mf[[response]], g)
    names(DATA) <- c("x", "y")
    y <- do.call("t_test", c(DATA, list(...)))
    y$data.name <- DNAME
    if(length(y$estimate) == 2L)
      names(y$estimate) <- paste("mean in group", levels(g))
    y
  }



#################################################
# Changed summary.lm print-out. Replaces "Residual standard error" by "s". Also extra line shift.
print.summary.lmm <- function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
                               signif.stars = getOption("show.signif.stars"), ...) 
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  cat(if (!is.null(x$w) && diff(range(x$w))) 
    "Weighted ", "Residuals:\n", sep = "")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L) 
      structure(apply(t(resid), 1L, quantile), dimnames = list(nam, 
                                                               dimnames(resid)[[2L]]))
    else {
      zz <- zapsmall(quantile(resid), digits + 1)
      structure(zz, names = nam)
    }
    print(rq, digits = digits, ...)
  }
  else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  }
  else {
    cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!\n")
  }
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    if (nsingular <- df[3L] - df[1L]) 
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                                                              colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
  }
  cat("\ns:", format(signif(x$sigma, 
                            digits)), "on", rdf, "degrees of freedom\n")
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared:", formatC(x$r.squared, digits = digits))
    cat(",\nAdjusted R-squared:", formatC(x$adj.r.squared, 
                                          digits = digits), "\nF-statistic:", formatC(x$fstatistic[1L], 
                                                                                      digits = digits), "on", x$fstatistic[2L], "and", 
        x$fstatistic[3L], "DF,  p-value:", format.pval(pf(x$fstatistic[1L], 
                                                          x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE), 
                                                       digits = digits), "\n")
  }
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1L) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2, 
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}


#################################################
## Property plots for relevant component analysis
plotprops <- function(Y,X, doscaleX=FALSE, docenterX=TRUE, ncomp, subset){
  n <- dim(X)[1]
  p <- dim(X)[2]
  ncomp <- min(ncomp,min(n,p))
  if(docenterX) ncomp <- ncomp-1
  if(ncomp<1)stop("Centering requires at least 2 components")
  if(missing(subset)) subset <- 1:n
  X <- scale(X[subset,], center=docenterX, scale=doscaleX)
  Y <- matrix(Y, ncol=1)[subset,,drop=F]
  svdres <- svd(X)
  eigval <- (svdres$d^2)/(svdres$d^2)[1]
  Z <- X%*%svdres$v
  covs <- cov(Y, Z)
  covs <- abs(covs)/max(abs(covs))
  par(mar=c(5.1, 4.1, 4.1, 4.1))
  plot(1:ncomp, eigval[1:ncomp], type="h", lwd=2, xlab="Component", ylab="Scaled eigenvalue", axes=FALSE, main="Property plot")
  points(1:ncomp, covs[1:ncomp], type="p", pch=20, cex=2, col=2)
  axis(1)
  axis(2,at=seq(0,1,0.1), labels=as.character(seq(0,1,0.1)))
  axis(4,at=seq(0,1,0.1), labels=as.character(seq(0,1,0.1)))
  mtext("Scaled covariance",side=4, line=3)
  box()
}


#################################################
## Tally of discrete variable
tally <- function(x){
  out  <- table(factor(x))
  perc <- 100*out/sum(out)
  cbind(Count=out,CumCount=cumsum(out),Percent=round(perc,2),CumPercent=round(cumsum(perc),2))
}


# Kommenter, flette inn generalTukey med effect
simple.glht <- function(mod, effect, corr = c("Tukey","Bonferroni","Fisher"), level = 0.95, ...) {
  if(missing(corr)){
    corr <- "Tukey"
  }
  random <- ifelse(is.null(mod$random),FALSE,TRUE)
  # mod <- aov(mod)
  if(random && corr!="Tukey")
    stop("Only Tukey correction supported for mixed models.")
  if(corr=="Tukey"){
    if(grepl(":", effect))
      stop("Only main effects supported for Tukey.")
    if(random && effect%in%mod$random$random)
      stop("Only fixed effects supported for mixed models.")
  }
  chkdots <- function(...) {
    
    lst <- list(...)
    if (length(lst) > 0) {
      warning("Argument(s) ", sQuote(names(lst)), " passed to ", sQuote("..."), 
              " are ignored", call. = TRUE)
    }
  }
  generalTukey <- function(mod,effect,random, ...){
    warn <- options("warn")
    options(warn=-1)
    if(grepl(":", effect)){
      pro  <- sub(":", "*", effect)
      prox <- sub(":", "_", effect)
      cola <- sub(":", "",  effect)
      spli <- unlist(strsplit(effect,":"))
      data <- model.frame(mod)
      eval(parse(text=paste("data$",prox," <- with(data, interaction(", paste(spli,collapse=",",sep=""),", sep=':'))")))
      mod <- update(mod, formula(paste(".~-", pro, "+",prox)), data=data)
      ret <- list()
      if(random){
        ret$res <- TukeyMix(mod,effect,level)
      } else {
        ret$res <- TukeyFix(mod,effect,level)
      }
      ret$model <- mod
    } else {
      ret <- list()
      if(random){
        ret$res <- TukeyMix(mod,effect,level)
        ret$model <- mod
      } else {
        if(corr == "Tukey"){
          # ret$res <- TukeyFix(mod,effect,level)
          ret <- eval(parse(text=paste("glht(mod, linfct=mcp(",effect,"='Tukey'),...)")))
          ret$model <- mod
        } else { # This will be overwritten.
          ret <- eval(parse(text=paste("glht(mod, linfct=mcp(",effect,"='Tukey'),...)")))
        }
      }
    }
    options(warn=warn$warn)
    ret
  }
  object <- generalTukey(mod,effect,random,...)
  
  chkdots(...)
  
  test <- switch(corr,
                 Tukey      = adjusted(),
                 Bonferroni = adjusted("bonferroni"),
                 Fisher     = univariate())
  calpha <- switch(corr,
                   Tukey      = adjusted_calpha(),
                   Bonferroni = adjusted_calpha("bonferroni"),
                   Fisher     = univariate_calpha())
  
  if(corr == "Tukey" && random){
    object$test <- 0
  } else {
    ts <- test(object)
    object$test <- ts
  }
  
  if(corr == "Tukey" && random){
    type <- "Tukey"
    object$type <- type
  } else {
    type <- attr(calpha, "type")
  }
  if (is.function(calpha))
    if(corr == "Tukey"){
      calpha <- 0
    } else {
      calpha <- calpha(object, level)
    }
  if (!is.numeric(calpha) || length(calpha) != 1)
    stop(sQuote("calpha"), " is not a scalar")
  if(corr == "Tukey" && random){
    error <- attr(object,"error")
  } else{
    error <- attr(calpha, "error")
  }
  attributes(calpha) <- NULL
  
  if(corr != "Tukey" || !random){
    betahat <- coef(object)
    ses <- sqrt(diag(vcov(object)))
    switch(object$alternative, "two.sided" = {
      LowerCL <- betahat - calpha * ses
      UpperCL <- betahat + calpha * ses
    }, "less" = {
      LowerCL <- rep(-Inf, length(ses))
      UpperCL <- betahat + calpha * ses
    }, "greater" = {
      LowerCL <- betahat + calpha * ses
      UpperCL <- rep( Inf, length(ses))
    })
    
    ci <- cbind(LowerCL, UpperCL)
    colnames(ci) <- c("lower", "upper")
    object$confint <- cbind(betahat, ci)
    colnames(object$confint) <- c("Estimate", "lwr", "upr")
    attr(object$confint, "conf.level") <- level
    attr(object$confint, "calpha") <- calpha
    attr(object$confint, "error") <- error
  } else {
    object$confint <- object$res[,c(2,1,3)]
    colnames(object$confint) <- c("Estimate","lwr","upr")
    attr(object$confint, "conf.level") <- level
  }
  
  if (is.null(type)) type <- "univariate"
  attr(object, "type") <- type
  attr(object, "corr") <- corr
  if(random){
    attr(object, "random") <- TRUE
  }
  if(corr == "Tukey"){
    object$focus  <- effect
    class(object) <- "summary.glht"
  } else {
    class(object) <- switch(class(ts), "mtest" = "summary.glht",
                            "gtest" = "summary.gtest")
  }
  class(object) <- c("simple.glht", "confint.glht", class(object), "glht")
  return(object)
}

# Print method
print.simple.glht <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\n\t", "Simultaneous Confidence Intervals and Tests for General Linear Hypotheses\n\n")
  if (!is.null(x$type)){
    cat("Multiple Comparisons of Means:", attr(x, "corr"), "Contrasts\n\n\n")
    level <- attr(x$confint, "conf.level")
  }
  attr(x$confint, "conf.level") <- NULL
  cat("Fit: ")
  if (isS4(x$model)) {
    print(x$model@call)
  } else {
    print(x$model$call)
  }
  cat("\n")
  error <- attr(x$confint, "error")
  if (!is.null(error) && error > .Machine$double.eps)
    digits <- min(digits, which.min(abs(1 / error - (10^(1:10)))))
  if(attr(x,'type') != "Tukey"){
    cat("Quantile =", round(attr(x$confint, "calpha"), digits))    
  } else {
    cat("Quantile =", round(attr(x$res, "quant"), digits), "\n")
    if(length(attr(x$res,"minSignDiff"))>0)
      cat("Minimum significant difference =", round(attr(x$res,"minSignDiff"),digits))
  }
  cat("\n")
  if (attr(x, "type") == "adjusted") {
    cat(paste(level * 100, 
              "% family-wise confidence level\n", sep = ""), "\n")
  } else {
    cat(paste(level * 100, 
              "% confidence level\n", sep = ""), "\n")
  }
  
  ### <FIXME>: compute coefmat in summary.glht for easier access???
  if(attr(x,'type') != "Tukey"){
    pq <- x$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    error <- attr(pq$pvalues, "error")
    pname <- switch(x$alternativ,
                    "less" = paste("Pr(<", ifelse(x$df == 0, "z", "t"), ")", sep = ""),
                    "greater" = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""),
                    "two.sided" = paste("Pr(>|", ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
    colnames(mtests) <- c("Estimate", "Std. Error",
                          ifelse(x$df == 0, "z value", "t value"), pname)
    type <- pq$type
    alt <- switch(x$alternative,
                  "two.sided" = "==", "less" = ">=", "greater" = "<=")
  } else {
    mtests <- x$res
    error <- 0
    alt <- "=="
    x$rhs <- 0
    type <- "single-step"
  }
  ### print p values according to simulation precision
  if (!is.null(error) && error > .Machine$double.eps) {
    sig <- which.min(abs(1 / error - (10^(1:10))))
    sig <- 1 / (10^sig)
  } else {
    sig <- .Machine$double.eps
  }
  cat("Linear Hypotheses:\n")
  #   alt <- switch(x$alternative,
  #                 "two.sided" = "==", "less" = ">=", "greater" = "<=")
  ### </FIXME>
  
  if(attr(x,'type') != "Tukey"){
    rownames(mtests) <- paste(rownames(mtests), alt, x$rhs)
    rownames(x$confint) <- paste(rownames(x$confint), alt, x$rhs)
    mtests <- cbind(x$confint,mtests[,2:4,drop=FALSE])
  }
  printCoefmat(mtests, digits = digits, 
               has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)
  
  switch(type, 
         "univariate" = cat("(Univariate p values reported)"),
         "single-step" = cat("(Adjusted p values reported -- single-step method)"),
         "Shaffer" = cat("(Adjusted p values reported -- Shaffer method)"),
         "Westfall" = cat("(Adjusted p values reported -- Westfall method)"),
         cat("(Adjusted p values reported --", type, "method)")
  )
  cat("\n\n")
  if (!is.balanced(x$model)) {
    cat("\nWARNING: Unbalanced data may lead to poor estimates\n")
  }
  invisible(x)                    
}

# Internal Tukey calculations for mixed models
TukeyMix <- function(mod, eff, level=0.95){
  object <- Anova(mod,type=3)
  if(!(eff%in%rownames(object$anova)))
    stop(paste(eff, ' not among model effects', sep=""))
  if(eff%in%object$random.effects)
    stop(paste(eff, ' is random', sep=""))
  data <- model.frame(mod)
  effVals <- mod$model[[eff]]
  effLevs <- sort(levels(effVals))
  weight  <- 2*eval(parse(text=paste("mean(1/xtabs(~",eff,",data=data))",sep="")))
  
  resp    <- model.response(data)
  means   <- tapply(resp,effVals,mean)
  mname   <- names(means)
  df      <- object$denom.df[[eff]]
  error   <- object$errors[match(eff,names(object$denom.df))]
  SE      <- sqrt(error*weight)
  if(df > 1){
    quant   <- qtukey(level,length(means),df)/sqrt(2)
  } else {
    warning('P-values for Tukey\'s HSD is not available for error df=1.')
    if(length(means) >= 2 & length(means) <= 20 & level %in% c(0.9,0.95,0.975,0.99)){
      quant <- qtukey1df[paste(length(means)),paste(level)]/sqrt(2)
    } else {
      quant <- NaN
    }
  }
  width   <- quant*SE
  
  n   <- choose(length(means),2)
  out <- data.frame(Lower=numeric(n), Diff=numeric(n), Upper=numeric(n), SE=numeric(n), T=numeric(n), 'P(>t)'=numeric(n))
  rownames(out) <- paste(1:n)
  k <- 1
  for(i in 1:(length(means)-1)){
    for(j in (i+1):length(means)){
      mij <- means[i]-means[j]
      p <- ptukey(abs(mij)/(SE/sqrt(2)),length(means),df,lower.tail=FALSE)
      lower <- mij-width
      upper <- mij+width
      out[k,] <- c(lower, mij, upper, SE, mij/SE, p)
      rownames(out)[k] <- paste(mname[i],'-',mname[j],sep="")
      k <- k+1
    }
  }
  colnames(out) <- c("Lower", "Center", "Upper", "Std.Err", "t value","P(>t)")
  attr(out,"minSignDiff") <- width
  attr(out,"means") <- means
  attr(out,"level") <- level
  attr(out,"effVals") <- effVals
  attr(out,"resp") <- resp
  attr(out,"quant") <- quant
  attr(out,"error") <- SE
  attr(out,"df") <- df
  class(out) <- c("TukeyMix","data.frame")
  out
}

# Internal Tukey calculations for fixed models
TukeyFix <- function(mod, eff, level=0.95){
  object <- Anova(mod,type=3)
  if(!(eff%in%rownames(object)))
    stop(paste(eff, ' not among model effects', sep=""))
  data <- model.frame(mod)
  effVals <- mod$model[[eff]]
  effLevs <- sort(levels(effVals))
  weight  <- 2*eval(parse(text=paste("mean(1/xtabs(~",eff,",data=data))",sep="")))
  
  resp    <- model.response(data)
  means   <- tapply(resp,effVals,mean)
  mname   <- names(means)
  n       <- dim(object)[1]
  df      <- object[n,'Df']
  error   <- object[n,'Sum Sq']/df
  SE      <- sqrt(error*weight)
  if(df > 1){
    quant   <- qtukey(level,length(means),df)/sqrt(2)
  } else {
    warning('P-values for Tukey\'s HSD is not available for error df=1.')
    if(length(means) >= 2 & length(means) <= 20 & level %in% c(0.9,0.95,0.975,0.99)){
      quant <- qtukey1df[paste(length(means)),paste(level)]/sqrt(2)
    } else {
      quant <- NaN
    }
  }
  width   <- quant*SE
  
  n   <- choose(length(means),2)
  out <- data.frame(Lower=numeric(n), Diff=numeric(n), Upper=numeric(n), SE=numeric(n), T=numeric(n), 'P(>t)'=numeric(n))
  rownames(out) <- paste(1:n)
  k <- 1
  for(i in 1:(length(means)-1)){
    for(j in (i+1):length(means)){
      mij <- means[i]-means[j]
      p <- ptukey(abs(mij)/(SE/sqrt(2)),length(means),df,lower.tail=FALSE)
      lower <- mij-width
      upper <- mij+width
      out[k,] <- c(lower, mij, upper, SE, mij/SE, p)
      rownames(out)[k] <- paste(mname[i],'-',mname[j],sep="")
      k <- k+1
    }
  }
  colnames(out) <- c("Lower", "Center", "Upper", "Std.Err", "t value","P(>t)")
  attr(out,"minSignDiff") <- width
  attr(out,"means") <- means
  attr(out,"level") <- level
  attr(out,"effVals") <- effVals
  attr(out,"resp") <- resp
  attr(out,"quant") <- quant
  attr(out,"error") <- SE
  attr(out,"df") <- df
  class(out) <- c("TukeyMix","data.frame")
  out
}


## CLD
cld.simple.glht <- function (object, alpha = 0.05, decreasing = TRUE, ...) {
  random <- ifelse(is.null(attr(object,'random')),FALSE,TRUE)
  
  if(attr(object,"type") != "Tukey"){
    class(object) <- class(object)[-1]
    ret <- cld(object)
    ret$object <- 0
    means <- sort(tapply(ret$y,ret$x,mean))
    if(decreasing){
      means <- rev(means)
      if(!is.null(dim(ret$mcletters$LetterMatrix))){
        p <- dim(ret$mcletters$LetterMatrix)[2]
        ret$mcletters$LetterMatrix <- ret$mcletters$LetterMatrix[,p:1]
      }
    }
    attr(ret$object,"means") <- means
    ret$lvl_order <- levels(ret$x)[order(means)]
  } else {
    # Catch failed Tukey
    if(is.nan(attributes(object$res)$quant)){
      warning('Tukey\'s HSD failed')
      ret <- 0
      attr(ret,"failed") <- TRUE
      class(ret) <- "cldMix"
      return(ret)
    }
    if(attr(object$res,"df") == 1){
      ret <- 0
      attr(ret,"df1") <- TRUE
      attr(ret,"failed") <- FALSE
      class(ret) <- "cldMix"
      return(ret)
    }
    object$test <- list()
    object$test$pvalues <- object$res[,6]
    class(object) <- class(object)[-1]
    ret <- cld(object)
    ret$object <- 0
    means <- sort(tapply(ret$y,ret$x,mean))
    if(decreasing){
      means <- rev(means)
      if(!is.null(dim(ret$mcletters$LetterMatrix))){
        p <- dim(ret$mcletters$LetterMatrix)[2]
        ret$mcletters$LetterMatrix <- ret$mcletters$LetterMatrix[,p:1]
      }
    }
    attr(ret$object,"means") <- means
    ret$lvl_order <- levels(ret$x)[order(means)]
    
  }
  class(ret) <- "cldMix"
  attr(ret, "alpha") <- alpha
  attr(ret,"failed") <- FALSE
  attr(ret,"df1")    <- FALSE
  ret
}

print.cldMix <- function(x, fill=TRUE, ...){
  # Graceful output if failed Tukey
  if(attr(x,"failed")){
    cat("No CLD available\n")
    return(invisible(NULL))
  }
  if(attr(x,"df1")){
    cat("No CLD available for single df error.\n")
    return(invisible(NULL))
  }
  
  object    <- x
  means     <- attr(object$object,"means")
  n         <- length(means)
  mname     <- names(means)
  letters   <- object$mcletters$LetterMatrix
  if(is.null(dim(letters))){
    lvl_order <- names(letters)
  } else {
    lvl_order <- rownames(letters) #object$lvl_order
  }
  morder    <- numeric(n)
  if(class(letters)=="logical"){
    letters <- as.matrix(letters)
  }
  G <- dim(letters)[2]
  LetterMatrix <- matrix("",nrow=dim(letters)[1],ncol=G)
  colnames(LetterMatrix) <- paste("G",1:G,sep="")
  
  for(i in 1:n){
    morder[i] <- which(lvl_order==mname[i])
    for(j in 1:G){
      if(letters[i,j])
        LetterMatrix[i,j] <- LETTERS[j]
    }
  }
  out <- data.frame("Mean"=means, LetterMatrix[morder,,drop=FALSE])
  if(fill){
    for(i in 2:dim(out)[2]){ # Loop over groups
      first <- match(LETTERS[i-1], out[,i])
      last  <- match(LETTERS[i-1], rev(out[,i]))
      out[first:(dim(out)[1]-last+1),i] <- LETTERS[i-1]
    }
  }
  cat("Tukey's HSD", paste("Alpha:", attr(object,"alpha")), "", sep="\n")
  print(out)
}


##########################
# CI for the grand mean of a linear model
CIgrandMean <- function(object, alpha=0.05){
  # Retrieve model
  opt <- options("contrasts")
  options(contrasts=c('contr.sum','contr.poly'))
  noRandom <- update(object)
  noRandom$random <- NULL
  model <- as.data.frame(Anova(noRandom, type='II', singular.ok=TRUE))
  if(!any("Mean Sq"%in%colnames(model))){
    model <- cbind(model[,1:2], model[,"Sum Sq"]/model[,"Df"], model[,3:dim(model)[2]])
    colnames(model)[3] <- "Mean Sq"
  }
  options(contrasts=opt$contrasts)
  
  n.eff <- dim(model)[1]-1
  model <- model[-(n.eff+1),]
  
  if(n.eff > 1){ # More than one model effect
    effects <- rownames(model)
    
    X <- matrix(0, nrow=n.eff, ncol=n.eff+1)
    X[,n.eff+1] <- 1
    for(i in 1:n.eff){
      X[i,i] <- 1
      splat <- strsplit(effects[i], ":")[[1]]
      if(length(splat) > 1){ # Interaction effect
        for(j in 1:length(splat)){
          X[i,match(splat[j], effects)] <- 1
        }
      }
    }
    C <- rref(X)[,n.eff+1]
  } else {
    X <- NULL
    C <- 1
  }
  MS <- model[,3]%*%C
  v  <- sum(model[,"Mean Sq"]*C)^2/sum((model[,"Mean Sq"]*C)^2/model[,"Df"])
  CI <- mean(model.response(model.frame(object))) + c(1,0,-1)*qt(alpha/2, v)*sqrt(MS/length(model.response(model.frame(object))))
  class(CI) <- "CIgm"
  attr(CI, "alpha") <- alpha
  CI
}
print.CIgm <- function(x, ...){
  cat(100*(1-attr(x,"alpha")), "% confidence interval for the grand mean (",x[2],")\n", sep="")
  cat("(", x[1], ", ", x[3], ")\n", sep="")
}
