step.gam <-
  function (object, scope, scale, direction = c("both", "backward", 
                                    "forward"), trace = TRUE, keep = NULL, steps = 1000, parallel=FALSE,...) 
{
 trace=as.numeric(trace)
 get.visit <- function(trial, visited){
    match(paste(trial,collapse=""),apply(visited,2,paste,collapse=""),FALSE)
  }
deviancelm <-
function(object, ...)
if(is.null(w <- object$weights)) sum(object$residuals^2) else sum(w * object$
		residuals^2)
 scope.char <- function(formula) {
    formula = update(formula, ~-1 + .)
    tt <- terms(formula)
    tl <- attr(tt, "term.labels")
    if (attr(tt, "intercept")) 
      c("1", tl)
    else tl
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
                         namc))
  }
  untangle.scope <- function(terms, regimens) {
    a <- attributes(terms)
    response <- deparse(a$variables[[2]])
    term.labels <- a$term.labels
    if (!is.null(a$offset)) {
      off1 <- deparse(a$variables[[a$offset]])
    }
    nt <- length(regimens)
    select <- integer(nt)
    for (i in seq(nt)) {
      j <- match(regimens[[i]], term.labels, 0)
      if (any(j)) {
        if (sum(j > 0) > 1) 
          stop(paste("The elements of a regimen", i, 
                     "appear more than once in the initial model", 
                     sep = " "))
        select[i] <- seq(j)[j > 0]
        term.labels <- term.labels[-sum(j)]
      }
      else {
        if (!(j <- match("1", regimens[[i]], 0))) 
          stop(paste("regimen", i, "does not appear in the initial model", 
                     sep = " "))
        select[i] <- j
      }
    }
    if (length(term.labels)) 
      term.labels <- paste(term.labels, "+")
    if (!is.null(a$offset)) 
      term.labels <- paste(off1, term.labels, sep = " + ")
    return(list(response = paste(response, term.labels, sep = " ~ "), 
                select = select))
  }
  make.step <- function(models, fit, scale, object) {
    chfrom <- sapply(models, "[[", "from")
    chfrom[chfrom == "1"] <- ""
    chto <- sapply(models, "[[", "to")
    chto[1]="<start>"
    chto[chto == "1"] <- ""
    dev <- sapply(models, "[[", "deviance")
    df <- sapply(models, "[[", "df.resid")
    ddev <- c(NA, diff(dev))
    ddf <- c(NA, diff(df))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                 "\nInitial Model:", deparse(as.vector(formula(object))), 
                 "\nFinal Model:", deparse(as.vector(formula(fit))), 
                 paste("\nScale: ", format(scale), "\n", sep = ""))
#    rowns=paste(chfrom,chto,sep=" -> ")
#    rowns[1]="<start>"
#    rowns=paste(seq(rowns)-1,rowns,sep=": ")
    aod <- data.frame(From=chfrom,To=chto, Df = ddf, 
                      Deviance = ddev, `Resid. Df` = df, `Resid. Dev` = dev, 
                      AIC = AIC, check.names = FALSE)
     aod <- as.anova(aod, heading)
     class(aod)=c("stepanova","data.frame")
     fit$anova=aod
    fit
  }
  direction <- match.arg(direction)
  if (missing(scope)) 
    stop("you must supply a scope argument to step.gam(); the gam.scope() function might be useful")
  if (!is.character(scope[[1]])) 
    scope <- lapply(scope, scope.char)
  response <- untangle.scope(object$terms, scope)
  form.y <- response$response
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  items <- response$select
  family <- family(object)
  Call <- object$call
  term.lengths <- sapply(scope, length)
  n.items <- length(items)
  visited <- matrix(items)
  form.vector <- character(n.items)
  for (i in seq(n.items)) form.vector[i] <- scope[[i]][items[i]]
  form <- deparse(object$formula)
  if (trace>0) 
    cat("Start: ", form)
  fit <- object
  n <- length(fit$fitted)
  if (missing(scale)) {
    famname <- family$family["name"]
    scale <- switch(famname, Poisson = 1, Binomial = 1, deviancelm(fit)/fit$df.resid)
  }
  else if (scale == 0) 
    scale <- deviancelm(fit)/fit$df.resid
  bAIC <- fit$aic
  if (trace>0) 
    cat("; AIC=", format(round(bAIC, 4)), "\n")
  models <- list(
                 list(deviance = deviance(fit), df.resid = fit$df.resid, 
                      AIC = bAIC, from = "", to = "")
                 )
  if (!is.null(keep))   {
    keep.list <- list(keep(fit,...))
    keep.it=TRUE}
  else keep.it=FALSE
  AIC <- bAIC + 1
  stepnum=0
  while (bAIC < AIC & steps > 0) {
    steps <- steps - 1
    stepnum=stepnum+1
    AIC <- bAIC
    form.list=NULL
###First some prelimenaries to see what formulas to try
    for (i in seq(n.items)) {
      if (backward) {
        trial <- items
        trial[i] <- trial[i] - 1
        if (trial[i] > 0 && !get.visit(trial,visited)) {
          visited<-cbind(visited,trial)
          tform.vector <- form.vector
          tform.vector[i] <- scope[[i]][trial[i]]
          form.list=c(form.list,list(list(trial=trial, form.vector=tform.vector, which=i)))
        }
      }
      if (forward) {
        trial <- items
        trial[i] <- trial[i] + 1
        if (trial[i] <= term.lengths[i] && !get.visit(trial,visited)){
          visited<-cbind(visited,trial)
          tform.vector <- form.vector
          tform.vector[i] <- scope[[i]][trial[i]]
          form.list=c(form.list,list(list(trial=trial, form.vector=tform.vector, which=i)))
        }
      }
    }
    if(is.null(form.list))break
### Now we are ready for the expensive loop
### Parallel is set up
#if(parallel&&require(foreach)){
if(parallel){
#   step.list=foreach(i=1:length(form.list),.inorder=FALSE,.packages="gam",.verbose=trace>1)%dopar%
   step.list=foreach(i=1:length(form.list),.inorder=FALSE,.verbose=trace>1)%dopar%
    {
      tform=paste(form.y, paste(form.list[[i]]$form.vector, collapse = " + "))
      update(object, eval(parse(text = tform)),trace = FALSE, ...)
    }
  }
### No parallel    
    else {
    step.list=as.list(sequence(length(form.list)))
    for(i in 1:length(form.list)){
      tform=paste(form.y, paste(form.list[[i]]$form.vector, collapse = " + "))
      step.list[[i]]=update(object, eval(parse(text = tform)),trace = FALSE, ...)
      if(trace>1)cat("Trial: ", tform,"; AIC=", format(round(step.list[[i]]$aic, 4)), "\n")
    }
}      
### end expensive loop
    taic.vec=sapply(step.list,"[[","aic")
    if(keep.it)  keep.list=c(keep.list, lapply(step.list,keep,...))
    bAIC=min(taic.vec)
    if (bAIC >= AIC | steps == 0) {
      if (keep.it) fit$keep <- re.arrange(keep.list)
      return(make.step(models, fit, scale, object))
    }
    else {
      o1=order(taic.vec)[1]
      fit=step.list[[o1]]
      form.list=form.list[[o1]]
      bwhich=form.list$which
      bfrom=form.vector[bwhich]
      form.vector=form.list$form.vector #this is the new one
      bto=form.vector[bwhich]
      if (trace>0) 
        cat(paste("Step:",stepnum,sep=""), deparse(fit$formula), "; AIC=", 
            format(round(bAIC, 4)), "\n")
      items <- form.list$trial
      models <- c(models,list(list(deviance = deviance(fit), df.resid = fit$df.resid, 
                           AIC = bAIC, from = bfrom, to = bto)))
    }
  }
}
