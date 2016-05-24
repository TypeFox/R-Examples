

print.cem.att <- function(x,full=FALSE,...){
 if(!is.null(x$tab)){
  cat("\n")
  print(x$tab)
  }


  if(!is.null(x$att.model)){    
   if(x$extrapolate)
    cat(sprintf("\n%s with extrapolation:\n\n",x$mod.type))
   else
    cat(sprintf("\n%s on CEM matched data:\n\n",x$mod.type))

   if(full){
    print(x$att.model)
   } else {
    tmp <- x$att.model[,x$treatment]
    cat(sprintf("SATT point estimate: %f (p.value=%f)\n95%% conf. interval: [%f, %f]\n\n", 
     tmp["Estimate"], tmp["p-value"], tmp["Estimate"] + tmp["Std. Error"]*qnorm(0.025), 
	 tmp["Estimate"] + tmp["Std. Error"]*qnorm(0.975)))
   }
  }  
}



plot.cem.att <- function(x, obj, data, vars=NULL, plot=TRUE, ecolors,...){
    #require(lattice, quietly=TRUE)
 if(!is.null(x$TE)){
  s <-  x$att.model["Std. Error", obj$treatment]
#  q1 <- quantile(x$TE, 0.25)
#  q3 <- quantile(x$TE, 0.75)
  q1 <- -s
  q3 <- s
  if(missing(ecolors))
	 ecolors <- c("blue","black","red")
  level <- rep("zero", length(x$TE))
  cols <-  rep(ecolors[2], length(x$TE))
  level[which(x$TE <= q1)] <- "negative"
  level[which(x$TE >= q3)] <- "positive"
  cols[which(x$TE <= q1)] <- ecolors[1]
  cols[which(x$TE >= q3)] <- ecolors[3]
  ord <- order(x$TE)
   if(x$extrapolate)
    title <- sprintf("\n%s with extrapolation\n\n",x$mod.type)
   else
    title <- sprintf("\n%s on CEM matched data\n\n",x$mod.type)
   if(plot)
	 p1 <- dotplot(x$TE[ord], main=title,xlab="Treatment Effect",ylab="CEM Strata",col=cols[ord],scales=list(y=list(draw=FALSE)))
  }

  stID <- names(x$TE)
  use.vars <- obj$vars
  if(!is.null(vars))
   use.vars <- vars
  tmp <- NULL
  if( x$extrapolate ){
    for(st in stID){
     idx <- which(obj$strata == st)
     tmp <- rbind(tmp, sapply( use.vars, function(i) mean(as.numeric(data[idx,i]))))
    } 
    tmp <- as.data.frame(tmp)       
  } else {
    for(st in stID){
     idx <- which(obj$mstrata == st)
     tmp <- rbind(tmp, sapply( use.vars, function(i) mean(as.numeric(data[idx,i]))))
    }
    tmp <- as.data.frame(tmp)     
  }
  rownames(tmp) <- stID
  tmp$level <- factor(level,levels=c("negative","zero","positive"),ordered=TRUE)   
  p21 <- NULL
  p22 <- NULL
  p23 <- NULL
  n.low <- length(which(level=="negative"))
  n.zero <- length(which(level=="zero"))
  n.high <- length(which(level=="positive"))
  if(plot){
	if(n.low>0)
	  p21 <- parallelplot(~ tmp[use.vars], subset=level=="negative", tmp, main="negative",alpha=1.5/n.low, col="blue")
    if(n.zero>0)
      p22 <- parallelplot(~ tmp[use.vars], subset=level=="zero", tmp, main="zero",alpha=1.5/n.zero, col="black")
    if(n.high>0)
      p23 <- parallelplot(~ tmp[use.vars], subset=level=="positive",  tmp, main="positive",alpha=1.5/n.high, col="red")
    plot(p1, split=c(1,1,1,2))
    if(!is.null(p21))
      plot(p21, split=c(1,2,3,2), newpage=FALSE)
    if(!is.null(p22))
      plot(p22, split=c(2,2,3,2), newpage=FALSE)
	if(!is.null(p23))
      plot(p23, split=c(3,2,3,2), newpage=FALSE)
  }
  names(level) <- stID	
  match( obj$mstrata, stID ) -> idx
  tmp1 <- level[idx]
  names(tmp1) <- rownames(obj$X)	
  tmp2 <- x$TE[idx]
  names(tmp2) <- rownames(obj$X)	
  names(cols) <- stID
  tmp3 <- cols[idx]	
  names(tmp3) <- rownames(obj$X)		
  return(invisible(list(st.egroup=level, st.TE=x$TE, st.ecol=cols, 
						 ind.egroup=tmp1, ind.TE=tmp2, ind.ecol=tmp3)))	
}


############### THE LME STUFF ##############


`att` <-
function (obj, formula, data, model="linear", extrapolate=FALSE,ntree=2000) 
{

 mod.type <- NULL
 mod <- NULL

 if((class(obj)[1] != "cem.match") && (class(obj)[1] != "cem.match.list"))
  stop("Argument `obj' must be a `cem.match' or `cem.match.list' object") 

 mods <- c("linear", "lm", "lme", "linear-RE", "logit", "logistic", "rf", "forest")
 if(!(model %in% mods))
  stop(sprintf("model `%s' not yet implemented",model))

 if(model %in% c("lm", "linear")){
  mod <- "lm"
  mod.type <- "Linear regression model"
 }

 if(model %in% c("lme", "linear-RE")){
  mod <- "lme"
  mod.type <- "Linear random effect model"
 }

 if(model %in% c("forest", "rf")){
  mod <- "randomForest"
  mod.type <- "Random forest model"
 }

 if(model %in% c("logit", "logistic")){
  mod <- "glm"
  mod.type <- "Logistic model"
 }
 
  random.cem <- NULL
  random.all <- NULL
  ranef.cem <- NULL
  ranef.all <- NULL
  reg.cem <- NULL
  reg.all <- NULL
  rf.cem <- NULL

  TE <- NULL
  #if(require(nlme))
   control= lmeControl(maxIter=1000,opt="optim")

  if(class(obj)[1] == "cem.match"){

  if(mod=="randomForest"){
      #requireNamespace("randomForest", quietly=TRUE)
   out <- do.call(mod, list(formula=formula, data=data, subset= obj$matched==TRUE & obj$groups==obj$g.names[1],ntree=ntree))
   response <- all.vars(formula)[attr(terms(formula),"response")]
   tmp.data <- data
   tmp.data[,obj$treatment] <- 0
   prd <- predict(out, tmp.data)
   TEi <- data[,response]-prd

   if(extrapolate){
    g <- function(i){
     idt <- which(obj$strata == i & obj$groups==obj$g.names[2])
	 mean(TEi[idt])
     }
    strata <- na.omit(unique(obj$strata))
    TE <- sapply(strata, g)
    names(TE) <- strata
    ww <- table(obj$strata, obj$groups)[,2] # hard coded
    w.coef <- weighted.mean(TE, ww, na.rm=TRUE)
    idx <- which(is.na(TE))
    if(length(idx)>0)
     TE <- TE[-idx]   
    rf.cem <- matrix(NA, 4,1)

    v0 <- var(data[obj$groups==obj$g.names[2],response])
    v1 <- var(prd[obj$groups==obj$g.names[2]])
	
	} else {
     g <- function(i){
      idt <- which(obj$mstrata == i & obj$groups==obj$g.names[2])
	  mean(TEi[idt])
     }
    mstrata <- na.omit(unique(obj$mstrata))
    TE <- sapply(mstrata, g)
    names(TE) <- mstrata
    ww <- table(obj$mstrata, obj$groups)[,2] # hard coded
    w.coef <- weighted.mean(TE, ww, na.rm=TRUE)
    idx <- which(is.na(TE))
    if(length(idx)>0)
     TE <- TE[-idx]   
    rf.cem <- matrix(NA, 4,1)

    v0 <- var(data[obj$matched==TRUE & obj$groups==obj$g.names[2],response])
    v1 <- var(prd[obj$matched==TRUE & obj$groups==obj$g.names[2]])
	}

   
	dimnames(rf.cem) <-  list(c("Estimate", "Std. Error", "t value", "p-value"), obj$treatment)
   
	rf.cem["Estimate", ] <- w.coef
	rf.cem["Std. Error",] <-  sqrt(( v1+v0 ) * sum(ww^2)/sum(ww)^2)
    rf.cem["t value",] <-  rf.cem["Estimate",]/rf.cem["Std. Error",]
    rf.cem["p-value",] <-  2*(1-pnorm(rf.cem["t value",]))
    att.model <- rf.cem


	out <- list(att.model = att.model, tab=obj$tab, treatment=obj$treatment, extrapolate=extrapolate, mod.type=mod.type, TE=TE)
    class(out) <- "cem.att"
    return(out)

   } 

# random effects models
  if(mod=="lme"){
      #require(nlme, quietly=TRUE)
   data1 <- data
   data1$allID <- factor(obj$strata)
   data1$cemID <- factor(obj$mstrata)
   data1$matched <- obj$matched
   matched <- obj$matched
   f1 <- formula(sprintf(" ~ %s | allID",obj$treatment))   
   f2 <- formula(sprintf(" ~ %s | cemID",obj$treatment))   

   if(extrapolate){     
  # all data
    rand.try <- try(lme(update(formula, ~ . -cemID -allID - matched), data=data1, 
	  random = f1, keep.data=FALSE, control=control), silent=TRUE)
    if(class(rand.try) == "try-error"){
     cat("\nCannot estimate the random effects model for the complete data set\n")  
	} else { 
     rand.all <- rand.try
     random.all <- t(summary(rand.all)$tTable)	 
	 rand.cf <- coef(rand.all)
     ww <- table(obj$strata, obj$groups)[,2] # hard coded

     cf.idx <- match(rownames(rand.cf), names(ww))
     rownames(random.all) <- c("Estimate", "Std. Error", "DF", "t value", "p-value")
	 w.coef <- apply(rand.cf, 2, function(x) weighted.mean( x, w=ww[cf.idx]))
     random.all["Estimate", ] <- w.coef
	 random.all["Std. Error",] <-  random.all["Std. Error",]*sqrt(sum(ww^2)/sum(ww)^2)
     random.all["t value",] <-  random.all["Estimate",]/random.all["Std. Error",]
     random.all["p-value",] <-  2*(1-pnorm(random.all["t value",]))
	 att.model <- random.all
	 TE <- rand.cf[,obj$treatment]
	 names(TE) <- rownames(rand.cf)
	}
   } else {
# CEM restricted
    rand.try <- try(lme(update(formula, ~ . -cemID -allID - matched), data=data1, random = f2, 
	 subset=matched==TRUE, keep.data=FALSE, control=control), silent=TRUE)
    if(class(rand.try) == "try-error"){
     cat("\nCannot estimate the random effects model for the CEM matched subsample\n")   
    } else {
    rand.cem <- rand.try
	random.cem <- t(summary(rand.cem)$tTable)
	rand.cf <- coef(rand.cem)
    ww <- table(obj$mstrata, obj$groups)[,2] # hard coded
    cf.idx <- match(rownames(rand.cf), names(ww))
    rownames(random.cem) <- c("Estimate", "Std. Error", "DF", "t value", "p-value")
	w.coef <- apply(rand.cf, 2, function(x) weighted.mean( x, w=ww[cf.idx]))
    random.cem["Estimate", ] <- w.coef
	random.cem["Std. Error",] <-  random.cem["Std. Error",]*sqrt(sum(ww^2)/sum(ww)^2)
    random.cem["t value",] <-  random.cem["Estimate",]/random.cem["Std. Error",]
    random.cem["p-value",] <-  2*(1-pnorm(random.cem["t value",]))
    ranef.cem <- random.effects(rand.cem)
	att.model <- random.cem
	TE <- rand.cf[,obj$treatment]
	names(TE) <- rownames(rand.cf)
    }
   }

   out <- list(att.model = att.model, tab=obj$tab, treatment=obj$treatment, extrapolate=extrapolate, mod.type=mod.type, TE=TE)

   class(out) <- "cem.att"
   return(out)
  }

  if(mod=="lm"){
   out <- do.call(mod, list(formula=formula, data=data, weights=obj$w))
   fit <- fitted(out)
   response <- all.vars(formula)[attr(terms(formula),"response")]
   tmp.data <- data
   tmp.data[,obj$treatment] <- 0
   prd <- predict(out, tmp.data)
   TEi <- data[,response]-prd
   g <- function(i){
    idt <- which(obj$mstrata == i & obj$groups==obj$g.names[2])
	mean(TEi[idt])
   }
   
   mstrata <- na.omit(unique(obj$mstrata))
   TE <- sapply(mstrata, g)
   names(TE) <- mstrata
   idx <- which(is.na(TE))
   if(length(idx)>0)
    TE <- TE[-idx]   

   reg.cem <- t(summary(out)$coefficients)
   rownames(reg.cem) <-  c("Estimate", "Std. Error",  "t value", "p-value")
   att.model <- reg.cem


   if(extrapolate){
    idx2 <- which(!obj$matched & obj$groups==obj$g.names[2])
    f3 <- formula(sprintf(" ~  . -%s ",obj$treatment)) 
    f4 <- update(formula, f3)
	tmp.data <- data[idx2,]
	tmp.data[obj$treatment] <- 0
    y1 <- predict(out, newdata=tmp.data)
    response <- all.vars(f4)[attr(terms(formula),"response")]
    tmp <- t.test(y1, data[idx2, response])
    att4 <- diff(tmp$estimate)
    att4.serr <- abs(diff(tmp$estimate)/tmp$statistic)
    p <- obj$tab[2,2]/obj$tab[1,2]
    att4 <-  (1-p) * att4 + p * reg.cem["Estimate",obj$treatment]
    reg.all <- matrix(NA, 4,1)
    rownames(reg.all) <- rownames(reg.cem)
    colnames(reg.all) <- obj$treatment
    reg.all["Estimate",1] <- att4 
    reg.all["Std. Error",1] <- sqrt( p^2 * reg.cem["Std. Error", obj$treatment]^2 +
      (1-p)^2 * att4.serr^2  )
    reg.all["t value",1] <- att4/att4.serr 
    reg.all["p-value",1] <- 2*(1-pnorm(att4/att4.serr))
    reg.cem <- NULL
	att.model <- reg.all
    TE <- NULL

	
    fit <- fitted(out)
    tmp.data <- data
	tmp.data[obj$treatment] <- 0
    prd <- predict(out, tmp.data)
	TEi <- data[,response]-prd
	 
    g <- function(i){
     idt <- which(obj$strata == i & obj$groups==obj$g.names[2])
	 if(length(idt)>0)
	  return(mean(TEi[idt]))
	 return(NA)
	}
   
   strata <- na.omit(unique(obj$strata))
   TE <- sapply(strata, g)
   names(TE) <- strata
   idx <- which(is.na(TE))
   TE <- TE[-idx]

   } 
    
   out <- list(att.model = att.model, tab=obj$tab, treatment=obj$treatment, extrapolate=extrapolate, mod.type=mod.type, TE=TE)
#  out <- list(att.model = att.model, reg.cem = reg.cem, random.cem=random.cem, random.all=random.all, tab=obj$tab, 
#   ranef.cem=ranef.cem, ranef.all=ranef.all, cem.model = out, treatment=obj$treatment,
#   reg.all = reg.all,   extrapolate=extrapolate, mod.type=mod.type)
  class(out) <- "cem.att"
  return(out)
 }

  if(mod=="glm"){
   out <- do.call(mod, list(formula=formula, data=data, weights=obj$w,family="binomial"))
   fit <- fitted(out)
   response <- all.vars(formula)[attr(terms(formula),"response")]
   tmp.data <- data
   tmp.data[,obj$treatment] <- 0
   prd <- predict(out, tmp.data,type="response")
   TEi <- data[,response]-prd
   g <- function(i){
    idt <- which(obj$mstrata == i & obj$groups==obj$g.names[2])
	mean(TEi[idt])
   }
   
   mstrata <- na.omit(unique(obj$mstrata))
   TE <- sapply(mstrata, g)
   names(TE) <- mstrata
   idx <- which(is.na(TE))
   if(length(idx)>0)
    TE <- TE[-idx]   

   reg.cem <- t(summary(out)$coefficients)
   rownames(reg.cem) <-  c("Estimate", "Std. Error",  "t value", "p-value")
   att.model <- reg.cem


   if(extrapolate){
    idx2 <- which(!obj$matched & obj$groups==obj$g.names[2])
    f3 <- formula(sprintf(" ~  . -%s ",obj$treatment)) 
    f4 <- update(formula, f3)
	tmp.data <- data[idx2,]
	tmp.data[obj$treatment] <- 0
    y1 <- predict(out, newdata=tmp.data,type="response")
    response <- all.vars(f4)[attr(terms(formula),"response")]
    tmp <- t.test(y1, data[idx2, response])
    att4 <- diff(tmp$estimate)
    att4.serr <- abs(diff(tmp$estimate)/tmp$statistic)
    p <- obj$tab[2,2]/obj$tab[1,2]
    att4 <-  (1-p) * att4 + p * reg.cem["Estimate",obj$treatment]
    reg.all <- matrix(NA, 4,1)
    rownames(reg.all) <- rownames(reg.cem)
    colnames(reg.all) <- obj$treatment
    reg.all["Estimate",1] <- att4 
    reg.all["Std. Error",1] <- sqrt( p^2 * reg.cem["Std. Error", obj$treatment]^2 +
      (1-p)^2 * att4.serr^2  )
    reg.all["t value",1] <- att4/att4.serr 
    reg.all["p-value",1] <- 2*(1-pnorm(att4/att4.serr))
    reg.cem <- NULL
	att.model <- reg.all
    TE <- NULL

	
    tmp.data <- data
	tmp.data[obj$treatment] <- 0
    prd <- predict(out, tmp.data,type="response")
	TEi <- data[,response]-prd
	 
    g <- function(i){
     idt <- which(obj$strata == i & obj$groups==obj$g.names[2])
	 if(length(idt)>0)
	  return(mean(TEi[idt]))
	 return(NA)
	}
   
   strata <- na.omit(unique(obj$strata))
   TE <- sapply(strata, g)
   names(TE) <- strata
   idx <- which(is.na(TE))
   TE <- TE[-idx]

   } 
    
   out <- list(att.model = att.model, tab=obj$tab, treatment=obj$treatment, extrapolate=extrapolate, mod.type=mod.type, TE=TE)
  class(out) <- "cem.att"
  return(out)
 }


 }
 if(class(data) != "list")
  stop("Argument `data' must be a list of `data.frame's")

 n.cems <- length(obj) - 1

 if(length(data) != n.cems)
  stop("lengths of `cem.match.list' object and `data' do not match")
  
 est <- vector(n.cems, mode="list")
 for(i in 1:n.cems){ 
  out <- att(obj[[i]], formula=formula, data = data[[i]], model=model, extrapolate=extrapolate)
  est[[i]] <- out
 }

 ncoef <- NCOL(est[[1]])
 
 qoi <- numeric(ncoef)
 seq <- numeric(ncoef)

 for(i in 1:n.cems){
  att.m <- est[[i]]$att.model
  qoi <- qoi + att.m["Estimate",]
  seq <- seq + att.m["Std. Error",]^2
 }

 
 qoi <- qoi/n.cems
 seq <- seq/n.cems


 s2 <- numeric(ncoef)
 for(i in 1:n.cems){
  att.m <- est[[i]]$att.model
  s2 <- s2 + (att.m["Estimate",] - qoi)^2
 }
 s2 <- s2/(n.cems-1)

 S <- sqrt(seq + s2*(1+1/n.cems))

 TE <- vector(n.cems, mode="list")
 for(i in 1:n.cems)
  TE[[i]] <- est[[i]]$TE

 att.model <- rbind(qoi, S)
 att.model <- rbind(att.model, att.model[1,]/att.model[2,]) # t value
 att.model <- rbind(att.model, 2*(1-pnorm(att.model[3,]))) # p-value
 if(colnames(att.model)[1] == "qoi")
  colnames(att.model)[1] <- obj[[1]]$treatment
 rownames(att.model) <-  c("Estimate", "Std. Error",  "t value", "p-value")
 
 out <- list(mult=est, att.model = att.model, treatment=obj[[1]]$treatment, extrapolate=extrapolate, mod.type=mod.type, TE=TE) 
 class(out) <- "cem.att"
 out
}


summary.cem.att <- function(object, ...) {
	if(class(object) == "cem.att"){
		cat("\nTreatment effect estimation for data:\n\n")
		print(object$tab)
		
		if(object$extrapolate)
		cat(sprintf("\n%s estimated on all data\n\nCoefficients:\n",object$mod.type))  
		else
		cat(sprintf("\n%s estimated on matched data only\n\nCoefficients:\n",object$mod.type))  
		mod <- t(object$att.model) 
		printCoefmat(mod,has.Pvalue=TRUE, ...)
	}
}

