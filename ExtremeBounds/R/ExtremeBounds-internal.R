# ExtremeBounds - Extreme Bounds Analysis in R
# Author: Marek Hlavac

# Note:
# v = focus variables (i.e., variables of interest)
# x = free variables (always included in the regression)
# z = doubtful variables ("cycled" through the regression in various combinations)
# k = how many doubtful variables should be included in the regression

.onAttach <- 
function(libname, pkgname) {
  packageStartupMessage("\nPlease cite as: \n")
  packageStartupMessage(" Hlavac, Marek (2015). ExtremeBounds: Extreme Bounds Analysis in R.")
  packageStartupMessage(" R package version 0.1.5.1. http://CRAN.R-project.org/package=ExtremeBounds \n")
}


.confidence.interval <-
function(beta.vector, se.vector, level) {
  if (is.null(names(beta.vector)) || is.null(names(se.vector))) { 
    out <- as.data.frame(cbind(beta.vector, se.vector))
  }
  else {
    out <- as.data.frame(cbind(beta.vector, se.vector[match(names(beta.vector), names(se.vector))]))
  }
  names(out) <- c("beta","se")
  
  # use the asymptotic normal distribution
  error <- qnorm((level+1)/2) * out$se
  out$ci.lower <- out$beta - error
  out$ci.upper <- out$beta + error
  
  return(out)
}

.ev <-
function(string.vector) {
  
  if (!is.character(string.vector)) { return(NULL) }
  if (is.null(string.vector)) { return(NULL) }
  
    
  out.vector <- c()
    
  for (j in 1:length(string.vector)) {
    string <- string.vector[j]
      
    s <- gsub(" ","",string, fixed=TRUE)
    s <- gsub("[^[:alnum:]._()]"," ",s,fixed=FALSE)
    s <- gsub("%in%","",s,fixed=TRUE)
    s <- strsplit(s, " ")[[1]]
      
    for (i in 1:length(s)) {
      first.bracket.pos <- regexpr("(", s[i], fixed=TRUE)
      if (first.bracket.pos != 0) { s[i] <- substr(s[i], first.bracket.pos+1, nchar(s[i])) }
    }
      
    s <- gsub("[()]","",s,fixed=FALSE)
    s <- s[is.na(suppressWarnings(as.numeric(s)))]
      
    out.vector <- c(out.vector, s)
  }
    
  return(out.vector)
}


################################## EBA ##################################
.eba.wrap <-
function(formula, data, y, free, doubtful, focus, k, mu, level, vif, exclusive, draws, reg.fun, se.fun, include.fun, weights, cl, ...) {
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
    if (!is.numeric(x)) { return(FALSE) }
    else { return(abs(x - round(x)) < tol) }
  }
  
  is.all.integers <-
  function(x) {
      if (!is.numeric(x)) { return(FALSE) }
      return (all(is.wholenumber(x))) 
  }
  
  is.all.integers.at.least <-
  function(x,at.least) {
    if (!is.na(match(NA,x))) { return(FALSE) }
    if (is.all.integers(x)) {
      return(min(x)>=at.least)
    }
    else {
      return(FALSE)
    }
  }
  
  # list of all k-member combinations of (z,v)
  get.combinations <-
  function(z, v, exclusive, k, draws=NULL) {
    
    message("\nGenerating combinations (2/4): ", appendLF=FALSE)
    
    max.k <- max(k)
    out <- matrix(NA, ncol=max.k, nrow=0)
    
    for (i in 1:length(k)) {
      combinations <- t(combn(z, k[i], simplify=T))
      
      combinations.resized <- matrix(NA, ncol=max.k, nrow=nrow(combinations))
      for (r in 1:nrow(combinations)) {
        for (c in 1:ncol(combinations)) {
          combinations.resized[r,c] <- combinations[r,c]
        }
      }
    
      out <- rbind(out, combinations.resized)
    }
    
    # only leave combinations that include the focus variable
    out <- out[apply(out, 1, FUN=function(x) (!all(is.na(match(v,x))))),]  
    
    # make sure that mutually exclusive variables are not included
    if (!is.null(exclusive)) {
      for (i in 1:length(exclusive)) {
        out <- out[apply(out, 1, FUN=function(x) (length(match(exclusive[[i]], x)[!is.na(match(exclusive[[i]],x))]))<=1),]
      }
    }
    
    if (is.null(nrow(out))) { out <- as.matrix(t(out)) }
    ncomb <- nrow(out)
    
    # randomly sample the requested number of draws
    if (!is.null(draws)) {
      
      if (draws > ncomb) { draws <- ncomb }
      which.rows <- sample(1:ncomb, size=draws, replace=FALSE)
      out <- out[which.rows,]
    }
    
    if (is.null(nrow(out))) { out <- as.matrix(t(out)) }
    nreg <- nrow(out)
    
    out.list <- list(out, ncomb, nreg)
    names(out.list) <- c("combinations","ncomb","nreg")
    
    if (nreg < ncomb) { 
      nreg.ratio = round(nreg/ncomb*100, 2)
      message("Estimate ",nreg," randomly sampled regressions out of ", ncomb," combinations (", nreg.ratio,"%).")
    }
    else {
      message("Estimate all ",nreg," combinations.")
    }
    
    return(out.list)
    
  }
  
  # returns highest value <= median, median itself, highest value >= median
  get.median <-
  function(x) {
    x <- x[!is.na(x)] # remove NAs
    if (length(x) != 0) {
      median <- median(x)
      median.lower <- max(x[x<=median]) 
      median.higher <- min(x[x>=median])
    }
    else {
      median <- median.lower <- median.higher <- NA
    }
  
    out <- c(median.lower, median, median.higher)
    names(out) <- c("median.lower","median","median.higher")
    return(out)
  }
  
  # mean wrapper that returns NA instead of NaN
  get.mean <-
  function(...) {
    value <- do.call(mean, list(...))
    if (!is.nan(value)) { return(value) }
    else { return(NA) }
  }
  
  # mean wrapper that returns NA instead of NaN
  get.weighted.mean <-
  function(...) {
    value <- do.call(weighted.mean, list(...))
    if (!is.nan(value)) { return(value) }
    else { return(NA) }
  }
  
  # population standard deviation
  stdev <-
  function(x) {
    if (length(x[!is.na(x)]) == 0) { return(NA) }
    
    N <- length(x)
    stdev <- sqrt((N-1)/N) * sd(x, na.rm=TRUE)
    return(stdev)
  }
  
  # figure out what R's variable labels will be + assign mu values
  variable.labels <-
  function(y, v, x, mu, mu.default, data) {
    
    out <- list()
    
    # put intercept first, and call it "free"
    out$name <- "1"
    out$label <- intercept.string
    out$type <- "free"
    if (!is.na(mu[intercept.string])) { out$mu <- mu[intercept.string] }
    else { out$mu <- mu.default }
    
    variable.source <- list(x, v)
    variable.types <- c("free","focus")
    
    # focus variables
    for (i in 1:length(variable.source)) {
      src <- variable.source[[i]]
      for (j in 1:length(src)) {
        if (src[j]!="1") {  # if not intercept (intercept has already been dealt with)
          reg.formula <- create.formula(y, src[j])
          reg.matrix <- model.matrix(reg.formula, data=data)
        
          lab <- colnames(reg.matrix)[colnames(reg.matrix)!=intercept.string]
        
          out$name <- c(out$name, rep(src[j], times=length(lab)))
          out$label <- c(out$label, lab)
          out$type <- c(out$type, rep(variable.types[i], times = length(lab)))
          
          # user-given mu cutoffs
          if (!is.na(mu[src[j]])) {
            out$mu <- c(out$mu, rep(mu[src[j]], times=length(lab)))
          }
          else {
            out$mu <- c(out$mu, rep(mu.default, times=length(lab)))
          }
          
        }
      } 
    }
    
    return(out)
  }
  
  get.vars <-
  function(lab) {
    
    # labels
    x.label <- lab$label[lab$type=="free"]
    v.label <- lab$label[lab$type=="focus"]
    
    # names
    x.name <- lab$name[lab$type=="free"]
    v.name <- lab$name[lab$type=="focus"]
    
    # types
    x.type <- lab$type[lab$type=="free"]
    v.type <- lab$type[lab$type=="focus"]
    
    # mu values
    x.mu <- lab$mu[lab$type=="free"]
    v.mu <- lab$mu[lab$type=="focus"]
    
    vars.labels <- c(x.label, v.label)
    vars.names <- c(x.name, v.name)
    vars.types <- c(x.type, v.type)
    vars.mu <- c(x.mu, v.mu)
    
    names(vars.labels) <- names(vars.names) <- names(vars.types) <- names(vars.mu) <- vars.labels
    
    out <- list(vars.labels, vars.names, vars.types, vars.mu)
    names(out) <- c("labels","names","types","mu")
    return(out)
  }
  
  create.formula <-
  function(lhs, rhs) {
    lhs.formula <- paste(lhs,"~ ",sep=" ")
    
    # remove the 1, if more than intercept in the function (aesthetic reasons)
    if ((length(rhs) > 1) && (rhs[1]=="1")) {
      rhs <- rhs[2:length(rhs)]
    }
    
    rhs.formula <- rhs[1]
    
    l <- length(rhs)
    if (l > 1) {
      for (i in 2:length(rhs)) {
        rhs.formula <- paste(rhs.formula, rhs[i], sep=" + ")
      }
    }
    
    string.formula <- paste(lhs.formula, rhs.formula, sep="")
    return( as.formula(string.formula ))
  }
  
  # include backticks
  backticks <- 
  function(s) {
    return( paste("`", s, "`", sep="") )
  }
  
  # remove backticks
  no.backticks <- 
  function(s) {
    s.out <- s
    if (substr(s.out,1,1)=="`") {
      s.out <- substr(s.out,2,nchar(s.out))
    }
      
    len <- nchar(s.out)
    if (substr(s.out,len,len)=="`") {
      s.out <- substr(s.out,1,len-1)
    }
      
    return(s.out)
  }
  
  # calculate the Variance Inflation Factor
  calculate.vif <-
  function(lm.object) {
    
    lm.vars <- colnames(model.matrix(lm.object))[-1]
    for (i in 1:length(lm.vars)) { lm.vars[i] <- backticks(lm.vars[i]) }
    
    temp.vif <- rep(NA, length=length(lm.vars))
    names(temp.vif) <- lm.vars
    
    for (i in 1:length(lm.vars)) {
      lhs <- lm.vars[i]
      rhs <- lm.vars[-i]

      if (length(rhs)==0) { rhs <- "1"}
      
      if (lhs != "1") {  # VIF only exists for non-intercept variables
        reg.beta.i.formula <- create.formula(lhs, rhs)
        reg.beta.i <- lm(reg.beta.i.formula, data=as.data.frame(model.matrix(lm.object)))
        rsq.beta.i <- summary(reg.beta.i)$r.squared
      
        temp.vif[lhs] <- 1/(1-rsq.beta.i)
      }
    }
    
    for (i in 1:length(temp.vif)) { names(temp.vif)[i] <- no.backticks(names(temp.vif)[i]) }
    
    return(temp.vif)
  }

  weights.lri <- 
  function(model.object, reg.fun, ...) {
    model.data <- model.frame(model.object)
    only.constant.formula <- create.formula(names(model.data)[1],"1")
    only.constant.object <- suppressMessages(do.call(reg.fun, list(formula=only.constant.formula, data=model.data, ...)))
      
    LL.o <- logLik(only.constant.object)
    LL.m <- logLik(model.object)
      
    LRI <- 1 - (LL.m / LL.o)
    return(LRI)
  }
  
  weights.r.squared <- 
  function(model.object) {
    if (!is.atomic(model.object)) {
      out <- summary(model.object)$r.squared
      if (!is.null(out)) { return(out) }
    }
    return(NA)
  }
  
  weights.adj.r.squared <- 
  function(model.object) {
    if (!is.atomic(model.object)) {
      out <- summary(model.object)$adj.r.squared
      if (!is.null(out)) { return(out) }
    }
    return(NA)
  }
  
  run.regressions <-
  function(y, v, x, combinations, data, level, vif.max, vars.labels, vars.names, vars.mu, reg.fun=lm, se.fun=NULL, include.fun=NULL, weights, indicate.quantiles, ...) {
    
    how.many.combinations <- nrow(combinations)
    how.many.vars.labels <- length(vars.labels)
    
    reg.formula <- list()
    beta <- se <- var <- t <- p <- ci.lower <- ci.upper <- vif <- nobs <- formula <- vif.satisfied <- include <- weight <- cdf.mu.generic <- matrix(NA, ncol=how.many.vars.labels, nrow=how.many.combinations)
    colnames(beta) <- colnames(se) <- colnames(var) <- colnames(t) <- colnames(p) <- colnames(ci.lower) <- colnames(ci.upper) <- colnames(vif) <- colnames(formula) <- colnames(vif.satisfied) <- colnames(include) <- colnames(weight) <- colnames(cdf.mu.generic) <- vars.labels
    
    # run all regressions
    message("\nEstimating regressions (3/4):")
    
    for (r in 1:how.many.combinations) {
      
      if (r %in% indicate.quantiles) {
        message(paste(as.character(r)," / ",as.character(how.many.combinations)," (",
                as.character(round((r/how.many.combinations)*100,2)),"%)", sep=""))
      }
      
      combinations.trim <- combinations[r,][!is.na(combinations[r,])]
      reg.formula[[r]] <- create.formula(y, c(x, combinations.trim))
      reg <- suppressMessages(do.call(reg.fun, list(formula=reg.formula[[r]], data=data, ...))) # regression object
      reg.summary <- summary(reg)        # regression summary object
      
      for (i in 1:how.many.vars.labels) {
        
        variable.label <- vars.labels[i]
        variable.mu <- vars.mu[i]
        
        if (variable.label %in% names(reg$coefficients)) {
          beta[r,i] <- reg.summary$coefficients[variable.label, 1]
          
          # adjust standard errors, if se.fun is non-NULL
          if (is.null(se.fun)) {
            se[r,i]  <- reg.summary$coefficients[variable.label, 2]
          }
          else if (is.function(se.fun)) { # user-given function for weights
            se.out <- suppressMessages(do.call(se.fun, list(reg)))
            se[r,i] <- se.out[variable.label]
          }
          else {
            se[r,i] <- NA
          }
          
          # variance is just the square of the standard error
          if (!is.na(se[r,i])) { var[r,i] <- (se[r,i])^2 } 
          else { var[r,i] <- NA }
          
          t[r,i]  <- reg.summary$coefficients[variable.label, 3]
          p[r,i]  <- reg.summary$coefficients[variable.label, 4]
          
          ci.temp <- suppressMessages(.confidence.interval(beta[r,i], se[r,i], level=level))
          ci.lower[r,i] <- ci.temp[variable.label,"ci.lower"]
          ci.upper[r,i] <- ci.temp[variable.label,"ci.upper"]
          
          vif[r,i] <- calculate.vif(reg)[variable.label] 
          nobs[r,i] <- length(reg$residuals)
          
          formula[r,i] <- Reduce(paste, deparse(reg.formula[[r]], width.cutoff = 500))
          
          # weights
          if (is.null(weights)) { weight[r,i] <- 1 }
          else if (is.function(weights)) { # user-given function for weights
            weight[r,i] <- suppressMessages(do.call(weights, list(reg)))
          }
          else {
            if (weights == "lri") { weight[r,i] <- weights.lri(reg, reg.fun, ...) }
            else if (weights == "r.squared") { weight[r,i] <- weights.r.squared(reg) }
            else if (weights == "adj.r.squared") { weight[r,i] <- weights.adj.r.squared(reg) }
            else { weight[r,i] <- NA }
          }
          
          # Sala-i-Martin, without assumption that betas are normally distributed across models 
          # get CDF(mu) for each model
          cdf.mu.generic[r,i] <- pnorm(variable.mu, mean = beta[r,i], sd = se[r,i], lower.tail=TRUE, log.p=FALSE)
          
          if (!is.null(vif.max)) { vif.satisfied[r,i] <- (vif[r,i] <= vif.max) }
          else { vif.satisfied[r,i] <- TRUE }
          
          if (!is.null(include.fun)) { include[r,i] <- as.logical(suppressMessages(do.call(include.fun, list(reg)))) }
          else { include[r,i] <- TRUE }
        }
        
      }
      
    }
    
    out <- list(beta, se, var, t, p, ci.lower, ci.upper, nobs, vif, vif.satisfied, formula, weight, cdf.mu.generic, include)
    names(out) <- c("beta","se","var","t","p","ci.lower","ci.upper", "nobs", "vif", "vif.satisfied", "formula", "weight", "cdf.mu.generic", "include")
    return(out)
  }
  
  
  eba.analysis <-
  function(coef, vars.labels, vars.types, vars.mu, r) {
    
    how.many.vars.labels <- length(vars.labels)
    
    # if both extreme bounds have the same sign --> robust, according to Leamer
    bounds.base <- as.data.frame(matrix(NA, nrow=how.many.vars.labels, ncol=6))
    colnames(bounds.base) <- c("type","leamer.lower","leamer.upper","leamer.robust","cdf.mu.normal","cdf.above.mu.normal")
    rownames(bounds.base) <- vars.labels
  
    
    bounds.base$type <- vars.types
    bounds.base$mu <- vars.mu
    
    # Leamer Analysis
    bounds.base$leamer.lower <- coef[["min.ci.lower"]]$ci.lower
    bounds.base$leamer.upper <- coef[["max.ci.upper"]]$ci.upper
    bounds.base$leamer.robust <- as.logical(( ((bounds.base$leamer.lower-bounds.base$mu) * (bounds.base$leamer.upper - bounds.base$mu)) > 0)) 
    
    ### Sala-i-Martin analysis
    
    # assuming betas are normally distributed across models
    mu.normal <- coef[["weighted.mean"]]$beta
    sigma.normal <- sqrt(coef[["weighted.mean"]]$var)  # standard deviation based on weighted mean of *variances*
    
    bounds.base$cdf.mu.normal <- pnorm(bounds.base$mu, mean=mu.normal, sd=sigma.normal, lower.tail=TRUE, log.p=FALSE)
    bounds.base$cdf.above.mu.normal <- 1 - bounds.base$cdf.mu.normal
    
    return(bounds.base)
  }
  
  run.eba <-
  function(y, v, x, combinations, data, level, vif.max, vars.labels, vars.names, vars.types, vars.mu, reg.fun, se.fun, include.fun, weights, cl, indicate.quantiles, ncomb, ...) {  
    
    r <- run.regressions(y, v, x, combinations, data, level, vif.max, vars.labels, vars.names, vars.mu, reg.fun, se.fun, include.fun, weights, indicate.quantiles,  ...)
    
    how.many.vars.labels <- length(vars.labels)
    
    stat.names <- c("beta","se","var","t","p","ci.lower","ci.upper","vif", "nobs", "formula", "weight", "cdf.mu.generic")
    points <- c("cdf.generic.unweighted","cdf.generic.weighted", "min", "max", "mean","weighted.mean","median","median.lower", "median.higher", "min.ci.lower", "max.ci.upper")
    how.many.points <- length(points)
    
    type <- beta <- se <- var <- t <- p <- ci.lower <- ci.upper <- vif <- nobs <- formula <- weight <- cdf.mu.generic <- matrix(NA, nrow=how.many.vars.labels, ncol=how.many.points)
    beta.significant <- beta.below.mu <- beta.above.mu  <- beta.significant.below.mu <- beta.significant.above.mu <- nreg.variable <- ncoef.variable <- rep(NA, times=how.many.vars.labels)
    
    rownames(type) <- rownames(beta) <- rownames(se) <- rownames(var) <- rownames(t) <- rownames(p) <- rownames(ci.lower) <- rownames(ci.upper) <- rownames(vif) <- rownames(nobs) <- rownames(formula) <- rownames(weight) <- rownames(cdf.mu.generic) <- vars.labels
    colnames(type) <- colnames(beta) <- colnames(se) <- colnames(var) <- colnames(t) <- colnames(p) <- colnames(ci.lower) <- colnames(ci.upper) <- colnames(vif) <- colnames(nobs) <- colnames(formula) <- colnames(weight) <- colnames(cdf.mu.generic) <- points
    names(beta.significant) <- names(beta.below.mu) <- names(beta.above.mu) <- names(beta.significant.below.mu) <- names(beta.significant.above.mu) <- names(nreg.variable)  <- names(ncoef.variable) <- vars.labels
    
    coef <- list()
    
    message("\nCalculating bounds (4/4): ", appendLF=FALSE)
    for (j in 1:how.many.points) {
      
      pt <- points[j]
      
      for (i in 1:how.many.vars.labels){
        
        variable.label <- vars.labels[i]
        variable.mu <- vars.mu[i]
        
        variable.type <- vars.types[i]
        type[variable.label, pt] <- variable.type
        
        r.adj <- list()
        r.unadj <- r[[stat.names[1]]][,i]
        
        for (q in 1:length(stat.names)) {
          replaced <- stat.names[q]
          
          # check for vif and intercept (no VIF exists for intercept) + check for include.fun
          if ((!is.null(vif.max)) && (vars.names[i]!="1")) { 
            replacement <- r[[replaced]][((r$vif.satisfied[,i]==TRUE) & (r$include[,i]==TRUE)),i] 
          }
          else { replacement <- r[[replaced]][(r$include[,i]==TRUE),i] }
            
          if (length(replacement) == 0) {
            replacement <- NA
          }
          
          r.adj[[replaced]] <- replacement
        }
        
        # get some info for Sala-i-Martin's EBA
        # - only one run necessary (hence j==1)
        if (j==1) {
          num.regs <- length(r.unadj[!is.na(r.unadj)])
          total.length <- length(r.adj$beta[!is.na(r.adj$beta)])
          
          number.not.significant <- length( r.adj$beta[(r.adj$ci.lower <= variable.mu) & (r.adj$ci.upper >= variable.mu) & (!is.na(r.adj$beta))] )
          beta.significant[i] <- (total.length - number.not.significant)/total.length
          
          if (is.nan(beta.significant[i])) { beta.significant[i] <- NA }
          
          beta.below.mu[i] <- length(r.adj$beta[(r.adj$beta < variable.mu) & (!is.na(r.adj$beta))]) / total.length
          beta.above.mu[i] <- length(r.adj$beta[(r.adj$beta > variable.mu) & (!is.na(r.adj$beta))]) / total.length
          
          if (is.nan(beta.below.mu[i])) { beta.below.mu[i] <- NA }
          if (is.nan(beta.above.mu[i])) { beta.above.mu[i] <- NA }
          
          beta.significant.below.mu[i] <- length(r.adj$beta[(r.adj$beta < variable.mu) & (r.adj$ci.upper < variable.mu) & (!is.na(r.adj$beta))]) / total.length
          beta.significant.above.mu[i] <- length(r.adj$beta[(r.adj$beta > variable.mu) & (r.adj$ci.lower > variable.mu) & (!is.na(r.adj$beta))]) / total.length
          
          if (is.nan(beta.significant.below.mu[i])) { beta.significant.below.mu[i] <- NA }
          if (is.nan(beta.significant.above.mu[i])) { beta.significant.above.mu[i] <- NA }
          
          nreg.variable[i] <- num.regs
          ncoef.variable[i] <- total.length
        }
  
        # get the minima, maxima, medians, CDF(0)
        if (pt == "min") { index <- which.min(r.adj$beta) }
        else if (pt == "max") { index <- which.max(r.adj$beta) }
        else if (pt == "median") { 
          value <- get.median(r.adj$beta)["median"]
          index <- which(r.adj$beta==value)[1] 
        }
        else if (pt == "median.lower")  {
          value <- get.median(r.adj$beta)["median.lower"]
          index <- which(r.adj$beta==value)[1] 
        }
        else if (pt == "median.higher") { 
          value <- get.median(r.adj$beta)["median.higher"]
          index <- which(r.adj$beta==value)[1] 
        }
        else if (pt == "min.ci.lower") { index <- which.min(r.adj$ci.lower) }
        else if (pt == "max.ci.upper") { index <- which.max(r.adj$ci.upper) }

        if (!(pt %in% c("mean","weighted.mean","cdf.generic.unweighted","cdf.generic.weighted"))) {
          if (length(index) == 0) {
            beta[variable.label, pt] <- NA
            se[variable.label, pt] <- NA
            var[variable.label, pt] <- NA
            t[variable.label, pt] <- NA
            p[variable.label, pt] <- NA
            ci.lower[variable.label, pt] <- NA
            ci.upper[variable.label, pt] <- NA
            vif[variable.label, pt] <- NA
            nobs[variable.label, pt] <- NA
            formula[variable.label, pt] <- NA
            weight[variable.label, pt] <- NA
            cdf.mu.generic[variable.label, pt] <- NA
          }
          else {
            beta[variable.label, pt] <- r.adj$beta[index]
            se[variable.label, pt] <- r.adj$se[index]
            var[variable.label, pt] <- r.adj$var[index]
            t[variable.label, pt] <- r.adj$t[index]
            p[variable.label, pt] <- r.adj$p[index]
            ci.lower[variable.label, pt] <- r.adj$ci.lower[index]
            ci.upper[variable.label, pt] <- r.adj$ci.upper[index]
            vif[variable.label, pt] <- r.adj$vif[index]
            nobs[variable.label, pt] <- r.adj$nobs[index]
            formula[variable.label, pt] <- r.adj$formula[index]
            weight[variable.label, pt] <- r.adj$weight[index]
            cdf.mu.generic[variable.label, pt] <- r.adj$cdf.mu.generic[index]
          }
        }
        # get means
        else if (pt == "mean") {  # unweighted mean
            beta[variable.label, pt] <- get.mean(r.adj$beta, na.rm=TRUE)
            se[variable.label, pt] <- get.mean(r.adj$se, na.rm=TRUE)
            var[variable.label, pt] <- get.mean(r.adj$var, na.rm=TRUE)
        }
        else if (pt == "weighted.mean") { # weighted mean
            beta[variable.label, pt] <- get.weighted.mean(x = r.adj$beta, w = r.adj$weight, na.rm=TRUE)
            se[variable.label, pt] <- get.weighted.mean(x = r.adj$se, w = r.adj$weight, na.rm=TRUE)
            var[variable.label, pt] <- get.weighted.mean(x = r.adj$var, w = r.adj$weight, na.rm=TRUE)
        }
        else if (pt == "cdf.generic.unweighted") {  # unweighted mean of CDF(mu)
            cdf.mu.generic[variable.label, pt] <- get.mean(r.adj$cdf.mu.generic, na.rm=TRUE)
        }
        else if (pt == "cdf.generic.weighted") {  # weighted mean of CDF(mu)
            cdf.mu.generic[variable.label, pt] <- get.weighted.mean(x = r.adj$cdf.mu.generic, w = r.adj$weight, na.rm=TRUE)
        }
      }
    
      
      # create data frame
      if (!(pt %in% c("mean","weighted.mean","cdf.generic.unweighted","cdf.generic.weighted"))) { 
        frame <- as.data.frame(matrix(NA,nrow=how.many.vars.labels, ncol=length(stat.names)))
        names(frame) <- stat.names 
      }
      else if (pt %in% c("mean", "weighted.mean")) { 
        frame <- as.data.frame(matrix(NA,nrow=how.many.vars.labels, ncol=4))
        names(frame) <- c("type","beta","se","var")
      }
      else if (pt %in% c("cdf.generic.unweighted","cdf.generic.weighted")) {
        frame <- as.data.frame(matrix(NA,nrow=how.many.vars.labels, ncol=3))
        names(frame) <- c("type","cdf.mu.generic", "cdf.above.mu.generic")
      }
      rownames(frame) <- vars.labels
      
      frame$type <- type[,pt]
      if (!(pt %in% c("cdf.generic.unweighted","cdf.generic.weighted"))) {
        frame$beta <- beta[,pt]
        frame$se <- se[,pt]
        frame$var <- var[,pt]
      }
      else {
        frame$cdf.mu.generic <- cdf.mu.generic[,pt]
        frame$cdf.above.mu.generic <- 1 - frame$cdf.mu.generic
      }
      
      if (!(pt %in% c("mean","weighted.mean","cdf.generic.unweighted","cdf.generic.weighted"))) {
        frame$t <- t[,pt]
        frame$p <- p[,pt]
        frame$ci.lower <- ci.lower[,pt]
        frame$ci.upper <- ci.upper[,pt]
        frame$vif <- vif[,pt]
        frame$nobs <- nobs[,pt]
        frame$formula <- formula[,pt]
        frame$weight <- weight[,pt]
        frame$cdf.mu.generic <- cdf.mu.generic[,pt]
      }
      
      coef[[pt]] <- frame
    }
    
    # run EBA: Leamer + Sala-i-Martin with normally distributed betas across models
    bounds <- eba.analysis(coef, vars.labels, vars.types, vars.mu, r)   # base for the bounds components
    bounds$beta.below.mu <- beta.below.mu
    bounds$beta.above.mu <- beta.above.mu
    bounds$beta.significant <- beta.significant
    bounds$beta.significant.below.mu <- beta.significant.below.mu
    bounds$beta.significant.above.mu <- beta.significant.above.mu
    bounds$cdf.mu.generic <- coef$cdf.generic.weighted$cdf.mu.generic
    bounds$cdf.above.mu.generic <- coef$cdf.generic.weighted$cdf.above.mu.generic
    
    bounds <- bounds[,c("type","mu",
                        "beta.below.mu","beta.above.mu",
                        "beta.significant","beta.significant.below.mu","beta.significant.above.mu",
                        "leamer.lower","leamer.upper","leamer.robust",
                        "cdf.mu.normal", "cdf.above.mu.normal",
                        "cdf.mu.generic","cdf.above.mu.generic")]
    
    # number of regressions estimated overall
    nreg <- nrow(r$beta)
        
    out <- list(bounds, cl, coef, vars.mu, level, ncomb, nreg, nreg.variable, ncoef.variable, r)
    names(out) <- c("bounds","call","coefficients","mu","level","ncomb","nreg","nreg.variable","ncoef.variable","regressions")
    class(out) <- "eba"
    
    message("Done.")
    return(out)
  } 
  
  # returns string without leading or trailing whitespace
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  
  # strsplit, but ignore anything inside brackets ( )
  str.split.outermost <-
  function(s, split.char) {
    
    inside.brackets <- 0
    all.splits <- c()
    current.split <- ""
    
    for (i in 1:nchar(s)) {
      c <- substr(s, i, i)  # current character
      if (c == "(") { inside.brackets <- inside.brackets + 1 }
      if (c == ")") { inside.brackets <- inside.brackets - 1 }
      
      if ((c == split.char) && (inside.brackets == 0)) {
        all.splits <- c(all.splits, current.split)
        current.split <- ""
      }
      else {
        current.split <- paste(current.split, c, sep = "")
      }
    }
    
    # add whatever remains
    if ((c != split.char) && (inside.brackets == 0)) {
      all.splits <- c(all.splits, current.split)
    }
    
    return(all.splits)
  }
  
  get.formula.term.labels <-
    function(Formula.object, rhs=NULL) {
      deparsed <- Reduce(paste, deparse(as.Formula(formula(as.Formula(Formula.object), lhs=0, rhs=rhs, width.cutoff=500))))
      subformula <- gsub("~","", deparsed, fixed=TRUE)
      return(trim(str.split.outermost(subformula, "+")))
    }
  
  ################################### 
  
  # set defaults
  intercept.string <- "(Intercept)"
  error.msg <- NULL
  mu.default <- 0
  
  # data frame
  if (!is.data.frame(data)) { error.msg <- c(error.msg,"Argument 'data' must contain a data frame.\n")}
  
  # there needs to be either a formula or y, free, doubtful, etc.
  if ((!is(formula, "formula")) && (!is.null(formula))) { error.msg <- c(error.msg,"Argument 'formula' must be NULL or a formula.\n")}
  
  if (is(formula, "formula")) {
      
    formula <- as.Formula(formula)
    
    y <- as.character(formula)[[2]]
    
    # if only one RHS of formula, assume everything is focus
    if (length(formula)[2] ==  1) {
      free <- NULL
      focus <- get.formula.term.labels(formula, rhs = 1)
      doubtful <- NULL
    }
    else if (length(formula)[2] ==  2) {
      free <- get.formula.term.labels(formula, rhs = 1)
      focus <- get.formula.term.labels(formula, rhs = 2)
      doubtful <- NULL
    } 
    else if (length(formula)[2] ==  3) {
      free <- get.formula.term.labels(formula, rhs = 1)
      focus <- get.formula.term.labels(formula, rhs = 2)
      doubtful <- get.formula.term.labels(formula, rhs = 3)
      doubtful <- unique(c(doubtful, focus))
    } 
    
    if (length(free) == 0) { free <- NULL }
    if (length(focus) == 0) { focus <- NULL }
    if (length(doubtful) == 0) { doubtful <- NULL }
      
  }
  

  # dependent variable
  if (!is.character(y)) { error.msg <- c(error.msg,"Argument 'y' must be a character string containing the name of the dependent variable.\n")}
  if (length(y)!=1)  { error.msg <- c(error.msg,"Argument 'y' must be of length 1.\n")}
  if (is.character(y) && (length(.ev(y))==1)) {
    if (!(.ev(y) %in% names(data)))  { error.msg <- c(error.msg,"Variable in argument 'y' must be in the data frame.\n")}
    if ((y %in% free) || (y %in% focus) || (y %in% doubtful)) { error.msg <- c(error.msg,"Variable in argument 'y' must not be among the free/focus/doubtful variables.\n")}
  }
  
  # free variables
  if ((!is.character(free)) && (!is.null(free))) { error.msg <- c(error.msg,"Argument 'free' must be NULL or a vector of character strings that contains the names of the free variables.\n")}
  if (is.character(free)) {
    if ((!is.null(free)) && (!all(.ev(free) %in% names(data)))) { error.msg <- c(error.msg,"All variables in argument 'free' must be in the data frame.\n")}
  }
  free <- unique(free)
  
  # doubtful variables
  if ((!is.character(doubtful)) && (!is.null(doubtful))) { error.msg <- c(error.msg,"Argument 'doubtful' must be NULL or a vector of character strings that contains the names of the doubtful variables.\n")}
  if (is.character(doubtful)) {
    if ((!is.null(doubtful)) && (!all(.ev(doubtful) %in% names(data)))) { error.msg <- c(error.msg,"All variables in argument 'doubtful' must be in the data frame.\n")}
  }
  doubtful <- unique(doubtful)
  
  # focus variables
  if ((!is.character(focus)) && (!is.null(focus))) { error.msg <- c(error.msg,"Argument 'focus' must be NULL or a vector of character strings that contains the names of the focus variables.\n")}
  if (is.character(focus)) {
    if ((!is.null(focus)) && (!all(.ev(focus) %in% names(data)))) { error.msg <- c(error.msg,"All variables in argument 'focus' must be in the data frame.\n")}
    if ((!is.null(focus)) && (!is.null(doubtful)) && (!all(focus %in% doubtful))) { error.msg <- c(error.msg,"The variables in argument 'focus' must be a subset of those in argument 'doubtful'.\n")}
  }
  focus <- unique(focus)
  
  # need at least one focus or doubtful variable
  if (is.null(focus) && (is.null(doubtful))) { error.msg <- c(error.msg,"At least one focus or doubtful variable must be specified.\n")}
  
  # if no focus variables are specified, assume interested in all doubtful variables, and vice versa
  if (is.null(focus)) { focus <- doubtful }
  if (is.null(doubtful)) { doubtful <- focus }
  
  # k = how many doubtful variables to include
  if (!is.all.integers.at.least(k,0))  { error.msg <- c(error.msg,"Argument 'k' must be a vector of non-negative integers.\n")}
  if ((length(doubtful)-1) < max(k)) { error.msg <- c(error.msg,"Argument 'k' is too high for the given number of doubtful variables.\n")}
  k <- unique(k)  # only unique values of k to prevent repetition
  
  # mu = named vector or a single number
  if (!is.numeric(mu)) { error.msg <- c(error.msg,"Argument 'mu' must be numeric.\n")}
  if (is.null(names(mu))) {  # if single non-named number, then this becomes mu-default
    if ((!is.vector(mu)) || (length(mu)!=1)) { error.msg <- c(error.msg,"Argument 'mu' must be a named vector or a single numeric value.\n")}
    else { mu.default <- mu }
  }
  else {
    if (!is.vector(mu))  { error.msg <- c(error.msg,"Argument 'mu' must be a named vector or a single numeric value.\n")}
  }
  
  # confidence level
  if (!is.numeric(level)) { error.msg <- c(error.msg, "Argument 'level' must be numeric.\n")}
  if (length(level)!=1)  { error.msg <- c(error.msg,"Argument 'level' must be of length 1.\n")}
  if (is.numeric(level) && (length(level)==1)) { 
    if (!((level >= 0) && (level<=1))) {
      error.msg <- c(error.msg, "Argument 'level' must be between 0 and 1.\n")
    }
  }
      
  # variance inflation factor
  if ((!is.numeric(vif)) && (!is.null(vif))) { error.msg <- c(error.msg, "Argument 'vif' must be NULL or numeric.\n")}
  if ((!is.null(vif)) && (length(vif)!=1))  { error.msg <- c(error.msg,"Argument 'vif' must be of length 1.\n")}
  
  # mutually exclusive variables
  if ((!is.list(exclusive)) && (!is(exclusive, "formula")) && (!is.null(exclusive)))  { error.msg <- c(error.msg, "Argument 'exclusive' must be NULL or a list of string vectors, or a Formula.\n") }
  
  if (is(exclusive, "formula")) {   # if given as a Formula, transform into list
    exclusive.formula <- as.Formula(exclusive)
    exclusive <- list()
    for (i in 1:(length(exclusive.formula)[2])) {
      exclusive[[i]] <- colnames(model.frame(exclusive.formula, lhs = 0, rhs = i, data = data))
    }
  }
  
  if (is.list(exclusive)) {
    
    error.exclusive.subset <- FALSE  # keeps track of which errors have already been announced
    error.exclusive.string <- FALSE
    
    for (i in 1:length(exclusive)) {
      if ((!is.character(exclusive[[i]])) && (!error.exclusive.string)) { 
        error.msg <- c(error.msg,"Argument 'exclusive' must contain string vectors (of variable names) as list components.\n")
        error.exclusive.string <- TRUE
      }
      if (is.character(exclusive[[i]])) {
        if ((!all(exclusive[[i]] %in% doubtful)) && (!error.exclusive.subset)) { 
          error.msg <- c(error.msg,"The variables in argument 'exclusive' must be a subset of those in argument 'doubtful'.\n")
          error.exclusive.subset <- TRUE
        }
      }
    }
  }
  
  # number of draws
  if ((!is.null(draws)) && (!is.all.integers.at.least(draws,1)))  { error.msg <- c(error.msg,"Argument 'draws' must be NULL or a positive integers.\n")}
  if ((!is.null(draws)) && (length(draws)!=1))  { error.msg <- c(error.msg,"Argument 'draws' must be of length 1.\n")}
   
  # regression function
  if (!is.function(reg.fun)) { error.msg <- c(error.msg,"Argument 'reg.fun' must be a function.\n")}
  
  # function for standard errors
  if ((!is.null(se.fun)) && (!is.function(se.fun))) { error.msg <- c(error.msg,"Argument 'se.fun' must be NULL or a function.\n")}
  
  # function for inclusion in the analysis
  if ((!is.null(include.fun)) && (!is.function(include.fun))) { error.msg <- c(error.msg,"Argument 'include.fun' must be NULL or a function.\n")}
  
  # weights
  if ((!is.null(weights)) && (!is.function(weights)) && (!(weights %in% c("lri","r.squared","adj.r.squared")))) { error.msg <- c(error.msg,"Argument 'weights' must be NULL, a function, or one of \"adj.r.squared\", \"lri\" or \"r.squared\".\n")}
  
  # stop execution if errors found
  if (!is.null(error.msg)) { 
    message(error.msg) 
    return(NULL)
  }
  
  # intercept is treated as a free variable
  if (!is.null(free)) { free <- c("1", free) }
  else { free <- "1" }
  
  ##################################
  
  message("ExtremeBounds: eba() performing analysis. Please wait.")
  
  # add 1 to k, because of the way the combination calculation is set up
  k <- k + 1
  
  # run things
  message("\nPreparing variables (1/4): ", appendLF=FALSE)
  lab <- variable.labels(y, focus, free, mu, mu.default, data)
  vars.labels <- get.vars(lab)$labels
  vars.names <- get.vars(lab)$names
  vars.types <- get.vars(lab)$types
  vars.mu <- get.vars(lab)$mu
  message("Done.")
  
  get.comb <- get.combinations(doubtful, focus, exclusive, k, draws)
  comb <- get.comb$combinations
  ncomb <- get.comb$ncomb
    
  indicate.quantiles <- quantile(1:nrow(comb), 
                                 c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
                                 type=1) # indicate progress at deciles
  
  out <- run.eba(y=y, v=focus, x=free, combinations=comb, data=data, 
                 level=level, vif.max=vif, vars.labels=vars.labels, 
                 vars.names=vars.names, vars.types=vars.types, vars.mu=vars.mu,
                 reg.fun=reg.fun, se.fun=se.fun, include.fun=include.fun, weights=weights, 
                 cl=cl, indicate.quantiles=indicate.quantiles, ncomb=ncomb, ...)  
  return(out)
}

################################## HIST.EBA ##################################
# histogram of EBA regression coefficients
.hist.eba.wrap <- function(x, variables, col, freq, main,
                           mu.show, mu.col, mu.lwd, mu.visible, 
                           density.show, density.col, density.lwd, density.args, 
                           normal.show, normal.col, normal.lwd, normal.weighted, 
                           xlim, ylim, cl, ...) {
  
  # what dimensions should we have for histogram display
  get.dimensions <- function(n) {
    square.root <- sqrt(n)
    ncols.dim <- ceiling(square.root)
    nrows.dim <- floor(square.root)
    if (n > (nrows.dim * ncols.dim)) { nrows.dim <- ncols.dim }
    
    out <- c(nrows.dim, ncols.dim)
    names(out) <- c("nrows","ncols")
    return(out)
  }
  
  # check that arguments are ok
  error.msg <- NULL
  if (!is.object(x)) { error.msg <- c(error.msg, "Argument 'x' must be an object of class \"eba\".\n")}
  else {
    if (class(x)!="eba") { error.msg <- c(error.msg, "Argument 'x' must be an object of class \"eba\".\n")}
  }
  
  if ((!is.character(variables)) && (!is.null(variables))) { error.msg <- c(error.msg,"Argument 'variables' must be NULL or a vector of character strings that contains variable names.\n")}
  if (is.character(variables)) {
    if ((!is.null(variables)) && (!all(variables %in% colnames(x$regressions$beta)))) { error.msg <- c(error.msg,"All variables in argument 'variables' must be in the model object.\n")}
  }
  
  if (!is.logical(freq)) { error.msg <- c(error.msg, "Argument 'freq' must be of type 'logical' (TRUE/FALSE).\n")}   
  if (length(freq)!=1) { error.msg <- c(error.msg, "Argument 'freq' must be of length 1.\n")}   
  
  if ((!is.character(main)) && (!is.null(main))) { error.msg <- c(error.msg,"Argument 'main' must be NULL or a vector of character strings that contains histogram titles.\n")}  
  
  if (!is.logical(mu.show)) { error.msg <- c(error.msg, "Argument 'mu.show' must be of type 'logical' (TRUE/FALSE).\n")}   
  if (length(mu.show)!=1) { error.msg <- c(error.msg, "Argument 'mu.show' must be of length 1.\n")}   
  
  if (!is.logical(mu.visible)) { error.msg <- c(error.msg, "Argument 'mu.visible' must be of type 'logical' (TRUE/FALSE).\n")}   
  if (length(mu.visible)!=1) { error.msg <- c(error.msg, "Argument 'mu.visible' must be of length 1.\n")}   
  
  if (!is.logical(density.show)) { error.msg <- c(error.msg, "Argument 'density.show' must be of type 'logical' (TRUE/FALSE).\n")}   
  if (length(density.show)!=1) { error.msg <- c(error.msg, "Argument 'density.show' must be of length 1.\n")}   
  
  if ((!is.list(density.args)) && (!is.null(density.args))) { error.msg <- c(error.msg,"Argument 'density.args' must be NULL or a list.\n")}
  
  if (!is.logical(normal.show)) { error.msg <- c(error.msg, "Argument 'normal.show' must be of type 'logical' (TRUE/FALSE).\n")}   
  if (length(normal.show)!=1) { error.msg <- c(error.msg, "Argument 'normal.show' must be of length 1.\n")}   
  
  if (!is.logical(normal.weighted)) { error.msg <- c(error.msg, "Argument 'normal.weighted' must be of type 'logical' (TRUE/FALSE).\n")}   
  if (length(normal.weighted)!=1) { error.msg <- c(error.msg, "Argument 'normal.weighted' must be of length 1.\n")}   
  
  if (!is.null(error.msg)) {
    message(error.msg) 
    return(NULL)
  }
  
  ##################################
  
  # if no variables are specified, all of the doubtful ones
  if (is.null(variables)) { variables <- colnames(x$regressions$beta) }
  
  # if no names in label vector, just assume that it labels variables from left to right
  if (!is.null(main)) {
    if (is.null(names(main))) {
      names(main) <- variables[1:length(main)]
    }
  } 
  
  ### draw the diagram
  how.many.variables <- length(variables)
  
  dimensions <- get.dimensions(how.many.variables)
  par(mfrow=c(dimensions["nrows"], dimensions["ncols"]), mar=c(2,2,2,2))
  
  histograms <- list()
  
  for (i in 1:how.many.variables) {
    var.label <- variables[i]
    print.var.label <- var.label
    if (!is.null(main)) {
      if (!is.na(main[var.label])) { print.var.label <- main[var.label] }
    }
      
    which.rows <- as.logical((x$regressions$vif.satisfied[,var.label]) & (x$regressions$include[,var.label]))
    x.values <- x$regressions$beta[which.rows,var.label]
    x.values <- x.values[!is.na(x.values)]
    
    mu.value <- x$mu[var.label]
    
    if (length(x.values[!is.na(x.values)]) == 0) { # in case there is nothing to plot
      histograms[[var.label]] <- h <- NULL
      names(plot(1, type="n", axes=F, xlab="", ylab="", main=print.var.label))
    }
    else {  # usual case when there is enough to plot
    
      # save histogram into a list
      histograms[[var.label]] <- h <- hist(x.values, plot=FALSE)
      
      # get default ylim that prevents cutting off
      if (freq == T) { ylim.max <- max(h$counts, na.rm=TRUE) }
      else { ylim.max <- max(h$density, na.rm=TRUE) }
      
      if (density.show) {
        density.object <- do.call(density, append(list(x=x.values), density.args))
        ylim.max <- max(max(density.object$y), ylim.max)
      }
      
      if (normal.show == TRUE) {
        if (normal.weighted == TRUE) {
          mu <- x$coefficients$weighted.mean[var.label,"beta"]
          sigma <- x$coefficients$weighted.mean[var.label,"se"]
        }
        else {
          mu <- x$coefficients$mean[var.label,"beta"]
          sigma <- x$coefficients$mean[var.label,"se"]
        }
        normal.max <- dnorm(mu, mean=mu, sd=sigma)
        ylim.max <- max(normal.max, ylim.max)
      }
      
      ylim.use <- c(0, ylim.max)
      if (!is.null(ylim)) { ylim.use <- ylim }
      
      ###
    
      if (mu.visible == TRUE) {
    
        min.x.value <- min(h$breaks, na.rm=TRUE)
        max.x.value <- max(h$breaks, na.rm=TRUE)
      
        if (mu.value <= min.x.value) { xlim.min <- mu.value } else { xlim.min <- min.x.value }
        if (mu.value >= max.x.value) { xlim.max <- mu.value } else { xlim.max <- max.x.value }
      
        xlim.use <- c(xlim.min, xlim.max)
        if (!is.null(xlim)) { xlim.use <- xlim }
      
        hist(x.values,
            col=col, freq = freq,
            ylab = "",
            xlim=xlim.use, ylim = ylim.use,
            main=print.var.label,...)
      } else {
        hist(x.values,
            col=col, freq = freq,
             ylab = "",
            ylim = ylim.use, main=print.var.label, ...)
      }
    
      if (mu.show == TRUE) { abline(v = mu.value, col=mu.col, lwd=mu.lwd) }
      if (density.show == TRUE) { lines(density.object, col=density.col, lwd=density.lwd) }
      if (normal.show == TRUE) { curve(dnorm(x, mean=mu, sd=sigma), col=normal.col, lwd=normal.lwd, add=TRUE) }
    }
  }
  
  out <- list(cl, histograms)
  names(out) <- c("call","histograms")
  class(out) <- "hist.eba"
  return(invisible(out))
}

################################## PRINT.EBA ##################################
.print.eba.wrap <- function(x, digits, ...) {
  
  # check that arguments are ok
  error.msg <- NULL
  if (!is.object(x)) { error.msg <- c(error.msg, "Argument 'x' must be an object of class \"eba\".\n")}
  else {
    if (class(x)!="eba") { error.msg <- c(error.msg, "Argument 'x' must be an object of class \"eba\".\n")}
  }
  
  if (!is.numeric(digits)) { error.msg <- c(error.msg, "Argument 'digits' must be of type 'numeric'.\n")}   
  if (length(digits)!=1) { error.msg <- c(error.msg, "Argument 'digits' must be of length 1.\n")}   
  
  if (!is.null(error.msg)) {
    message(error.msg) 
    return(NULL)
  }
  
  ##################################
  
  mu <- x$mu
  unique.mu <- unique(mu)
  
  beta.coefficients <- x$coefficients$weighted.mean[,c("type","beta","se")]
  beta.coefficients <- cbind(beta.coefficients, x$coefficients$min[,c("beta","se")])
  beta.coefficients <- cbind(beta.coefficients, x$coefficients$max[,c("beta","se")])
  beta.coefficients[,2:7] <- round(beta.coefficients[,2:7], digits)
  colnames(beta.coefficients) <- c("Type","Coef (Wgt Mean)","SE (Wgt Mean)","Min Coef","SE (Min Coef)","Max Coef","SE (Max Coef)")
  
  if (length(unique.mu) == 1) {
    beta.coefficients.distrib <- x$bounds[,c("type","beta.below.mu","beta.above.mu","beta.significant","beta.significant.below.mu","beta.significant.above.mu")]
    beta.coefficients.distrib[,2:6] <- round(100 * beta.coefficients.distrib[,2:6], digits) # multiply by 100 to get percentages
    colnames(beta.coefficients.distrib) <- c("Type",paste("Pct(beta < ", unique.mu,")", sep=""), paste("Pct(beta > ", unique.mu,")", sep=""), 
                                           paste("Pct(significant != ",unique.mu,")",sep=""),paste("Pct(signif & beta < ", unique.mu,")", sep=""), paste("Pct(signif & beta > ", unique.mu,")", sep=""))
    
    leamer.output <- x$bounds[,c("type","leamer.lower","leamer.upper","leamer.robust")]
    
    leamer.output$leamer.lower <- round(leamer.output$leamer.lower, digits)
    leamer.output$leamer.upper <- round(leamer.output$leamer.upper, digits)
    leamer.output$robust.fragile[leamer.output$leamer.robust==TRUE] <- "robust"
    leamer.output$robust.fragile[leamer.output$leamer.robust==FALSE] <- "fragile"
    leamer.output <- leamer.output[,-4]
    
    colnames(leamer.output) <- c("Type","Lower Extreme Bound","Upper Extreme Bound",paste("Robust/Fragile? (mu = ",unique.mu,")", sep=""))
    
    sala.i.martin.output <- x$bounds[,c("type","cdf.mu.normal","cdf.above.mu.normal", "cdf.mu.generic","cdf.above.mu.generic")]
    sala.i.martin.output[,2:5] <- round(100 * sala.i.martin.output[,2:5], digits)
    colnames(sala.i.martin.output) <- c("Type",
                                        paste("N: CDF(beta <= ",unique.mu,")", sep=""),paste("N: CDF(beta > ",unique.mu,")", sep=""),
                                        paste("G: CDF(beta <= ",unique.mu,")", sep=""),paste("G: CDF(beta > ",unique.mu,")", sep=""))
    
    
  }
  else {
    beta.coefficients.distrib <- x$bounds[,c("type","mu","beta.below.mu","beta.above.mu","beta.significant","beta.significant.below.mu","beta.significant.above.mu")]
    beta.coefficients.distrib[,3:7] <- round(100 * beta.coefficients.distrib[,3:7], digits) # multiply by 100 to get percentages
    colnames(beta.coefficients.distrib) <- c("Type","mu","Pct(beta < mu)", "Pct(beta > mu)","Pct(significant != mu)","Pct(signif & beta < mu)", "Pct(signif & beta > mu)")
    
    leamer.output <- x$bounds[,c("type","mu","leamer.lower","leamer.upper","leamer.robust")]
    
    leamer.output$leamer.lower <- round(leamer.output$leamer.lower, digits)
    leamer.output$leamer.upper <- round(leamer.output$leamer.upper, digits)
    leamer.output$robust.fragile[leamer.output$leamer.robust==TRUE] <- "robust"
    leamer.output$robust.fragile[leamer.output$leamer.robust==FALSE] <- "fragile"
    leamer.output <- leamer.output[,-5]
    
    colnames(leamer.output) <- c("Type","mu","Lower Extreme Bound","Upper Extreme Bound","Robust/Fragile?")
    
    sala.i.martin.output <- x$bounds[,c("type","mu","cdf.mu.normal","cdf.above.mu.normal", "cdf.mu.generic","cdf.above.mu.generic")]
    sala.i.martin.output[,3:6] <- round(100 * sala.i.martin.output[,3:6], digits)
    colnames(sala.i.martin.output) <- c("Type","mu",
                                        "N: CDF(beta <= mu)", "N: CDF(beta > mu)",
                                        "G: CDF(beta <= mu)", "G: CDF(beta > mu)")
  }
  
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Confidence level: ",x$level,"\n",sep="")
  cat("Number of combinations: ", x$ncomb,"\n",sep="")
  cat("Regressions estimated: ", x$nreg," (", round((x$nreg/x$ncomb*100),2),"% of combinations)\n\n",sep="")
  cat("Number of regressions by variable:\n\n")
  print.default(format(x$nreg.variable, digits = digits), print.gap = 1L, quote=FALSE, ...)
  
  cat("\nNumber of coefficients used by variable:\n\n")
  print.default(format(x$ncoef.variable, digits = digits), print.gap = 1L, quote=FALSE, ...)
  
  cat("\nBeta coefficients:\n\n")
  print.data.frame(beta.coefficients, row.names=TRUE,  print.gap=2L, ...)
  
  cat("\nDistribution of beta coefficients:\n\n")
  print.data.frame(beta.coefficients.distrib, row.names=TRUE,  print.gap=2L, ...)
  
  cat("\nLeamer's Extreme Bounds Analysis (EBA):\n\n")
  print.data.frame(leamer.output, row.names=TRUE,  print.gap=2L, ...)
  
  cat("\nSala-i-Martin's Extreme Bounds Analysis (EBA):\n")
  cat("- Normal model (N): beta coefficients assumed to be distributed normally across models\n")
  cat("- Generic model (G): no assumption about the distribution of beta coefficients across models\n\n")
  
  print.data.frame(sala.i.martin.output, row.names=TRUE, print.gap=2L, ...)
  
  cat("\n")
  invisible(x)
}