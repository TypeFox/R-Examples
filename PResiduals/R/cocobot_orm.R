cocobot.orm <- function(formula, data, link.x=c("logistic", "probit", "cauchit", "loglog", "cloglog"),
                        link.y=c("logistic", "probit", "cauchit", "loglog", "cloglog"),
                        subset, na.action=getOption('na.action'), 
                        fisher=FALSE,conf.int=0.95) {
  
  # Construct the model frames for x ~ z and y ~ z
  F1 <- Formula(formula)
  Fx <- formula(F1, lhs=1)
  Fy <- formula(F1, lhs=2)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  # We set xlev to a benign non-value in the call so that it won't get partially matched
  # to any variable in the formula. For instance a variable named 'x' could possibly get
  # bound to xlev, which is not what we want.
  mf$xlev <- integer(0) 
  mf[[1L]] <- as.name("model.frame")
  
  
  mx <- my <- mf
  
  # NOTE: we add the opposite variable to each model frame call so that
  # subsetting occurs correctly. Later we strip them off.
  mx[["formula"]] <- Fx
  yName <- all.vars(Fy[[2]])[1]
  mx[[yName]] <- Fy[[2]]
  
  my[["formula"]] <- Fy
  xName <- all.vars(Fx[[2]])[1]
  my[[xName]] <- Fx[[2]]
  
  mx <- eval(mx, parent.frame())
  mx[[paste('(',yName,')',sep='')]] <- NULL
  
  my <- eval(my, parent.frame())
  my[[paste('(',xName,')',sep='')]] <- NULL
  
  data.points <- nrow(mx)
  
  # Construct the model matrix z
  mxz <- model.matrix(attr(mx,'terms'),mx) 
  zzint <- match("(Intercept)", colnames(mxz), nomatch = 0L)
  if(zzint > 0L) {
    mxz <- mxz[, -zzint, drop = FALSE]
  }
  
  myz <- model.matrix(attr(my,'terms'),my) 
  zzint <- match("(Intercept)", colnames(myz), nomatch = 0L)
  if(zzint > 0L) {
    myz <- myz[, -zzint, drop = FALSE]
  }
  
  ## return value
  ans <- list(
    TS=list(),
    fisher=fisher,
    conf.int=conf.int,
    data.points=data.points
  )
  
  score.xz <- orm.scores(y=model.response(mx), X=mxz, link=link.x)
  score.yz <- orm.scores(y=model.response(my), X=myz, link=link.y)
  ts = corTS(score.xz$presid, score.yz$presid,
             score.xz$dl.dtheta, score.yz$dl.dtheta,
             as.matrix(score.xz$d2l.dtheta.dtheta), as.matrix(score.yz$d2l.dtheta.dtheta),
             score.xz$dpresid.dtheta, score.yz$dpresid.dtheta,fisher)
  ts.label = "PResid vs. PResid"
  
  ans$TS$TS <-  list( ts=ts$TS, var=ts$var.TS, pval=ts$pval.TS,
                      label = ts.label)
  ans <- structure(ans, class="cocobot")
  
  # Apply confidence intervals
  for (i in seq_len(length(ans$TS))){
    ts_ci <- getCI(ans$TS[[i]]$ts,ans$TS[[i]]$var,ans$fisher,conf.int)
    ans$TS[[i]]$lower <- ts_ci[1]
    ans$TS[[i]]$upper <- ts_ci[2]
  }
  
  ans
  
  
}

# ##### example
# ##### same result with cobot, if both X and Y are ordinal variables
# generate.data = function(alphax, betax, alphay, betay, eta, N) {
#   z = rnorm(N,0,1)
#   x = y = numeric(N)
#   
#   ## px is an N x length(alphax) matrix.
#   ## Each row has the TRUE cummulative probabilities for each subject.
#   px = (1 + exp(- outer(alphax, betax*z, "+"))) ^ (-1)
#   aa = runif(N)
#   for(i in 1:N)
#     x[i] = sum(aa[i] > px[,i])
#   x = as.numeric(as.factor(x))
#   ## x = x+1 may have category gaps if there are small probability categories.
#   
#   py = (1 + exp(- outer(alphay, betay*z+eta[x], "+"))) ^ (-1)
#   aa = runif(N)
#   for(i in 1:N)
#     y[i] = sum(aa[i] > py[,i])
#   y = as.numeric(as.factor(y))
#   ## y = y+1 may have category gaps if there are small probability categories.
#   
#   return(list(x=x, y=y, z=z))
# }


# N = 500
# alphay = c(-1, 0, 1)
# betay = -.5
# alphax = c(-1, 0, 1, 2)
# betax = 1
# eta = rep(0,5)
# data = generate.data(alphax, betax, alphay, betay, eta, N)
# cobot(as.factor(x)|as.factor(y) ~z, data=data)
# cocobot.orm(x|y~z, data=data)
# cobot(as.factor(x)|as.factor(y) ~z, data=data, link="probit")
# cocobot.orm(x|y~z, data=data, link.x="probit", link.y="probit") 
# 
# 
# ####### when X and Y are continuous
# set.seed(1)
# n <- 500
# z <- rnorm(n, 0, 1)
# x.latent <- rnorm(n, z+0.5, 1)
# x <-  exp(x.latent) + x.latent^3
# y.latent <- rnorm(n, x.latent+0.5*z)
# y <- exp(y.latent)
# #cobot(as.factor(x) |as.factor(y) ~z)
# Sys.time()
# cocobot.orm(x|y~z, link.x="logistic", link.y="logistic")
# Sys.time()
# cocobot.orm(x|y~z, link.x="probit", link.y="probit")
