#######################################################
##           Compound Poisson GLMM                   ##
## Author: Wayne Zhang, actuary_zhang@hotmail.com    ##
#######################################################

cpglmm <- function(formula, link = "log", data, weights, offset,
                  subset, na.action, inits = NULL, 
                  contrasts = NULL, control = list(),
                  basisGenerators = c("tp", "bsp", "sp2d"),
                  optimizer = "nlminb", doFit = TRUE, nAGQ = 1) {
    
  call <- expand.call(match.call())  
  if (missing(data)) 
    data <- environment(formula)   
  link.power <- make.link.power(link)
  # identify smooth terms 
  formula <- eval(call$formula)
  tf <- terms.formula(formula, specials = eval(call$basisGenerators, 
        parent.frame(2)))
  n.f <- length(unlist(attr(tf, "specials")))
  # create model frame and get factor list  
  if (n.f) {
    call2 <- as.list(call)[-1]
    call2 <- call2[-match(c("link", "inits", "control", 
                            "optimizer", "doFit", "nAGQ"), names(call2), 0L)]
    setup <- do.call(frFL, as.list(call2))
    fr <- setup$m$fr 
    FL <- setup$m$FL
  } else {  
    fr <- lmerFrames(call, formula, contrasts)
    FL <- lmerFactorList(formula, fr, 0L, 0L)
  }
   
  # set control parameters
  ctr <- do.call(cplm.control, control)
  FL$dims["mxit"] <- ctr$max.iter
  FL$dims["mxfn"] <- ctr$max.fun
  # get dims 
  dm <- mkZt(FL, NULL)
  dm$dd["verb"] <- ctr$trace
  # update nAGQ in dd
  if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
  if (nAGQ %% 2 == 0) nAGQ <- nAGQ + 1L # reset nAGQ to be an odd number
  dm$dd["nAGQ"] <- as.integer(nAGQ)
  AGQlist <- .Call("cpglmm_ghq", nAGQ)
  M1 <- length(levels(dm$flist[[1]]))
  n <- ncol(dm$Zt)
  q <- dm$dd[["q"]]
  d <- dm$dd[["p"]]
  # default offset and prior wts
  if (is.null(fr$wts) || length(fr$wts) == 0)  
    fr$wts <- as.double(rep(1, n)) 
  if (is.null(fr$off) || length(fr$off) == 0)
    fr$off <- as.double(rep(0, n)) 
  if (M1 >= n) {
        msg1 <- "Number of levels of a grouping factor for the random effects\n"
        msg3 <- "n, the number of observations"
        if (dm$dd["useSc"]) 
            stop(msg1, "must be less than ", msg3)
        else if (M1 == n) 
            message(msg1, "is *equal* to ", msg3)
    }
  
  # initial values
  if (!is.null(inits)){    
    check.inits.cpglmm(inits, dm$dd['p'], dm$dd['nt'])
  } else {
    # generate starting values
    inits <- cpglm.init(fr, link.power)
  }
  names(inits$beta) <- names(fr$fixef)
  
  # input cpglmm class for optimization 
  # muEta and var are not initialized so that lmer is fitted when PQL.init is TRUE
  ans <- new(Class = "cpglmm", env = new.env( ), nlmodel = (~I(x))[[2]], 
             frame = fr$mf, call = call, flist = dm$flist, 
             Zt = dm$Zt, X = fr$X, y = as.numeric(fr$Y), 
             pWt = fr$wts, offset = fr$off,
             Gp = unname(dm$Gp), dims = dm$dd, 
             ST = dm$ST, A = dm$A, Cm = dm$Cm, L = dm$L, 
             Cx = rep(1.0, length((dm$A)@x)),  
             # was Cx = dm$A)@x which broke in R3.0.x. 
             # Referecen to the same value. Need duplicate in C level
             deviance = dm$dev, fixef = inits$beta, ranef = numeric(q), 
             u = numeric(q), eta = numeric(n), 
             mu = numeric(n), resid = numeric(n), 
             muEta = numeric(n), var = numeric(n),
             sqrtXWt = as.matrix(numeric(n)), sqrtrWt = numeric(n), 
             RZX = matrix(0, q, d), RX = matrix(0, d, d), 
             ghx = AGQlist[[1]], ghw = AGQlist[[2]], 
             p = inits$p, phi = inits$phi, link.power = as.double(link.power), 
             bound.p = ctr$bound.p, formula = formula, contrasts = contrasts,
             model.frame = fr$mf, inits = inits, vcov = matrix(0, d, d), smooths = list())
  # return cpglmm object if model fitting is turned off
  if (!doFit)  return(ans)
  
  # run optimization
  if (optimizer == "nlminb") {
    invisible(.Call("cpglmm_optimize", ans))
    if (ans@dims[["cvg"]] > 6) 
        warning(convergenceMessage(ans@dims[["cvg"]]))
  } else {    
    # function to be fed to optimizers (return deviance)
    # parm: theta + beta + log(phi) + p
    cpglmm_dev <- function(parm){             
      .Call("cpglmm_update_dev", ans, parm) 
    }
    # initial values for theta, beta, log(phi), p, 
    parm <- c(.Call("cpglmm_ST_getPars", ans), ans$fixef, log(ans$phi), ans$p) 
    parm <- unname(parm)
    # set bounds
    n.parm <- length(parm)
    lower <- rep(-Inf, n.parm)
    upper <- rep(Inf, n.parm)
    # reset bounds for diagonal elements of theta
    lower[1:sum(sapply(ans$ST, ncol))] <- 0
    # reset bounds for p
    lower[n.parm] <- ans$bound.p[1]
    upper[n.parm] <- ans$bound.p[2]

    # run optimizations
    rslt <- cplm_optim(parm, cpglmm_dev, lower = lower, upper = upper, 
                   control = ctr, optimizer = optimizer)
    ans@dims[["cvg"]] <- as.integer(rslt$convergence)
    if (rslt$convergence)  warning(rslt$message)
    
    # update ans using the found optima 
    invisible(.Call("cpglmm_update_dev", ans, rslt$par))    
  }
  # update random effects 
  invisible(.Call("cpglmm_update_ranef", ans))  
  # update sigmaML to be used in postVar
  dev <- ans@deviance 
  dev['sigmaML'] <- sqrt(ans@phi)
  ans@deviance <- dev 
  # update vcov in ans
  invisible(.Call("cpglmm_update_RX", ans))
  ans@vcov <- vcov(ans)
  # add smooth terms 
  if (n.f) {
    ans@smooths <- indsF(ans, setup$fct, setup$fctterm)
    vars <- unname(sapply(ans@smooths, function(tt) {
      pos <- grep("^x\\d?$", names(tt), perl = TRUE)
      sapply(pos, function(x) as.character(tt[[x]]))
      }))
    vars <- as.character(vars)
    ans@frame <- data.frame(ans@frame, base::subset(data, select = vars)) 
  }
  ans
  
}        












