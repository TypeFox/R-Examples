SemiParSampleSel <- function(formula, data = list(), weights = NULL, subset = NULL, 
                             start.v = NULL, start.theta = NULL,
                             BivD = "N", margins = c("probit","N"), fp = FALSE, infl.fac = 1,  
                             rinit = 1, rmax = 100, iterlimsp = 50, pr.tolsp = 1e-6, bd = NULL,
                             parscale){ 
 
 #if(margins[2] == "G") stop("For tested version of Gamma case, check next release.")
 #if(margins[2] %in% c("P", "NB", "D", "PIG", "S")) stop("Check next release for final tested version of discrete case.")
 
  if(margins[2]=="GA") margins[2] <- "G" 
  if(margins[1]== "probit") margins[1] <- "N"  
 
  marg2  <- c("N", "G", "P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", "LG", 
            "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", 
            "ZALG", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG")
  marg2D <- c("P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", "LG", 
            "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", 
            "ZALG", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG")  
  copul  <- c("N", "FGM", "F", "AMH", "C0", "C90", "C180", "C270","J0", "J90", "J180", "J270","G0", "G90", "G180", "G270")
  
  
 
  if(!(margins[1] %in% c("N","probit"))) stop("Error in selection equation's margin name. It can currently only be probit or equivalently N.")
  if(!(margins[2] %in% marg2) ) stop("Error in outcome equation's margin name. It should be one of: N, GA, P, NB, D, PIG, S, BB, BI, GEOM, LG, NBII, WARING, YULE, ZIBB, ZABB, ZABI, ZIBI, ZALG, ZANBI, ZINBI, ZAP, ZIP, ZIP2, ZIPIG.")
  if(!(BivD %in% copul)) stop("Error in parameter BivD value. It should be one of: N, FGM, F, AMH, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270.")
   
  
  if(length(formula) > 2 && margins[2] %in% c("N", "G", "NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")){ if(length(formula)!=4) stop("You need to specify four equations.") } 
  if(length(formula) > 2 && margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")){ if(length(formula)!=3) stop("You need to specify three equations.") } 
  if(length(formula) > 2 && margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG")){ if(length(formula)!=5) stop("You need to specify five equations.") } 
    
  ########### 
   
  qu.mag <- sp <- nu <- gam2.1 <- KendTau <- gamlss.fit <- edf <- edf1 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- sp.theta <- NULL  
  sp4 <- gp4 <- gam4 <- X4 <- NULL  
  sp5 <- gp5 <- gam5 <- X5 <- NULL    
  
  X3.d2 <- X4.d2 <- X5.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- 0
  
  if(length(formula) > 2){
  
  f3t <- try(formula[[3]][[3]], silent = TRUE)  
  if(class(f3t)!="try-error") stop("The third equation does not require a response.")
  
  	if(length(formula) > 3){
  
  f4t <- try(formula[[4]][[3]], silent = TRUE)  
  if(class(f4t)!="try-error") stop("The fourth equation does not require a response.")  
  
  		if(length(formula) > 4){
  
  f5t <- try(formula[[5]][[3]], silent = TRUE)  
  if(class(f5t)!="try-error") stop("The fifth equation does not require a response.")    
  
  					}
  				}
  }
  
  ############
  # Data check 
  ############  
  
  ig <- interpret.gam(formula)
  mf <- match.call(expand.dots = FALSE)
  
  if( length(formula) == 2 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[2]]$response))
  if( length(formula) == 3 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[3]]$pred.names, ig[[2]]$response))
  if( length(formula) == 4 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[3]]$pred.names, ig[[4]]$pred.names, ig[[2]]$response))
  if( length(formula) == 5 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[3]]$pred.names, ig[[4]]$pred.names, ig[[5]]$pred.names, ig[[2]]$response))
  if( length(formula) == 6 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[3]]$pred.names, ig[[4]]$pred.names, ig[[5]]$pred.names, ig[[6]]$pred.names, ig[[2]]$response))
  
  fake.formula <- paste(ig[[1]]$response, "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(ig$fake.formula)
  mf$formula <- fake.formula 
  mf$start.v <- mf$method <- mf$start.theta <- mf$BivD <- mf$margins <- mf$infl.fac <- mf$fp <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$pr.tolsp <- mf$parscale <- mf$bd  <- NULL
  
  mf$drop.unused.levels <- TRUE 
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())

  indS <- as.logical(data[,ig[[1]]$response])==FALSE 
  indS <- ifelse( is.na(indS), FALSE, indS) 
  data[indS, ig[[2]]$response] <- ifelse( is.na(data[indS, ig[[2]]$response]), 0, data[indS, ig[[2]]$response]) 
  data <- na.omit(data)
   
  if(is.null(weights)) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
    
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
    
  #########################################
  # Naive models needed for starting values 
  #########################################



  gam1 <- eval(substitute(gam(formula.eq1, binomial(link="probit"), gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
  y1 <- gam1$y
  inde <- y1 > 0
  
  if(margins[2]=="N")               gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, data=data, weights=weights, subset=inde),                               list(weights=weights,inde=inde)))
  if(margins[2]=="G")               gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, data=data, weights=weights, subset=inde, family=Gamma(link = "log")),   list(weights=weights,inde=inde)))
  
  if(margins[2] %in% c("P", "NB", "D", "PIG", "S", "GEOM", "NBII", "WARING", "YULE", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG", "LG", "ZALG", "BB", "BI", "ZIBB", "ZABB", "ZABI", "ZIBI")) {
    gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, data=data, weights=weights, subset=inde, family=poisson(link = "log")), list(weights=weights,inde=inde)))
  }
  

  ##############
  # Data Objects
  ##############

  ##############
  # first two equations
  ##############
  
  n <- length(y1)
  X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]; 
  X2.d2 <- length(coef(gam2))
  X2 <- matrix(0, n, X2.d2, dimnames = list(c(1:n),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2)
  if (margins[2]=="LG") { y2 <- rep(1,n) } else { y2 <- rep(0,n) }; y2[inde] <- gam2$y
  l.sp1 <- length(gam1$smooth); if(l.sp1 != 0) sp1 <- gam1$sp else sp1 <- NULL 
  l.sp2 <- length(gam2$smooth); if(l.sp2 != 0) sp2 <- gam2$sp else sp2 <- NULL 
  gp1 <- gam1$nsdf
  gp2 <- gam2$nsdf
  
  #################  
  
  if(length(formula) == 3) theta <- rnorm(n); nad1 <- "theta"  
  
  if(length(formula) == 4){
                      sigma <- rnorm(n); nad1 <- "sigma"  
                      theta <- rnorm(n); nad2 <- "theta"  
  }
  
  if(length(formula) == 5){
                      sigma <- rnorm(n); nad1 <- "sigma"  
                      nu    <- rnorm(n); nad2 <- "nu"                        
                      theta <- rnorm(n); nad3 <- "theta"  
  }  
  
  #################
  
  if(length(formula) > 2){
  formula.eq3 <- formula[[3]] 
  formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") )   
  gam3 <- eval(substitute(gam(formula.eq3, weights=weights, data=data))) 
  X3 <- model.matrix(gam3)
  X3.d2 <- llM <- dim(X3)[2]
  l.sp3 <- length(gam3$sp)
  if(l.sp3 != 0) sp3 <- sp.theta <- gam3$sp else sp3 <- NULL 
  environment(gam3$formula) <- environment(gam2$formula)
  gp3 <- gam3$nsdf
  a.theta <- coef(gam3)
   
                        }
                        
  #################  
  
  if(length(formula) > 3){
  formula.eq4 <- formula[[4]] 
  formula.eq4 <- as.formula( paste(nad2,"~",formula.eq4[2],sep="") )   
  gam4 <- eval(substitute(gam(formula.eq4, weights=weights, data=data))) 
  X4 <- model.matrix(gam4)
  X4.d2 <- llM <- dim(X4)[2]
  l.sp4 <- length(gam4$sp)
  if(l.sp4 != 0) sp4 <- sp.theta <- gam4$sp else sp4 <- NULL 
  environment(gam4$formula) <- environment(gam2$formula)
  gp4 <- gam4$nsdf
  a.theta <- coef(gam4)
                        }   
                        
  #################  
  
  if(length(formula) > 4){
  formula.eq5 <- formula[[5]]  
  formula.eq5 <- as.formula( paste(nad3,"~",formula.eq5[2],sep="") )   
  gam5 <- eval(substitute(gam(formula.eq5, weights=weights, data=data))) 
  X5 <- model.matrix(gam5)
  X5.d2 <- llM <- dim(X5)[2]
  l.sp5 <- length(gam5$sp)
  if(l.sp5 != 0) sp5 <- sp.theta <- gam5$sp else sp5 <- NULL 
  environment(gam5$formula) <- environment(gam2$formula)
  gp5 <- gam5$nsdf
  a.theta <- coef(gam5)
                        }   
                        
  #################    
  
  
  
  
  
  if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0) && fp==FALSE){ 
                         qu.mag  <- S.m(gam1, gam2, gam3, gam4, gam5, l.sp1, l.sp2, l.sp3, l.sp4, l.sp5, margins, eq1 = "yes")
                         sp  <- c(sp1, sp2, sp3, sp4, sp5)
                         qu.mag1 <- S.m(gam1, gam2, gam3, gam4, gam5, l.sp1, l.sp2, l.sp3, l.sp4, l.sp5, margins, eq1 = "no")                                                                       
                         sp1 <- c(sp2, sp3, sp4, sp5) 
                         
                                                                              }

  #################


  if(length(formula) < 3){

  ########################################################
  # Approximate Heckmann's procedure for linear predictors # this will have to be removed eventually
  ########################################################

  p.g1 <- predict(gam1)
  imr <- data$imr <- dnorm(p.g1)/pnorm(p.g1)
  data$bd <- bd
  
  formula.eq2.1 <- update.formula(formula.eq2, ~. + imr)
  
  
  if (margins[2] %in% c("P", "NB", "D", "PIG", "S", "NBII", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG")) 
        { formula.eq2.1 <- update.formula(formula.eq2, round((. + mean(.[.>=0]))/2) ~. + imr)  }
  else if (margins[2] %in% c("BB", "BI")) 
        { formula.eq2.1 <- update.formula(formula.eq2,   cbind(round(bd*(. + 0.5)/(bd + 1)), bd - round(bd*(. + 0.5)/(bd + 1))) ~. + imr)      }
  else if (margins[2] %in% c("GEOM", "YULE", "WARING", "LG", "ZALG")) 
        { formula.eq2.1 <- update.formula(formula.eq2, ~. + imr)  }
  else if (margins[2] %in% c("ZABI", "ZIBI", "ZIBB", "ZABB")) 
        { formula.eq2.1 <- update.formula(formula.eq2, cbind(., bd - .) ~. + imr)  }
  
  
  if(margins[2]=="N")               gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=infl.fac, data=data, weights=weights, subset=inde),                               list(weights=weights,inde=inde)))
  if(margins[2]=="G")               gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=infl.fac, data=data, weights=weights, subset=inde, family=Gamma(link = "log")),   list(weights=weights,inde=inde)))
  
  if(margins[2] %in% c("P", "NB", "D", "PIG", "S", "GEOM", "NBII", "WARING", "YULE", "ZANBI", "ZINBI", "ZAP", "ZIP", "ZIP2", "ZIPIG", "LG", "ZALG")) { 
    gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=infl.fac, data=data, weights=weights, subset=inde, family=poisson(link = "log")), list(weights=weights,inde=inde)))
  }
  
  if(margins[2] %in% c("BB", "BI", "ZIBB", "ZABB", "ZABI", "ZIBI")) { 
    gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=infl.fac, data=data, weights=weights, subset=inde, family=binomial(link = "logit")), list(weights=weights,inde=inde)))
  }
  
  
  sigma <- sqrt(mean(residuals(gam2.1, type = "deviance")^2)+mean(imr[inde]*(imr[inde]+p.g1[inde]))*gam2.1$coef["imr"]^2)[[1]]
  names(sigma) <- "sigma"
  co  <- (gam2.1$coef["imr"]/sigma)[[1]] # this is an approximate correlation

  if(margins[2]=="G") { sigma <- k <- (summary(gam2)$dispersion)^(-1); names(k) <- "sigma" }
  
  #################
  # Starting values
  #################

  a.theta <- st.theta.star(start.theta, co, BivD); names(a.theta) <- "theta.star"

  if(margins[2] %in% c("N", "NB", "PIG", "WARING", "BB", "NBII")) { 
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],log(sigma),a.theta)
        names(start.v)[length(start.v)-1] <- "log.sigma"
  }
  
  if(margins[2] %in% c("ZIP", "ZAP", "ZIP2", "ZALG", "ZIBI", "ZABI")) { 
        sigma.star.start1 <- sum(gam2.1$y==0)/length(gam2.1$y)
        sigma.star.start1 <- ifelse(sigma.star.start1 > 0.999, 0.999, sigma.star.start1)
        sigma.star.start1 <- ifelse(sigma.star.start1 < 0.001, 0.001, sigma.star.start1)
        sigma.star.start1 <- qlogis(sigma.star.start1)
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], sigma.star.start1, a.theta)
        names(start.v)[length(start.v)-1] <- "logi.sigma"
  }
  
  if(margins[2]=="G") {
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],log(k), a.theta)
        names(start.v)[length(start.v)-1] <- "log.shape"
  }
  
  
  if(margins[2] %in% c("P", "YULE", "GEOM", "BI", "LG")) {
          if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],a.theta)
  }
  
  
  if(margins[2] %in% c("D")) { 
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], log(sigma), qlogis(0.51), a.theta)
        names(start.v)[length(start.v)-2] <- "log.sigma"
        names(start.v)[length(start.v)-1] <- "logi.nu" 
  }
  
  if(margins[2] %in% c("ZINBI")) { 
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], log(sigma), qlogis(((sum(gam2$y == 0)/length(gam2$y)) + 0.01)/2), a.theta)
        names(start.v)[length(start.v)-2] <- "log.sigma"
        names(start.v)[length(start.v)-1] <- "logi.nu" 
  }
  
  if(margins[2] %in% c("ZANBI")) { 
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], log(sigma), qlogis(max(sum(gam2$y == 0)/length(gam2$y), 0.01)), a.theta)
        names(start.v)[length(start.v)-2] <- "log.sigma"
        names(start.v)[length(start.v)-1] <- "logi.nu" 
  }
  
  if(margins[2] %in% c("ZIPIG")) {
        nu.star.start1 <- sum(gam2.1$y==0)/length(gam2.1$y)
        nu.star.start1 <- ifelse(nu.star.start1 > 0.999, 0.999, nu.star.start1)
        nu.star.start1 <- ifelse(nu.star.start1 < 0.001, 0.001, nu.star.start1)
        nu.star.start1 <- qlogis(nu.star.start1)
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], log(sigma), nu.star.start1, a.theta)
        names(start.v)[length(start.v)-2] <- "log.sigma"
        names(start.v)[length(start.v)-1] <- "logi.nu" 
  }
  
  if(margins[2] %in% c("ZIBB", "ZABB")) {
        nu.star.start1 <- sum(gam2.1$y==0)/length(gam2.1$y)
        nu.star.start1 <- ifelse(nu.star.start1 > 0.999, 0.999, nu.star.start1)
        nu.star.start1 <- ifelse(nu.star.start1 < 0.001, 0.001, nu.star.start1)
        nu.star.start1 <- qlogis(nu.star.start1)
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], log(sigma), nu.star.start1, a.theta)
        names(start.v)[length(start.v)-2] <- "log.sigma"
        names(start.v)[length(start.v)-1] <- "logi.nu" 
  }
  
    if(margins[2] == "S") { 
          if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], log(sigma), -0.5 , a.theta)
          names(start.v)[length(start.v)-2] <- "log.sigma"
          names(start.v)[length(start.v)-1] <- "nu"  
    }


start.v1 <- start.v[-c(1:X1.d2, length(start.v))]    
    
 }else {
 
 start.v <- start.v <- c( coef(gam1), coef(gam2), coef(gam3), coef(gam4), coef(gam5) )   
 start.v1 <- start.v[-c(1:X1.d2, (length(start.v) - (llM - 1)):length(start.v))]  
 
 }
  

  ###################
  # Fitting procedure
  ###################
  
    if(missing(parscale)) parscale <- 1   
  
    if(margins[2] %in% c("N", "G") ) funcop <- ghss else funcop <- ghssD  
    
    respvec <- cbind(y1, y2)
    
    
    VC <- list(X1 = X1, 
               X2 = X2,
               X3 = X3,
               X4 = X4,
               X5 = X5,
               X1.d2 = X1.d2, 
               X2.d2 = X2.d2,
               X3.d2 = X3.d2,
               X4.d2 = X4.d2,
               X5.d2 = X5.d2,
               gp1 = gp1, 
               gp2 = gp2,
               gp3 = gp3,
               gp4 = gp4,
               gp5 = gp5,
               n = n,
               bd = bd,
               l.sp1 = l.sp1, 
               l.sp2 = l.sp2,
               l.sp3 = l.sp3,
               l.sp4 = l.sp4,
               l.sp5 = l.sp5,
               infl.fac = infl.fac,
               weights = weights,
               fp = fp,
               BivD = BivD, 
               margins = margins, parscale = parscale)
  
    
    if(margins[2] %in% marg2D  ){ 
    
    gamlss.fit <- fit.SemiParSampleSel(funcop = ghssDuniv, start.v = start.v1, dat = respvec,
                           rinit = rinit, rmax = rmax, VC = VC, qu.mag = qu.mag1, sp = sp1,
                           iterlimsp = iterlimsp, pr.tolsp = pr.tolsp, univariate = TRUE)
                           
    start.v <- c(coef(gam1), gamlss.fit$fit$argument, a.theta)    
    sp <- c(sp1, gamlss.fit$sp, sp.theta)
    
    gam2$coefficients <- gamlss.fit$fit$argument[1:X2.d2]
    gam2$Vp <- gamlss.fit$magpp$Vb[1:X2.d2,1:X2.d2]
    gam2$sig2 <- 1
    gam2$edf <- gamlss.fit$magpp$edf[1:X2.d2]
                                        
                                        
    if(length(formula) > 3){
            
            if(margins[2] %in% c("NB", "PIG", "D", "S", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG")){ 
                                        gam3$coefficients <- gamlss.fit$fit$argument[X2.d2 + (1:X3.d2)]
	                                gam3$Vp <- gamlss.fit$magpp$Vb[X2.d2 + (1:X3.d2), X2.d2 + (1:X3.d2)]
	                                gam3$sig2 <- 1
                                        gam3$edf <- gamlss.fit$magpp$edf[X2.d2 + (1:X3.d2)]  
                                                        }
                                                        
            if(margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG")){ 
                                        gam4$coefficients <- gamlss.fit$fit$argument[X2.d2 + X3.d2 + (1:X4.d2)]
	                                gam4$Vp <- gamlss.fit$magpp$Vb[X2.d2 + X3.d2 + (1:X4.d2), X2.d2 + X3.d2 + (1:X4.d2)]
	                                gam4$sig2 <- 1
                                        gam4$edf <- gamlss.fit$magpp$edf[X2.d2 + X3.d2 + (1:X4.d2)]  
                                           }                                                        
                                               
                           } 
    
  } 
  

   fit <- fit.SemiParSampleSel(funcop = funcop, start.v = start.v, dat = respvec, rinit = rinit, rmax = rmax, 
                    VC = VC, qu.mag = qu.mag, sp = sp,
                    iterlimsp = iterlimsp, pr.tolsp = pr.tolsp, univariate = FALSE)


  ############################
  # Post-estimation quantities # this will have to be a separate function eventually
  ############################

  epsilon <- sqrt(.Machine$double.eps)
  He <- fit$fit$hessian
  He.eig <- eigen(He, symmetric=TRUE)
  if(min(He.eig$values) < epsilon) He.eig$values[which(He.eig$values < epsilon)] <- 0.0000001
  Vb <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values),He.eig$vectors)  
  
  
  if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0) && fp==FALSE){
  
      HeSh <- He - fit$fit$S.h
      F <- Vb%*%HeSh # diag(SemiParFit$magpp$edf)   # this could be taken from magic as well
      F1 <- fit$magpp$edf1                          # needed for testing
      R <- fit$bs.mgfit$R                           # needed for testing
      Ve <- F%*%Vb                                  # diag(SemiParFit$magpp$Ve) and diag(SemiParFit$magpp$Vb) but need to be careful with dispersion parameters
                                            
  }else{ 
  
  HeSh <- He
  Ve <- Vb
  F <- F1 <- diag(rep(1,dim(Vb)[1]))
  R <- fit$bs.mgfit$R
  
} 

  t.edf <- sum(diag(F))
  
  param <- fit$fit$argument[-c(1:X1.d2)]
  param <- param[1:X2.d2]
  X2s   <- try(predict.gam(gam2, newdata = data, type="lpmatrix"), silent = TRUE)
  if(class(X2s)=="try-error") stop("Check that the range of the covariate values of the selected sample is the same as that of the complete dataset.") 
  
  eta2 <- X2s%*%param
  
  
  dimnames(fit$fit$hessian)[[1]] <- dimnames(fit$fit$hessian)[[2]] <- dimnames(Vb)[[1]] <- dimnames(Vb)[[2]] <- dimnames(HeSh)[[1]] <- dimnames(HeSh)[[2]] <- dimnames(F)[[1]] <- dimnames(F)[[2]] <- dimnames(He)[[1]] <- dimnames(He)[[2]] <- names(fit$fit$argument)   

  
if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0) ){

  edf <- edf1 <- list(0,0,0,0,0)
        
     for(i in 1:5){

       if(i==1) {mm <- l.sp1; if(mm==0) next}
       if(i==2) {mm <- l.sp2; if(mm==0) next} 
       if(i==3) {mm <- l.sp3; if(mm==0) next} 
       if(i==4) {mm <- l.sp4; if(mm==0) next} 
       if(i==5) {mm <- l.sp5; if(mm==0) break}        

          for(k in 1:mm){

              if(i==1){ gam <- gam1; ind <- (gam$smooth[[k]]$first.para):(gam$smooth[[k]]$last.para) } 
              if(i==2){ gam <- gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + X1.d2 } 
              if(i==3){ gam <- gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + X1.d2 + X2.d2 } 
              if(i==4){ gam <- gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + X1.d2 + X2.d2 + X3.d2 } 
              if(i==5){ gam <- gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + X1.d2 + X2.d2 + X3.d2 + X4.d2 } 
              
	      edf[[i]][k]  <- sum(diag(F)[ind])
	      edf1[[i]][k] <- sum(F1[ind])
                        }
                  }
         
  if(l.sp1!=0) names(edf[[1]]) <- names(edf1[[1]]) <- names(gam1$sp)  
  if(l.sp2!=0) names(edf[[2]]) <- names(edf1[[2]]) <- names(gam2$sp)   
  if(l.sp3!=0) names(edf[[3]]) <- names(edf1[[3]]) <- names(gam3$sp)  
  if(l.sp4!=0) names(edf[[4]]) <- names(edf1[[4]]) <- names(gam4$sp) 
  if(l.sp5!=0) names(edf[[5]]) <- names(edf1[[5]]) <- names(gam5$sp) 

}  
  
  
  
  ###############################################
  # Transforming theta back to its original scale 
  ###############################################

  KendTau <- theta <- sigma <- k <- phi <- nu <- NULL
  
  if(length(formula) < 3){

  theta <- theta.tau(BivD = BivD, theta.star = fit$fit$argument["theta.star"])
  KendTau <- theta[,2] # we probably do not want kendall's tau in output
  theta   <- theta[,1]

  names(theta) <- names(KendTau) <- NULL

  # dispersion: variance for normal case 
  # dispersion: 1/k for gamma case

  if(margins[2] %in% c("N", "NB", "PIG", "D", "S", "BB", "NBII", "WARING", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") ) { 
    sigma <- exp(fit$fit$argument["log.sigma"]); names(sigma) <- k <- NULL; phi <- sigma^2 
  }
  
  if(margins[2] %in% c("ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") ) { 
    sigma <- plogis(fit$fit$argument["logi.sigma"]); names(sigma) <- k <- NULL; phi <- sigma^2 
    }

  if(margins[2]=="G") { sigma <- k <- exp(fit$fit$argument["log.shape"]); names(k) <- names(sigma) <- NULL; phi <- k^{-1} }   

  if(margins[2] == "S"){nu <- fit$fit$argument["nu"]; names(nu) <- NULL}
  
  if(margins[2] %in% c("D", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") ){
    nu <- plogis(fit$fit$argument["logi.nu"]); names(nu) <- NULL }

  theta.a <- theta
  sigma.a <- sigma
  phi.a   <- phi
  #k.a     <- k
  nu.a    <- nu

  }
  
  
  if(length(formula) > 2){
  
  sigma <- k <- phi <- nu <- NULL
  
  etatheta <- fit$fit$etatheta 
  etasqv   <- fit$fit$etasqv 
  etak     <- fit$fit$etak
  etak     <- fit$fit$etak  
  
  theta <- theta.tau(BivD = BivD, theta.star = fit$fit$etatheta)
  KendTau <- theta[,2] 
  theta   <- as.matrix(theta[,1])
  
  if(margins[2] %in% c("N", "NB", "PIG", "D", "S", "BB", "NBII", "WARING", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") ) { 
    sigma <- exp(fit$fit$etasqv); phi <- sigma^2
  } 
  
  if(margins[2] %in% c("ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") ) { 
    sigma <- plogis(fit$fit$etasqv); phi <- sigma^2
  } 

  if(margins[2]=="G") { sigma <- k <- exp(fit$fit$etak); phi <- k^{-1} }   

  if(margins[2] == "S") nu <- fit$fit$etanu
  
  if(margins[2] %in% c("D", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") ) nu <- plogis(fit$fit$etanu)

  
  # need several averages here ... for summary function

  theta.a <- mean(as.numeric(theta))
  sigma.a <- mean(as.numeric(sigma))
  phi.a   <- mean(as.numeric(phi))
# k.a     <- mean(as.numeric(k))
  nu.a    <- mean(as.numeric(nu))  

  }  
  
 
  ###########################################################################
  # Final object
  ###########################################################################

       L <- list(fit = fit$fit, gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5,  gam2.1 = gam2.1, 
                 coefficients = fit$fit$argument, weights = weights, sp = fit$sp, 
                 iter.if = fit$iter.if, iter.sp = fit$iter.sp, iter.inner = fit$iter.inner, 
                 start.v = start.v, formula = formula,
                 phi = phi, sigma = sigma, nu = nu, theta = theta, tau = KendTau, 
                 n = n, n.sel = length(gam2$y), bd = bd,
                 X1 = X1, X2 = X2, X1.d2 = X1.d2, X2.d2 = X2.d2, 
                 X3 = X3, X3.d2 = X3.d2, X4 = X4, X4.d2 = X4.d2, X5 = X5, X5.d2 = X5.d2, 
                 l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, l.sp5 = l.sp5, 
                 He = He, HeSh = HeSh, Vb = Vb, F = F, 
                 BivD = BivD, margins = margins,
                 t.edf = t.edf, bs.mgfit = fit$bs.mgfit, 
                 conv.sp = fit$conv.sp, wor.c = fit$wor.c,  
                 eta1 = fit$fit$eta1, eta2 = eta2, # more ETAs here? 
                 y1 = y1, y2 = y2, logLik = -fit$fit$l, fp = fp,
                 gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5,               
                 X2s = X2s, magpp = fit$magpp,
                 sigma.a = sigma.a, phi.a = phi.a, theta.a = theta.a, nu.a = nu.a, gamlss.fit = gamlss.fit,
                 edf = edf, edf11=edf1, R = R, Ve = Ve)

  class(L) <- "SemiParSampleSel"

  L

}
