copulaSampleSel <- function(formula, data = list(), weights = NULL, subset = NULL,
                             BivD = "N", margins = c("probit","N"), gamlssfit = FALSE,
                             fp = FALSE, infl.fac = 1, 
                             rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t"){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  stop("Check next release for tested version of this function.")
  
  i.rho <- sp <- qu.mag <- qu.mag1 <- qu.mag2 <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss <- inde <- spgamlss <- NULL  
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- 0
  ngc <- 2
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- NULL
  sp1 <- sp2 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- NULL  
  sp4 <- gp4 <- gam4 <- X4 <- NULL  
  sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- NULL  
  sp7 <- gp7 <- gam7 <- X7 <- NULL   
  dataset <- NULL
  X2s <- X3s <- X4s <- X5s <- NULL 
  
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","AMH","FGM")
  scc  <- c("C0", "C180", "J0", "J180", "G0", "G180")
  sccn <- c("C90", "C270", "J90", "J270", "G90", "G270")
  m2   <- c("N","GU","rGU","LO","LN","WEI","iG","GA","GAi","BE","FISK")
  m3   <- c("DAGUM","SM")
  
  bl   <- c("probit", "logit", "cloglog")

  if(!is.list(formula)) stop("You must specify a list of equations.")

  l.flist <- length(formula)
    
  if(!(BivD %in% opc)) stop("Error in parameter BivD value. It should be one of: N, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270, F, AMH, FGM.")
  if(!(extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of: t, pC or sED.")
  
  if(!(margins[1] %in% bl) ) stop("Error in first margin value. It can be: probit, logit, cloglog.")
  if(!(margins[2] %in% c(m2,m3)) ) stop("Error in second margin value. It can be: N, GU, rGU, LO, LN, WEI, iG, GA, GAi, DAGUM, SM, BE, FISK.")  
  
  if(l.flist > 2 && margins[2] %in% m2){ if(l.flist!=4) stop("You need to specify four equations.") } 
  if(l.flist > 2 && margins[2] %in% m3){ if(l.flist!=5) stop("You need to specify five equations.") }  

  if( l.flist > 4 && margins[2] %in% m2)     stop("The chosen model can not have more than four equations.")
  if( l.flist > 5 && margins[2] %in% m3)     stop("The chosen model can not have more than five equations.")

 #######################################################################################  
 # formula check  
 #######################################################################################  
  
    if(l.flist > 2){
    
    f3t <- try(formula[[3]][[3]], silent = TRUE)  
    if(class(f3t)!="try-error") stop("The third equation does not require a response.")
    
    	if(l.flist > 3){
    	
    f4t <- try(formula[[4]][[3]], silent = TRUE)  
    if(class(f4t)!="try-error") stop("The fourth equation does not require a response.")  
    
    		if(l.flist > 4){
    		
    f5t <- try(formula[[5]][[3]], silent = TRUE)  
    if(class(f5t)!="try-error") stop("The fifth equation does not require a response.")    
    
    			if(l.flist > 5){
    			
    f6t <- try(formula[[6]][[3]], silent = TRUE)  
    if(class(f6t)!="try-error") stop("The sixth equation does not require a response.")      			
    			
    			                        }    
    					}
    				}    				
                            }

 #######################################################################################  

  mf  <- match.call(expand.dots = FALSE)
  or2 <- as.character(formula[[2]][2])
  ig  <- interpret.gam(formula)

  if( l.flist == 2 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    pred.n <- union(v1,c(v2,or2))
                    }
  
  if( l.flist == 3 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names 
    
    pred.n <- union(v1,c(v2,v3,or2))
                    } 
  
  if( l.flist == 4 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names  
    v4 <- ig[[4]]$pred.names  
    pred.n <- union(v1,c(v2,v3,v4,or2))
                    } 
  
  if( l.flist == 5 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names 
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names 
    pred.n <- union(v1,c(v2,v3,v4,v5,or2))
                    }   
  
  if( l.flist == 6 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names  
    v6 <- ig[[6]]$pred.names
    pred.n <- union(v1,c(v2,v3,v4,v5,v6,or2))
                    }    
  
  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  mf$BivD <- mf$margins <- mf$fp <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
  mf$drop.unused.levels <- TRUE 
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()
        
     indS <- data[, v1[1]]    
     indS[is.na(indS)] <- 0   
     indS <- as.logical(indS)  
     data[indS == FALSE, v2[1]] <- 0.01  
     data <- na.omit(data)

  if(is.null(weights)) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
  
  ct  <- data.frame( c(opc),
                    c(1:14,55,56) 
                     )
  cta <- data.frame( c(opc),
                     c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,55,56) 
                     )                   
  nC  <- ct[which(ct[,1]==BivD),2]
  nCa <- cta[which(cta[,1]==BivD),2]  
  
 ##############################################################  
 # Equation 1
 ##############################################################  
   
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link = margins[1]), gamma=infl.fac, weights=weights, 
                              data=data),list(weights=weights))) 

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 
  if(l.sp1 != 0) sp1 <- gam1$sp
 
 ##############################################################  
 # Equation 2
 ##############################################################   
 # if( margins[2] %in% c("LN") ) formula.eq2 <- update(formula.eq2, (log(.) + mean(log(.)))/2 ~ . )  

 
    inde <- as.logical(y1) 

    formula.eq2r <- formula.eq2   
    y2 <- y2.test <- data[, v2[1]]
    
    if( v2[1] != as.character(formula.eq2r[2]) ) y2.test <- try(data[, as.character(formula.eq2r[2])], silent = TRUE) 
    if(class(y2.test) == "try-error") stop("Please check the syntax for the response of the second equation.")     
    
    if(margins[2] %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK") && min(y2.test[inde], na.rm = TRUE) <= 0) stop("The response of the second margin must be positive.")    
    if(margins[2] %in% c("BE") && (min(y2.test[inde], na.rm = TRUE) <= 0 || max(y2.test[inde], na.rm = TRUE) >= 1) ) stop("The response of the second margin must be in the interval (0,1).")
    
    if( margins[2] %in% c("N","LO","GU","rGU","GAi") )           formula.eq2 <- update(formula.eq2, (. + mean(.))/2 ~ . )  
    if( margins[2] %in% c("LN") )                                formula.eq2 <- update(formula.eq2, (log(.) + mean(log(.)))/2 ~ . )  
    if( margins[2] %in% c("iG","GA","DAGUM","SM","FISK") )       formula.eq2 <- update(formula.eq2, log((. + mean(.))/2) ~ . )    
    if( margins[2] %in% c("WEI") )                               formula.eq2 <- update(formula.eq2, log( exp(log(.) + 0.5772/(1.283/sqrt(var(log(.)))))  ) ~ . )     
    if( margins[2] %in% c("BE") )                                formula.eq2 <- update(formula.eq2, qlogis((. + mean(.))/2) ~ . )    
  
    gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data, subset=inde),list(weights = weights, inde = inde)))
    
    ######
    # TEST
    ######
    X2s <- try(predict.gam(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X2s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    ######
    
    gam2$formula <- formula.eq2r  
    names(gam2$model)[1] <- as.character(formula.eq2r[2])
    X2.d2 <- length(coef(gam2))    
    X2 <- model.matrix(gam2) 
    
    y2 <- y2.test # just in case the above condition was true
    y2 <- y2[inde] 
    if( margins[2] %in% c("LN") ) y2 <- log(y2) 
    
    n.sel <- sum(as.numeric(inde)) 
    l.sp2 <- length(gam2$sp)
    if(l.sp2 != 0) sp2 <- gam2$sp     
    
    cy <- 1 - y1
    
##############################################################  
##############################################################
  
  gp1 <- gam1$nsdf 
  gp2 <- gam2$nsdf
  
##############################################################
# Starting values for dependence parameter
##############################################################

ass.s <- 0.01 

if(BivD %in% scc)  ass.s <-  abs(ass.s)   
if(BivD %in% sccn) ass.s <- -abs(ass.s) 

if(!(BivD %in% c("AMH","FGM"))) i.rho <- BiCopTau2Par(family = nCa, tau = ass.s)
if(  BivD %in% c("AMH","FGM") ) i.rho <- BiCopTau2Par(family = 1, tau = ass.s)

if(BivD %in% c("N","AMH","FGM")) i.rho <- atanh( i.rho )
if(BivD == "F") i.rho <- ifelse( abs(i.rho) < 0.0000001, 0.0000001, i.rho ) 
if(!(BivD %in% c("N","AMH","FGM","F"))) i.rho <- abs(i.rho)

if(BivD %in% c("C0","C180","C90","C270"))                            i.rho <-  log(i.rho)   
if(BivD %in% c("J0","J180","G0","G180","J90","J270","G90","G270"))   i.rho <-  log(i.rho - 1)   

names(i.rho) <- "theta.star"   

                  
##############################################################
# Starting values for whole parameter vector
##############################################################
                      
        par.est <- try( resp.check(y2, margin = margins[2], plots = FALSE, print.par = TRUE), silent = TRUE)
        
        if(class(par.est)=="try-error") {
        
 		if( margins[2] %in% c("N","LN") )             log.sig2 <- log(var(y2))  
		if( margins[2] %in% c("LO") )                 log.sig2 <- log( 3*var(y2)/pi^2 )   
		if( margins[2] %in% c("iG") )                 log.sig2 <- log( var(y2)/mean(y2)^3 )      
		if( margins[2] %in% c("GU","rGU") )           log.sig2 <- log(6*var(y2)/pi^2)   
		if( margins[2] %in% c("WEI") )                log.sig2 <- log( (1.283/sqrt(var(log(y2))))^2 )              
		if( margins[2] %in% c("GA","GAi") )           log.sig2 <- log(var(y2)/mean(y2)^2)           
        	if( margins[2] %in% c("DAGUM","SM","FISK") )  log.sig2 <- log(sqrt(2)) # log(0.01) #         # 0.1     
        	if( margins[2] %in% c("BE"))                  log.sig2 <- qlogis( var(y2)/( mean(y2)*(1-mean(y2)) )  )        
        	
        
        } else log.sig2 <- par.est[2]
        
        names(log.sig2) <- "sigma2.star"
        
        
        
        if( margins[2] %in% c("DAGUM","SM") ){
        
              if( margins[2] %in% c("DAGUM","SM") ){
        	if(class(par.est)=="try-error") log.nu <- log(1) else log.nu <- par.est[3]
                                         }
        	names(log.nu) <- "nu.star"
                                              }
        

        
       		if(margins[2] %in% m2){
       
       		                     start.v <- c(coef(gam1), coef(gam2), log.sig2, i.rho) 
                                    start.v1 <- c(            coef(gam2), log.sig2       ) 
                           
       		                      }                    
        
       		if(margins[2] %in% m3){
       
       		                     start.v <- c(coef(gam1), coef(gam2), log.sig2, log.nu, i.rho) 
       		                    start.v1 <- c(            coef(gam2), log.sig2, log.nu       ) 
                           
       				      }        
         

##############################################################  
    
    if(l.flist == 4){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    nad1 <- "sigma2" 
    nad2 <- "theta" 
    formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2,"~",formula.eq4[2],sep="") )        
    
    set.seed(1)
    sigma2 <- rnorm(n, log.sig2, 0.001) 
    theta  <- rnorm(n, i.rho, 0.001)           
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, subset=inde) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, subset=inde)  
    
    ######
    # TEST
    ######
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X3s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    X4s <- try(predict.gam(gam4, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X4s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    ######    
    
    
    
    l.sp3 <- length(gam3$sp)    
    l.sp4 <- length(gam4$sp)    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, subset=inde); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, subset=inde); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]    

    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    
    
    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
  
    start.v  <- c(coef(gam1), coef(gam2), coef(gam3), coef(gam4) )
    start.v1 <- c(            coef(gam2), coef(gam3)             ) 
    
    
    
  }  
  
  
    if(l.flist == 5){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]     
    nad1 <- "sigma2" 
    nad2 <- "nu" 
    nad3 <- "theta"     
    formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad3,"~",formula.eq5[2],sep="") )  
    
    set.seed(1)
    sigma2 <- rnorm(n, log.sig2, 0.001) 
    nu     <- rnorm(n, log.nu, 0.001)         
    theta  <- rnorm(n, i.rho, 0.001)       
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc, subset=inde) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc, subset=inde)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc, subset=inde)  
    
    ######
    # TEST
    ######
    X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X3s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    X4s <- try(predict.gam(gam4, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X4s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    X5s <- try(predict.gam(gam5, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
    if(class(X5s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?copulaSampleSel for more information.") 
    ######       
    
    
    l.sp3 <- length(gam3$sp)    
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)    
  
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, subset=inde); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, subset=inde); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, subset=inde); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
  
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]     

    if(l.sp3 != 0) sp3 <- gam3$sp
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    
    if(l.sp4 != 0) sp4 <- gam4$sp
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    
    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf      
    
    start.v  <- c(coef(gam1), coef(gam2), coef(gam3), coef(gam4), coef(gam5) )
    start.v1 <- c(            coef(gam2), coef(gam3), coef(gam4)             )     
    
  }   
  
  
  
  
  
##########################################################
  
l.gam1 <- length(coef(gam1))
l.gam2 <- length(coef(gam2))
l.gam3 <- length(coef(gam3))
l.gam4 <- length(coef(gam4))
l.gam5 <- length(coef(gam5))
l.gam6 <- length(coef(gam6))
l.gam7 <- length(coef(gam7))
  
  if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0) && fp==FALSE ){ 
  
                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7)
                 qu.mag <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
                               l.sp1, l.sp2, l.sp3, l.sp4, l.sp5, l.sp6, l.sp7, 
                               l.gam1, l.gam2, l.gam3, l.gam4, l.gam5, l.gam6, l.gam7) 
       
  }     

  
  
  
  if(gamlssfit == TRUE){
  
  
  
  
  if(margins[2] %in% m3 && (l.sp2!=0 || l.sp3!=0 || l.sp4!=0)){  

                                          spgamlss <- c(sp2, sp3, sp4)
                                          qu.mag1 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
                                                         0, l.sp2, l.sp3, l.sp4, 0, 0, 0,     
                                                         0, l.gam2, l.gam3, l.gam4, 0, 0, 0)     
                                                         
                                                              }
                                                              
  if(margins[2] %in% m2 && (l.sp2!=0 || l.sp3!=0)){  

                                          spgamlss <- c(sp2, sp3)
                                          qu.mag1 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
                                                         0, l.sp2, l.sp3, 0, 0, 0, 0,     
                                                         0, l.gam2, l.gam3, 0, 0, 0, 0)     
                                                         
                                                              }                                                              
  
  
  }
  
  

##########################################################


if(missing(parscale)) parscale <- 1   


  respvec <- list(y1 = y1,
                  y2 = y2,
                  y1.y2 = y1.y2, 
                  y1.cy2 = y1.cy2, 
                  cy1.y2 = cy1.y2, 
                  cy1.cy2 = cy1.cy2, 
                  cy1 = cy1,
                  cy = cy, univ = 0)
  
  respvec1 <- respvec; respvec1$univ <- 1  
  
  
  VC <- list(X1 = X1, 
             X2 = X2, 
             X3 = X3,
             X4 = X4, 
             X5 = X5, 
             X6 = X6,
             X7 = X7, 
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             X3.d2 = X3.d2,
             X4.d2 = X4.d2,
             X5.d2 = X5.d2,
             X6.d2 = X6.d2,       
             X7.d2 = X7.d2, 
             gp1 = gp1, 
             gp2 = gp2,
             gp3 = gp3,
             gp4 = gp4, 
             gp5 = gp5,
             gp6 = gp6,  
             gp7 = gp7, nCa = nCa,
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6,    
             l.sp7 = l.sp7, 
             infl.fac = infl.fac,
             weights = weights,
             fp = fp,
             hess = NULL,
             Model = "B", gamlssfit = gamlssfit,
             end = end,
             BivD = BivD,
             nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "NO", ccss = "yes", m2 = m2, m3 = m3, bl = bl, inde = inde,
             X2s = X2s, X3s = X3s, X4s = X4s, X5s = X5s, triv = FALSE) # original n only needed in SemiParBIVProbit.fit
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################

  if(margins[2] %in% m2  ) {func.opt <- bprobgHsContSS  ; func.optUniv <- bprobgHsContUniv}   
  if(margins[2] %in% m3  ) {func.opt <- bprobgHsCont3SS ; func.optUniv <- bprobgHsContUniv3}   
  
  ##########################################################################################################################
  ##########################################################################################################################


  if(gamlssfit == TRUE){ 
  
  
  
  gamlss <- SemiParBIVProbit.fit(func.opt = func.optUniv, start.v = start.v1, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec1, VC = VC, sp = spgamlss, qu.mag = qu.mag1)
                         
    lk <- -gamlss$fit$l
    attr(lk, "nobs") <- n
    attr(lk, "df") <- sum(gamlss$magpp$edf)
    class(lk) <- "logLik"
    gamlss$lk <- lk                         
                         
                         
                         
  # new starting values                       
                         
  if( l.flist == 2 ) start.v <- c(coef(gam1), gamlss$fit$argument, i.rho)
  if( l.flist == 4 ) start.v <- c(coef(gam1), gamlss$fit$argument, coef(gam4))
  if( l.flist == 5 ) start.v <- c(coef(gam1), gamlss$fit$argument, coef(gam5))
  
  if(l.sp2 != 0) sp2 <- gamlss$sp[1:l.sp2] 
  if(l.sp3 != 0) sp3 <- gamlss$sp[l.sp2 + (1:l.sp3)] 
  if( margins[2] %in% m3 ){ if(l.sp4 != 0) sp4 <- gamlss$sp[l.sp2 + l.sp3 + (1:l.sp4)] } 
  
  sp <- c(sp1, sp2, sp3, sp4, sp5, sp6)
  
  gam2$coefficients <- gamlss$fit$argument[1:X2.d2]
  gam2$Vp <- gamlss$magpp$Vb[1:X2.d2,1:X2.d2]
  gam2$sig2 <- 1
  gam2$edf <- gamlss$magpp$edf[1:X2.d2]
                                          
                                          
      if((l.flist - 2) > 1){

          gam3$coefficients <- gamlss$fit$argument[X2.d2 + (1:X3.d2)]
  	  gam3$Vp <- gamlss$magpp$Vb[X2.d2 + (1:X3.d2), X2.d2 + (1:X3.d2)]
  	  gam3$sig2 <- 1
          gam3$edf <- gamlss$magpp$edf[X2.d2 + (1:X3.d2)]  
          
          
          if((l.flist - 2) > 2){ 
          
                       gam4$coefficients <- gamlss$fit$argument[X2.d2 + X3.d2 + (1:X4.d2)]
	       	       gam4$Vp <- gamlss$magpp$Vb[X2.d2 + X3.d2 + (1:X4.d2), X2.d2 + X3.d2 + (1:X4.d2)]
	       	       gam4$sig2 <- 1
                       gam4$edf <- gamlss$magpp$edf[X2.d2 + X3.d2 + (1:X4.d2)]  

                                 }
                                 
                    #       if(l.flist > 4){ 
                    #
                    #         gam5$coefficients <- gamlss$fit$argument[X2.d2 + X3.d2 + X4.d2 + (1:X5.d2)]
	       	    #   gam5$Vp <- gamlss$magpp$Vb[X2.d2 + X3.d2 + X4.d2 + (1:X5.d2), X2.d2 + X3.d2 + X4.d2 + (1:X5.d2)]
	       	    #   gam5$sig2 <- 1
                    #   gam5$edf <- gamlss$magpp$edf[X2.d2 + X3.d2 + X4.d2 + (1:X5.d2)]  
                    #
                    #             }                
          
                                                                            
                           } 
                           
                               
  gamlss$fit$eta2 <- X2s%*%gamlss$fit$argument[1:X2.d2]
  
  if(!is.null(X3)) gamlss$fit$sigma2 <- esp.tr( X3s%*%gamlss$fit$argument[(X2.d2+1):(X2.d2+X3.d2)] , margins[2])$vrb  
  if(!is.null(X4)) gamlss$fit$nu     <- esp.tr( X4s%*%gamlss$fit$argument[(X2.d2 + X3.d2 + 1):(X2.d2 + X3.d2 + X4.d2)], margins[2])$vrb                           
                           
                           
  
  }  
  
  
  ##########################################################################################################################
  ##########################################################################################################################

  SemiParFit <- SemiParBIVProbit.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag) 
                                            
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- copulaSampleSel.fit.post(SemiParFit = SemiParFit, 
                                            VC = VC, 
                                            qu.mag = qu.mag, gam1 = gam1, gam2 = gam2, gam3 = gam3,
                                            gam4 = gam4, gam5 = gam5, gam6 = gam6)
                                            
  SemiParFit <- SemiParFit.p$SemiParFit # useful for SS models, eta2 calculatons etc.
 
  y2.m <- y2  
  if(margins[2] == "LN") y2.m <- exp(y2) # y2.m[inde] <- exp(y2[inde])

  ##########################################################################################################################


if(gc.l == TRUE) gc()

rm(respvec1)


  ##########################################################################################################################


e.v <- min(eigen(SemiParFit$fit$hessian, symmetric=TRUE, only.values = TRUE)$values)
gradi <- round(max(abs(SemiParFit$fit$gradient)),1)

me1 <- "Largest absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?copulaSampleSel."

if(gradi > 0.1 && e.v <= 0){ warning(me1, call. = FALSE); warning(paste(me2,"\n",me3), call. = FALSE)} 
if(gradi > 0.1 && e.v > 0)   warning(paste(me1,"\n",me3), call. = FALSE)
if(gradi < 0.1 && e.v <= 0)  warning(paste(me2,"\n",me3), call. = FALSE)


  ##########################################################################################################################

L <- list(fit = SemiParFit$fit, dataset = dataset, formula = formula,
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7,  
          coefficients = SemiParFit$fit$argument,  iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,l.sp7 = l.sp7, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta = SemiParFit.p$theta, 
          theta.a = SemiParFit.p$theta.a,
          OR = SemiParFit.p$OR, 
          GM = SemiParFit.p$GM,    
          n = n, n.sel = n.sel, 
          X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4 = X4, X5 = X5, X6 = X6, X7 = X7, X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2,           
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf11 = SemiParFit.p$edf11,   # what is this for?
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,
          p11 = SemiParFit$fit$p11, p10 = SemiParFit$fit$p10, p01 = SemiParFit$fit$p01, p00 = SemiParFit$fit$p00, 
          p1 = SemiParFit$fit$p1, p2 = SemiParFit$fit$p2,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, etad = SemiParFit$fit$etad,
          etas = SemiParFit$fit$etas, etan = SemiParFit$fit$etan,
          y1 = y1, y2 = y2.m, 
          BivD = BivD, margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC = nC, hess = TRUE, 
          respvec = respvec, inde = inde, 
          qu.mag = qu.mag, 
          sigma2 = SemiParFit.p$sigma2, sigma2.a = SemiParFit.p$sigma2.a,
          nu = SemiParFit.p$nu,     nu.a = SemiParFit.p$nu.a,
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, 
          X2s = X2s, 
          X3s = X3s,
          X4s = X4s,
          X5s = X5s,
          p1n=SemiParFit.p$p1n , p2n = SemiParFit.p$p2n, 
          VC = VC, Model = NULL, magpp = SemiParFit$magpp,
          gamlss = gamlss, gamlssfit = gamlssfit, Cont = "NO", tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, 
          l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE)

class(L) <- c("copulaSampleSel", "SemiParBIVProbit")

L

}

