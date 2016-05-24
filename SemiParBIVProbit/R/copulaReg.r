copulaReg <- function(formula, data = list(), weights = NULL, subset = NULL,  
                             BivD = "N", margins = c("N","N"), gamlssfit = FALSE,
                             fp = FALSE, infl.fac = 1, 
                             rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t"){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  i.rho <- sp <- qu.mag <- qu.mag1 <- qu.mag2 <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss1 <- gamlss2 <- inde <- spgamlss1 <- spgamlss2 <- NULL  
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- 0
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- NULL
  ngc <- 2
  sp1 <- sp2 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- NULL  
  sp4 <- gp4 <- gam4 <- X4 <- NULL  
  sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- NULL   
  sp7 <- gp7 <- gam7 <- X7 <- NULL    
  
  l.flist <- length(formula)
  
  if(!is.list(formula)) stop("You must specify a list of equations.")

  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","AMH","FGM")
  scc  <- c("C0", "C180", "J0", "J180", "G0", "G180")
  sccn <- c("C90", "C270", "J90", "J270", "G90", "G270")
  m2   <- c("N","GU","rGU","LO","LN","WEI","iG","GA","GAi","BE","FISK")
  m3   <- c("DAGUM","SM")
  bl   <- c("probit", "logit", "cloglog", "cauchit") 

  if(!(BivD %in% opc)) stop("Error in parameter BivD value. It should be one of: N, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270, F, AMH, FGM.")
  if(!(extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of: t, pC or sED.")
  
  if(!(margins[1] %in% c(m2,m3)) ) stop("Error in first margin value. It can be: N, GU, rGU, LO, LN, WEI, iG, GA, GAi, DAGUM, SM, BE, FISK.")  
  if(!(margins[2] %in% c(m2,m3)) ) stop("Error in second margin value. It can be: N, GU, rGU, LO, LN, WEI, iG, GA, GAi, DAGUM, SM, BE, FISK.")  
  
  if(l.flist > 2 && margins[1] %in% m2 && margins[2] %in% m2){ if(l.flist!=5) stop("You need to specify five equations.") } 
  if(l.flist > 2 && margins[1] %in% m2 && margins[2] %in% m3){ if(l.flist!=6) stop("You need to specify six equations.") } 
  if(l.flist > 2 && margins[1] %in% m3 && margins[2] %in% m2){ if(l.flist!=6) stop("You need to specify six equations.") } 
  if(l.flist > 2 && margins[1] %in% m3 && margins[2] %in% m3){ if(l.flist!=7) stop("You need to specify seven equations.") } 
  
  
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
    			
    
        		      if(l.flist > 6){
        			
    f7t <- try(formula[[7]][[3]], silent = TRUE)  
    if(class(f7t)!="try-error") stop("The seventh equation does not require a response.")   
    			
    			                               }
    			
    			
    			                        }    
    					}
    				}    				
                            }

 #######################################################################################  
  
    mf  <- match.call(expand.dots = FALSE)
    or1 <- as.character(formula[[1]][2])
    or2 <- as.character(formula[[2]][2])

ig <- interpret.gam(formula)

  if( l.flist == 2 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    pred.n <- union(v1,c(v2,or1,or2))
                    }
  
  if( l.flist == 3 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names 
    
    pred.n <- union(v1,c(v2,v3,or1,or2))
                    } 
  
  if( l.flist == 4 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names  
    v4 <- ig[[4]]$pred.names  
    pred.n <- union(v1,c(v2,v3,v4,or1,or2))
                    } 
  
  if( l.flist == 5 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names 
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names 
    pred.n <- union(v1,c(v2,v3,v4,v5,or1,or2))
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
    pred.n <- union(v1,c(v2,v3,v4,v5,v6,or1,or2))
                    }  

  if( l.flist == 7 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names  
    v6 <- ig[[6]]$pred.names
    v7 <- ig[[7]]$pred.names
    pred.n <- union(v1,c(v2,v3,v4,v5,v6,v7,or1,or2))
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
 
  n <- dim(data)[1]
        
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
   
    formula.eq1r <- formula.eq1   
    y1 <- y1.test <- data[, v1[1]]
       
    if( v1[1] != as.character(formula.eq1r[2]) ) y1.test <- try(data[, as.character(formula.eq1r[2])], silent = TRUE)
    if(class(y1.test) == "try-error") stop("Please check the syntax for the response of the first equation.") 

          
    if(margins[1] %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK") && min(y1.test, na.rm = TRUE) <= 0) stop("The response of the first margin must be positive.")
    if(margins[1] %in% c("BE") && (min(y1.test, na.rm = TRUE) <= 0 || max(y1.test, na.rm = TRUE) >= 1) ) stop("The response of the first margin must be in the interval (0,1).")
        
    if( margins[1] %in% c("N","LO","GU","rGU","GAi") )           formula.eq1 <- update(formula.eq1, (. + mean(.))/2 ~ . )  
    if( margins[1] %in% c("LN") )                                formula.eq1 <- update(formula.eq1, (log(.) + mean(log(.)))/2 ~ . )  
    if( margins[1] %in% c("iG","GA","DAGUM","SM","FISK") )       formula.eq1 <- update(formula.eq1, log((. + mean(.))/2) ~ . )    
    if( margins[1] %in% c("WEI") )                               formula.eq1 <- update(formula.eq1, log( exp(log(.) + 0.5772/(1.283/sqrt(var(log(.)))))  ) ~ . )     
    if( margins[1] %in% c("BE") )                                formula.eq1 <- update(formula.eq1, qlogis((. + mean(.))/2) ~ . )    
  
  
    gam1         <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
    gam1$formula <- formula.eq1r  
    names(gam1$model)[1] <- as.character(formula.eq1r[2])
     
    y1 <- y1.test # just in case the above condition was true
    if( margins[1] %in% c("LN") ) y1 <- log(y1) 
    
    X1 <- model.matrix(gam1)
    X1.d2 <- dim(X1)[2]
    l.sp1 <- length(gam1$sp)
    if(l.sp1 != 0) sp1 <- gam1$sp
    gp1 <- gam1$nsdf 
    
 ##############################################################
 # Equation 2 
 ##############################################################  

    formula.eq2r <- formula.eq2   
    y2 <- y2.test <- data[, v2[1]]
       
    if( v2[1] != as.character(formula.eq2r[2]) ) y2.test <- try(data[, as.character(formula.eq2r[2])], silent = TRUE)  
    if(class(y2.test) == "try-error") stop("Please check the syntax for the response of the second equation.") 
 
    if(margins[2] %in% c("LN","WEI","WEI2","iG","GA","GAi","DAGUM","SM","FISK") && min(y2.test, na.rm = TRUE) <= 0) stop("The response of the second margin must be positive.")    
    if(margins[2] %in% c("BE") && (min(y2.test, na.rm = TRUE) <= 0 || max(y2.test, na.rm = TRUE) >= 1) ) stop("The response of the second margin must be in the interval (0,1).")
    
    if( margins[2] %in% c("N","LO","GU","rGU","GAi") )            formula.eq2 <- update(formula.eq2, (. + mean(.))/2 ~ . )  
    if( margins[2] %in% c("LN") )                                 formula.eq2 <- update(formula.eq2, (log(.) + mean(log(.)))/2 ~ . )  
    if( margins[2] %in% c("iG","GA","DAGUM","SM","FISK") )        formula.eq2 <- update(formula.eq2, log((. + mean(.))/2) ~ . )    
    if( margins[2] %in% c("WEI") )                                formula.eq2 <- update(formula.eq2, log( exp(log(.) + 0.5772/(1.283/sqrt(var(log(.)))))  ) ~ . )     
    if( margins[2] %in% c("BE") )                                 formula.eq2 <- update(formula.eq2, qlogis((. + mean(.))/2) ~ . )    
  
 
    gam2         <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
    gam2$formula <- formula.eq2r  
    names(gam2$model)[1] <- as.character(formula.eq2r[2])    
    
    y2 <- y2.test 
    if( margins[2] %in% c("LN") ) y2 <- log(y2) 
     
    X2 <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    l.sp2 <- length(gam2$sp)
    if(l.sp2 != 0) sp2 <- gam2$sp 
    gp2 <- gam2$nsdf    
    

##############################################################
# Starting values for dependence parameter
##############################################################

# cor(u1, u2, method = "spearman") where u1 and u2 are on (0,1)
# could possibly improve starting value for theta? Not sure ...


res1 <- residuals(gam1)
res2 <- residuals(gam2)
                  
ass.s <- cor(res1, res2, method = "kendall")

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
         
############         
## first eq.         
############
        
        par.est <- try( resp.check(y1, margin = margins[1], plots = FALSE, print.par = TRUE), silent = TRUE)
        
        if(class(par.est)=="try-error") {
        
 		if( margins[1] %in% c("N","LN") )    log.sig2.1 <- log(var(y1))  
		if( margins[1] %in% c("LO") )        log.sig2.1 <- log( 3*var(y1)/pi^2 )   
		if( margins[1] %in% c("iG") )        log.sig2.1 <- log( var(y1)/mean(y1)^3 )      
		if( margins[1] %in% c("GU","rGU") )  log.sig2.1 <- log(6*var(y1)/pi^2)   
		if( margins[1] %in% c("WEI") )       log.sig2.1 <- log( (1.283/sqrt(var(log(y1))))^2 )              
		if( margins[1] %in% c("GA","GAi") )          log.sig2.1 <- log(var(y1)/mean(y1)^2)           
        	if( margins[1] %in% c("DAGUM","SM","FISK") ) log.sig2.1 <- log(sqrt(2))  # log(0.01) #  log(sqrt(2))       # 0.1  
        	if( margins[1] %in% c("BE"))                 log.sig2.1 <- qlogis( var(y1)/( mean(y1)*(1-mean(y1)) )  )        
        	
        
                                        } else log.sig2.1 <- par.est[2]
        
        names(log.sig2.1) <- "sigma2.1.star"
        
        if( margins[1] %in% m3 ){
        
              if( margins[1] %in% c("DAGUM","SM") ){
        	if(class(par.est)=="try-error") log.nu.1 <- log(1) else log.nu.1 <- par.est[3]
                                         }
              
        	names(log.nu.1) <- "nu.1.star"
        }
        

#############         
## second eq.         
#############            
   
        par.est <- try( resp.check(y2, margin = margins[2], plots = FALSE, print.par = TRUE), silent = TRUE)
        
        if(class(par.est)=="try-error") {
        
 		if( margins[2] %in% c("N","LN") )    log.sig2.2 <- log(var(y2))  
		if( margins[2] %in% c("LO") )        log.sig2.2 <- log( 3*var(y2)/pi^2 )   
		if( margins[2] %in% c("iG") )        log.sig2.2 <- log( var(y2)/mean(y2)^3 )      
		if( margins[2] %in% c("GU","rGU") )  log.sig2.2 <- log(6*var(y2)/pi^2)   
		if( margins[2] %in% c("WEI") )       log.sig2.2 <- log( (1.283/sqrt(var(log(y2))))^2 )              
		if( margins[2] %in% c("GA","GAi") )          log.sig2.2 <- log(var(y2)/mean(y2)^2)           
        	if( margins[2] %in% c("DAGUM","SM","FISK") ) log.sig2.2 <- log(sqrt(2)) # log(0.01) #      # 0.1  
        	if( margins[2] %in% c("BE"))                 log.sig2.1 <- qlogis( var(y2)/( mean(y2)*(1-mean(y2)) )    )        
        	
        
                                       } else log.sig2.2 <- par.est[2]
        
        names(log.sig2.2) <- "sigma2.2.star"
        
        if( margins[2] %in% m3 ){
        
              if( margins[2] %in% c("DAGUM","SM") ){
        	if(class(par.est)=="try-error") log.nu.2 <- log(1) else log.nu.2 <- par.est[3]
                                         }

        	names(log.nu.2) <- "nu.2.star"
        }
        

#############         
## start val.         
#############  

       
       		if(margins[1] %in% m2 && margins[2] %in% m2){
       
       			             start.v <- c(coef(gam1), coef(gam2), log.sig2.1, log.sig2.2, i.rho) 
       		                    start.v1 <- c(coef(gam1),             log.sig2.1                   ) 
                                    start.v2 <- c(            coef(gam2),             log.sig2.2       ) 
                           
       		} 
       		                      
       		                     
       		if(margins[1] %in% m3 && margins[2] %in% m3){
       
       		                       start.v <- c(coef(gam1), coef(gam2), log.sig2.1, log.sig2.2, log.nu.1, log.nu.2, i.rho) 
       		                      start.v1 <- c(coef(gam1),             log.sig2.1,             log.nu.1                 ) 
       		                      start.v2 <- c(            coef(gam2),             log.sig2.2,           log.nu.2       ) 
                           
       					}
       					
       					
       		if(margins[1] %in% m2 && margins[2] %in% m3){
       
       		                     start.v <- c(coef(gam1), coef(gam2), log.sig2.1, log.sig2.2, log.nu.2, i.rho) 
       		                    start.v1 <- c(coef(gam1),             log.sig2.1                             ) 
       		                    start.v2 <- c(            coef(gam2),             log.sig2.2, log.nu.2       ) 
                           
       					}  
       					
       					
       		if(margins[1] %in% m3 && margins[2] %in% m2){
       
       		                     start.v <- c(coef(gam1), coef(gam2), log.sig2.1, log.sig2.2, log.nu.1, i.rho) 
       		                    start.v1 <- c(coef(gam1),             log.sig2.1,             log.nu.1       ) 
       		                    start.v2 <- c(            coef(gam2),             log.sig2.2                 ) 
                           
       					     	
       				}	

##############################################################  
##############################################################  
  
#if(l.flist > 2) seqq <- seq(-0.005, 0.005, length.out = n) 
  
  

    if(l.flist > 2 && margins[1] %in% m2 && margins[2] %in% m2){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]     
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad3   <- "theta"     
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(  nad3,"~",formula.eq5[2],sep="") )   
    
    set.seed(1)
    sigma2.1 <- rnorm(n, log.sig2.1, 0.001) # rep(log.sig2.1,n) + seqq
    sigma2.2 <- rnorm(n, log.sig2.2, 0.001) # rep(log.sig2.2,n) #+ seqq
    theta    <- rnorm(n, i.rho, 0.001)      # rep(i.rho, n)     #+ seqq
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc)     
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)  
    
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
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
    start.v1 <- c(coef(gam1),             coef(gam3)                         )  
    start.v2 <- c(            coef(gam2),             coef(gam4)             )     
    
  }   
  
  
    if(l.flist > 2 && margins[1] %in% m3 && margins[2] %in% m3){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]   
    formula.eq6 <- formula[[6]] 
    formula.eq7 <- formula[[7]]     
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad.1 <- "nu.1" 
    nad.2 <- "nu.2"     
    nad3 <- "theta"     
    
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad.1,"~", formula.eq5[2],sep="") ) 
    formula.eq6 <- as.formula( paste(nad.2,"~", formula.eq6[2],sep="") )     
    formula.eq7 <- as.formula( paste(  nad3,"~",formula.eq7[2],sep="") )      
    
    sigma2.1 <- rep(log.sig2.1,n) #+ seqq
    sigma2.2 <- rep(log.sig2.2,n) #+ seqq
    nu.1     <- rep(log.nu.1,n)   #+ seqq
    nu.2     <- rep(log.nu.2,n)   #+ seqq        
    theta    <- rep(i.rho, n)     #+ seqq
    
    set.seed(1)
    sigma2.1 <- rnorm(n, log.sig2.1, 0.001) # rep(log.sig2.1,n) + seqq
    sigma2.2 <- rnorm(n, log.sig2.2, 0.001) # rep(log.sig2.2,n) #+ seqq
    nu.1     <- rnorm(n, log.nu.1, 0.001)   #rep(log.nu.1,n)   #+ seqq
    nu.2     <- rnorm(n, log.nu.2, 0.001)   #rep(log.nu.2,n)   #+ seqq        
    theta    <- rnorm(n, i.rho, 0.001)      # rep(i.rho, n)     #+ seqq    
    rm(list=".Random.seed", envir=globalenv()) 
    
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc)     
    gam6 <- gam(formula.eq6, data = data, gamma = ngc)    
    gam7 <- gam(formula.eq7, data = data, gamma = ngc)    
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)   
    l.sp6 <- length(gam6$sp)    
    l.sp7 <- length(gam7$sp)  
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
 
    if(l.sp6 != 0){
    ngc <- 2
    while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }   

    if(l.sp7 != 0){
    ngc <- 2
    while( any(round(summary(gam7)$edf, 1) > 1) ) {gam7 <- gam(formula.eq7, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                      
    
    
    
    
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]  
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2]
    X7 <- model.matrix(gam7)
    X7.d2 <- dim(X7)[2]       


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf    
    

    if(l.sp6 != 0) sp6 <- gam6$sp 
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf   
    

    if(l.sp7 != 0) sp7 <- gam7$sp 
    environment(gam7$formula) <- environment(gam2$formula)
    gp7 <- gam7$nsdf       
    
    start.v  <- c(coef(gam1), coef(gam2), coef(gam3), coef(gam4), coef(gam5), coef(gam6), coef(gam7) )
    start.v1 <- c(coef(gam1),             coef(gam3),             coef(gam5)                         )  
    start.v2 <- c(            coef(gam2),             coef(gam4),             coef(gam6)             )     
    
  }    
  
  
    if(l.flist > 2 && margins[1] %in% m2 && margins[2] %in% m3){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]   
    formula.eq6 <- formula[[6]]    
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad.2 <- "nu.2"     
    nad3 <- "theta"     
    
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad.2,"~", formula.eq5[2],sep="") )     
    formula.eq6 <- as.formula( paste(  nad3,"~",formula.eq6[2],sep="") ) 
    
    set.seed(1)
    sigma2.1 <- rnorm(n, log.sig2.1, 0.001) # rep(log.sig2.1,n) + seqq
    sigma2.2 <- rnorm(n, log.sig2.2, 0.001) # rep(log.sig2.2,n) #+ seqq
    nu.2     <- rnorm(n, log.nu.2, 0.001)   #rep(log.nu.2,n)   #+ seqq        
    theta    <- rnorm(n, i.rho, 0.001)      # rep(i.rho, n)     #+ seqq   
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc)     
    gam6 <- gam(formula.eq6, data = data, gamma = ngc)   
    l.sp3 <- length(gam3$sp)    
    l.sp4 <- length(gam4$sp)    
    l.sp5 <- length(gam5$sp)    
    l.sp6 <- length(gam6$sp)  
    
    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
 
    if(l.sp6 != 0){
    ngc <- 2
    while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }   

     
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]  
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2]      


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf    
    

    if(l.sp6 != 0) sp6 <- gam6$sp 
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf   
         
    start.v  <- c(coef(gam1), coef(gam2), coef(gam3), coef(gam4), coef(gam5), coef(gam6) )
    start.v1 <- c(coef(gam1),             coef(gam3)                                     )  
    start.v2 <- c(            coef(gam2),             coef(gam4), coef(gam5)             )     
    
  }    
    
  

    if(l.flist > 2 && margins[1] %in% m3 && margins[2] %in% m2){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
    formula.eq5 <- formula[[5]]   
    formula.eq6 <- formula[[6]]    
    nad2.1 <- "sigma2.1" 
    nad2.2 <- "sigma2.2" 
    nad.1 <- "nu.1"     
    nad3 <- "theta"     
    
    formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
    formula.eq5 <- as.formula( paste(nad.1,"~", formula.eq5[2],sep="") )     
    formula.eq6 <- as.formula( paste(  nad3,"~",formula.eq6[2],sep="") )     
    
    set.seed(1)
    sigma2.1 <- rnorm(n, log.sig2.1, 0.001) # rep(log.sig2.1,n) + seqq
    sigma2.2 <- rnorm(n, log.sig2.2, 0.001) # rep(log.sig2.2,n) #+ seqq
    nu.1     <- rnorm(n, log.nu.1, 0.001)   #rep(log.nu.1,n)   #+ seqq
    theta    <- rnorm(n, i.rho, 0.001)      # rep(i.rho, n)     #+ seqq   
    rm(list=".Random.seed", envir=globalenv()) 
    
    gam3 <- gam(formula.eq3, data = data, gamma = ngc) 
    gam4 <- gam(formula.eq4, data = data, gamma = ngc)   
    gam5 <- gam(formula.eq5, data = data, gamma = ngc)     
    gam6 <- gam(formula.eq6, data = data, gamma = ngc)    
    l.sp3 <- length(gam3$sp)   
    l.sp4 <- length(gam4$sp)   
    l.sp5 <- length(gam5$sp)   
    l.sp6 <- length(gam6$sp)  
    

    if(l.sp3 != 0){
    ngc <- 2
    while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   } 
                   
    if(l.sp4 != 0){
    ngc <- 2
    while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    

    if(l.sp5 != 0){
    ngc <- 2
    while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }                    
 
    if(l.sp6 != 0){
    ngc <- 2
    while( any(round(summary(gam6)$edf, 1) > 1) ) {gam6 <- gam(formula.eq6, data = data, gamma = ngc + 1); ngc <- ngc + 1; if(ngc > 5) break}  
                   }   

     
        
         
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]  
    X5 <- model.matrix(gam5)
    X5.d2 <- dim(X5)[2]  
    X6 <- model.matrix(gam6)
    X6.d2 <- dim(X6)[2]      


    if(l.sp3 != 0) sp3 <- gam3$sp 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    

    if(l.sp4 != 0) sp4 <- gam4$sp 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
    

    if(l.sp5 != 0) sp5 <- gam5$sp 
    environment(gam5$formula) <- environment(gam2$formula)
    gp5 <- gam5$nsdf    
    

    if(l.sp6 != 0) sp6 <- gam6$sp 
    environment(gam6$formula) <- environment(gam2$formula)
    gp6 <- gam6$nsdf   
         
    start.v  <- c(coef(gam1), coef(gam2), coef(gam3), coef(gam4), coef(gam5), coef(gam6) )
    start.v1 <- c(coef(gam1),             coef(gam3),             coef(gam5)             )  
    start.v2 <- c(            coef(gam2),             coef(gam4)                         )     
    
  }    
      
  

##########################################################
# SPs and penalties
##########################################################
  
l.gam1 <- length(coef(gam1))
l.gam2 <- length(coef(gam2))
l.gam3 <- length(coef(gam3))
l.gam4 <- length(coef(gam4))
l.gam5 <- length(coef(gam5))
l.gam6 <- length(coef(gam6))
l.gam7 <- length(coef(gam7))  
  
  
  if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0) && fp==FALSE ){ 
  
                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7)
                 qu.mag <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
                               l.sp1, l.sp2, l.sp3, l.sp4, l.sp5, l.sp6, l.sp7, 
                               l.gam1, l.gam2, l.gam3, l.gam4, l.gam5, l.gam6, l.gam7 ) 
                               
  }
  
  
  
if(gamlssfit == TRUE){


  if( margins[1] %in% m2 && margins[2] %in% m2){   
 
             if(l.sp1!=0 || l.sp3!=0) {
               	spgamlss1 <- c(sp1, sp3)
               	qu.mag1 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
               	               l.sp1, 0, l.sp3, 0, 0, 0, 0, 
               	               l.gam1, 0, l.gam3, 0, 0, 0, 0)  
                                      }   
             
             if(l.sp2!=0 || l.sp4!=0) {
               	spgamlss2 <- c(sp2, sp4)
               	qu.mag2 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
               	               0, l.sp2, 0, l.sp4, 0, 0, 0,
               	               0, l.gam2, 0, l.gam4, 0, 0, 0)    
                                      }     
                                      
                                               }


  if( margins[1] %in% m3 && margins[2] %in% m3){   
 
             if(l.sp1!=0 || l.sp3!=0 || l.sp5!=0) {
               	spgamlss1 <- c(sp1, sp3, sp5)
               	qu.mag1 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
               	               l.sp1, 0, l.sp3, 0, l.sp5, 0, 0, 
               	               l.gam1, 0, l.gam3, 0, l.gam5, 0, 0)  
                                      }   
             
             if(l.sp2!=0 || l.sp4!=0 || l.sp6!=0) {
               	spgamlss2 <- c(sp2, sp4, sp6)
               	qu.mag2 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
               	               0, l.sp2, 0, l.sp4, 0, l.sp6, 0,
               	               0, l.gam2, 0, l.gam4, 0, l.gam6, 0)    
                                      }     
                                      
                                               }
                                               
  if( margins[1] %in% m2 && margins[2] %in% m3){   
 
             if(l.sp1!=0 || l.sp3!=0) {
               	spgamlss1 <- c(sp1, sp3)
               	qu.mag1 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
               	               l.sp1, 0, l.sp3, 0, 0, 0, 0,
               	               l.gam1, 0, l.gam3, 0, 0, 0, 0)  
                                      }   
             
             if(l.sp2!=0 || l.sp4!=0 || l.sp5!=0) {
               	spgamlss2 <- c(sp2, sp4, sp5)
               	qu.mag2 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
               	               0, l.sp2, 0, l.sp4, l.sp5, 0, 0,
               	               0, l.gam2, 0, l.gam4, l.gam5, 0, 0)    
                                      }     
                                      
                                               }                                               

  if( margins[1] %in% m3 && margins[2] %in% m2){   
 
             if(l.sp1!=0 || l.sp3!=0 || l.sp5!=0) {
               	spgamlss1 <- c(sp1, sp3, sp5)
               	qu.mag1 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
               	               l.sp1, 0, l.sp3, 0, l.sp5, 0, 0, 
               	               l.gam1, 0, l.gam3, 0, l.gam5, 0, 0 )  
                                      }   
             
             if(l.sp2!=0 || l.sp4!=0) {
               	spgamlss2 <- c(sp2, sp4)
               	qu.mag2 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
               	               0, l.sp2, 0, l.sp4, 0, 0, 0,
               	               0, l.gam2, 0, l.gam4, 0, 0, 0)    
                                      }     
                                      
                                               }

}
  
##########################################################
##########################################################

if(missing(parscale)) parscale <- 1   


  respvec <- list(y1 = y1,
                  y2 = y2,
                  y1.y2 = NULL, 
                  y1.cy2 = NULL, 
                  cy1.y2 = NULL, 
                  cy1.cy2 = NULL, 
                  cy1 = NULL,
                  cy = NULL, univ = 0)
                  
  respvec2 <- respvec3 <- respvec
  
  respvec2$univ <- 2 
  respvec3$univ <- 3
  
  
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
             gp7 = gp7,
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6, 
             l.sp7 = l.sp7, 
             infl.fac = infl.fac,
             weights = weights,
             fp = fp, gamlssfit = gamlssfit,
             hess = NULL,
             Model = "CC",
             end = end,
             BivD = BivD, nCa = nCa,
             nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "YES", ccss = "no", m2 = m2, m3 = m3, bl = bl, triv = FALSE) # original n only needed in SemiParBIVProbit.fit
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################

  if(margins[1] %in% m2 && margins[2] %in% m2) {func.opt  <- bcont  ; func.opt1 <- func.opt2 <- bprobgHsContUniv }  
  if(margins[1] %in% m3 && margins[2] %in% m3) {func.opt  <- bcont3 ; func.opt1 <- func.opt2 <- bprobgHsContUniv3 } 
  if(margins[1] %in% m2 && margins[2] %in% m3) {func.opt  <- bcont23; func.opt1 <- bprobgHsContUniv;  func.opt2 <- bprobgHsContUniv3 }
  if(margins[1] %in% m3 && margins[2] %in% m2) {func.opt  <- bcont32; func.opt1 <- bprobgHsContUniv3; func.opt2 <- bprobgHsContUniv  } 
  
  ##########################################################################################################################
  ##########################################################################################################################
  # GAMLSS fit
  ##########################################################################################################################
  ##########################################################################################################################



if(gamlssfit == TRUE){ 


  gamlss1 <- SemiParBIVProbit.fit(func.opt = func.opt1, start.v = start.v1, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec2, VC = VC, sp = spgamlss1, qu.mag = qu.mag1) 

  gamlss2 <- SemiParBIVProbit.fit(func.opt = func.opt2, start.v = start.v2, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec3, VC = VC, sp = spgamlss2, qu.mag = qu.mag2)    
 
  #gamlss1$t.edf <- sum(gamlss1$magpp$edf)
  #gamlss2$t.edf <- sum(gamlss2$magpp$edf)
  #gamlss1$n <- gamlss1$n <- n 
  #gamlss1$lk <- -gamlss1$fit$l
  #gamlss2$lk <- -gamlss2$fit$l
  
  
    lk <- -gamlss1$fit$l
    attr(lk, "nobs") <- n
    attr(lk, "df") <- sum(gamlss1$magpp$edf)
    class(lk) <- "logLik"
    gamlss1$lk <- lk
    

    lk <- -gamlss2$fit$l
    attr(lk, "nobs") <- n
    attr(lk, "df") <- sum(gamlss2$magpp$edf)
    class(lk) <- "logLik"
    gamlss2$lk <- lk
      
  
                         
  #########################   
  # updated starting values 
  #########################
  
  
  if(margins[1] %in% m2 && margins[2] %in% m2 && l.flist == 2){
  
  b1 <- gamlss1$fit$argument[1:X1.d2]
  s1 <- gamlss1$fit$argument[X1.d2+1]
  b2 <- gamlss2$fit$argument[1:X2.d2]
  s2 <- gamlss2$fit$argument[X2.d2+1]  
 
  start.v  <- c(b1, b2, s1, s2, i.rho) 
  
  if( l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  }
  
  if(margins[1] %in% m2 && margins[2] %in% m2 && l.flist > 2){
  
  b1 <- gamlss1$fit$argument[1:X1.d2]
  s1 <- gamlss1$fit$argument[(X1.d2+1):(X1.d2+X3.d2)]
  b2 <- gamlss2$fit$argument[1:X2.d2]
  s2 <- gamlss2$fit$argument[(X2.d2+1):(X2.d2+X4.d2)]  
 
  start.v  <- c(b1, b2, s1, s2, coef(gam5)) 
  
  if( l.sp1 != 0 ) sp1 <- gamlss1$sp[1:l.sp1]
  if( l.sp3 != 0 ) sp3 <- gamlss1$sp[(l.sp1 + 1):(l.sp1 + l.sp3)]
  
  if( l.sp2 != 0 ) sp2 <- gamlss2$sp[1:l.sp2]
  if( l.sp4 != 0 ) sp4 <- gamlss2$sp[(l.sp2 + 1):(l.sp2 + l.sp4)]  

  }  
  
  #
  #
  
  if(margins[1] %in% m3 && margins[2] %in% m3 && l.flist == 2){
  
  b1 <- gamlss1$fit$argument[1:X1.d2]
  s1 <- gamlss1$fit$argument[X1.d2+1]
  n1 <- gamlss1$fit$argument[X1.d2+2]
  b2 <- gamlss2$fit$argument[1:X2.d2]
  s2 <- gamlss2$fit$argument[X2.d2+1] 
  n2 <- gamlss2$fit$argument[X2.d2+2]    
 
  start.v  <- c(b1, b2, s1, s2, n1, n2, i.rho) 
  
  if( l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( l.sp2 != 0 ) sp2 <- gamlss2$sp  
  
  
  }  
  
  if(margins[1] %in% m3 && margins[2] %in% m3 && l.flist > 2){
  
  b1 <- gamlss1$fit$argument[1:X1.d2]
  s1 <- gamlss1$fit$argument[(X1.d2+1):(X1.d2+X3.d2)]
  n1 <- gamlss1$fit$argument[(X1.d2+X3.d2+1):(X1.d2+X3.d2+X5.d2)]  
  b2 <- gamlss2$fit$argument[1:X2.d2]
  s2 <- gamlss2$fit$argument[(X2.d2+1):(X2.d2+X4.d2)]
  n2 <- gamlss2$fit$argument[(X2.d2+X4.d2+1):(X2.d2+X4.d2+X6.d2)]  
 
  start.v  <- c(b1, b2, s1, s2, n1, n2, coef(gam7)) 
  
  if( l.sp1 != 0 ) sp1 <- gamlss1$sp[1:l.sp1]
  if( l.sp3 != 0 ) sp3 <- gamlss1$sp[(l.sp1 + 1):(l.sp1 + l.sp3)]
  if( l.sp5 != 0 ) sp5 <- gamlss1$sp[(l.sp1 + l.sp3 + 1):(l.sp1 + l.sp3 + l.sp5)]
  
  if( l.sp2 != 0 ) sp2 <- gamlss2$sp[1:l.sp2]
  if( l.sp4 != 0 ) sp4 <- gamlss2$sp[(l.sp2 + 1):(l.sp2 + l.sp4)]  
  if( l.sp6 != 0 ) sp6 <- gamlss2$sp[(l.sp2 + l.sp4 + 1):(l.sp2 + l.sp4 + l.sp6)]  
  
  }     
  
  #
  #
  
  if(margins[1] %in% m2 && margins[2] %in% m3 && l.flist == 2){
  
  b1 <- gamlss1$fit$argument[1:X1.d2]
  s1 <- gamlss1$fit$argument[X1.d2+1]
  b2 <- gamlss2$fit$argument[1:X2.d2]
  s2 <- gamlss2$fit$argument[X2.d2+1] 
  n2 <- gamlss2$fit$argument[X2.d2+2]    
 
  start.v  <- c(b1, b2, s1, s2, n2, i.rho) 
  
  if( l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( l.sp2 != 0 ) sp2 <- gamlss2$sp
  
  }  
  
  if(margins[1] %in% m2 && margins[2] %in% m3 && l.flist > 2){
  
  b1 <- gamlss1$fit$argument[1:X1.d2]
  s1 <- gamlss1$fit$argument[(X1.d2+1):(X1.d2+X3.d2)] 
  b2 <- gamlss2$fit$argument[1:X2.d2]
  s2 <- gamlss2$fit$argument[(X2.d2+1):(X2.d2+X4.d2)]
  n2 <- gamlss2$fit$argument[(X2.d2+X4.d2+1):(X2.d2+X4.d2+X5.d2)]  
 
  start.v  <- c(b1, b2, s1, s2, n2, coef(gam6)) 
  
  
  if( l.sp1 != 0 ) sp1 <- gamlss1$sp[1:l.sp1]
  if( l.sp3 != 0 ) sp3 <- gamlss1$sp[(l.sp1 + 1):(l.sp1 + l.sp3)]
    
  if( l.sp2 != 0 ) sp2 <- gamlss2$sp[1:l.sp2]
  if( l.sp4 != 0 ) sp4 <- gamlss2$sp[(l.sp2 + 1):(l.sp2 + l.sp4)]  
  if( l.sp5 != 0 ) sp5 <- gamlss2$sp[(l.sp2 + l.sp4 + 1):(l.sp2 + l.sp4 + l.sp5)]  
  
  }     
  
  #
  #  
  
  if(margins[1] %in% m3 && margins[2] %in% m2 && l.flist == 2){
    
    b1 <- gamlss1$fit$argument[1:X1.d2]
    s1 <- gamlss1$fit$argument[X1.d2+1]
    n1 <- gamlss1$fit$argument[X1.d2+2]
    b2 <- gamlss2$fit$argument[1:X2.d2]
    s2 <- gamlss2$fit$argument[X2.d2+1] 
   
    start.v  <- c(b1, b2, s1, s2, n1, i.rho) 
    
  if( l.sp1 != 0 ) sp1 <- gamlss1$sp
  if( l.sp2 != 0 ) sp2 <- gamlss2$sp    
    
    }  
    
  if(margins[1] %in% m3 && margins[2] %in% m2 && l.flist > 2){
    
    b1 <- gamlss1$fit$argument[1:X1.d2]
    s1 <- gamlss1$fit$argument[(X1.d2+1):(X1.d2+X3.d2)]
    n1 <- gamlss1$fit$argument[(X1.d2+X3.d2+1):(X1.d2+X3.d2+X5.d2)]  
    b2 <- gamlss2$fit$argument[1:X2.d2]
    s2 <- gamlss2$fit$argument[(X2.d2+1):(X2.d2+X4.d2)]
   
    start.v  <- c(b1, b2, s1, s2, n1, coef(gam6)) 
   
  if( l.sp1 != 0 ) sp1 <- gamlss1$sp[1:l.sp1]
  if( l.sp3 != 0 ) sp3 <- gamlss1$sp[(l.sp1 + 1):(l.sp1 + l.sp3)]
  if( l.sp5 != 0 ) sp5 <- gamlss1$sp[(l.sp1 + l.sp3 + 1):(l.sp1 + l.sp3 + l.sp5)]
    
  if( l.sp2 != 0 ) sp2 <- gamlss2$sp[1:l.sp2]
  if( l.sp4 != 0 ) sp4 <- gamlss2$sp[(l.sp2 + 1):(l.sp2 + l.sp4)]  
   
    }     
    
  #
  #
  
sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7)  
  
  

}


  ##########################################################################################################################
  ##########################################################################################################################

  SemiParFit <- SemiParBIVProbit.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag) 
                                            
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- copulaReg.fit.post(SemiParFit = SemiParFit, 
                                            VC = VC, qu.mag = qu.mag, gam1 = gam1, gam2 = gam2, gam3 = gam3,
                                            gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7)                                     
 
 
  y1.m <- y1; if(margins[1] == "LN") y1.m <- exp(y1) 
  y2.m <- y2; if(margins[2] == "LN") y2.m <- exp(y2)

  SemiParFit <- SemiParFit.p$SemiParFit  

  ##########################################################################################################################

if(gc.l == TRUE) gc()

rm(respvec2, respvec3)

  ##########################################################################################################################


e.v <- min(eigen(SemiParFit$fit$hessian, symmetric=TRUE, only.values = TRUE)$values)
gradi <- round(max(abs(SemiParFit$fit$gradient)),1)

me1 <- "Largest absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?copulaReg."

#if(gradi > 0.1) warning(me1, call. = FALSE)
#if(e.v <= 0)    warning(me2, call. = FALSE)

if(gradi > 0.1 && e.v <= 0){ warning(me1, call. = FALSE); warning(paste(me2,"\n",me3), call. = FALSE)} 
if(gradi > 0.1 && e.v > 0)   warning(paste(me1,"\n",me3), call. = FALSE)
if(gradi < 0.1 && e.v <= 0)  warning(paste(me2,"\n",me3), call. = FALSE)


  ##########################################################################################################################



L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, gamlss1 = gamlss1, gamlss2 = gamlss2, formula = formula,        
          edf11 = SemiParFit.p$edf11,   ## what is this for? ##
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7,  
          coefficients = SemiParFit$fit$argument, iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta = SemiParFit.p$theta, 
          theta.a = SemiParFit.p$theta.a,  
          sigma21 = SemiParFit.p$sigma21, sigma22 = SemiParFit.p$sigma22, 
          sigma21.a = SemiParFit.p$sigma21.a, sigma22.a = SemiParFit.p$sigma22.a,
          nu1 = SemiParFit.p$nu1, nu2 = SemiParFit.p$nu2, 
          nu1.a = SemiParFit.p$nu1.a, nu2.a = SemiParFit.p$nu2.a,
          X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7,
          X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2,            
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, edf7 = SemiParFit.p$edf7,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, edf1.7 = SemiParFit.p$edf1.7, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, 
          etad=SemiParFit$fit$etad, etas1 = SemiParFit$fit$etas1, etas2 = SemiParFit$fit$etas2,
          y1 = y1.m, y2 = y2.m, 
          BivD = BivD, margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC = nC, 
          respvec = respvec, hess = TRUE,
          qu.mag = qu.mag, 
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, 
          VC = VC, magpp = SemiParFit$magpp,
          gamlssfit = gamlssfit, Cont = "YES",
          tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE)

class(L) <- c("copulaReg","SemiParBIVProbit")

L

}

