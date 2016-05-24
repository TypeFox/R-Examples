guide <- function(){
  cat("Hello! This function will guide you through Bayesian background estimation procedure in 7 simple steps. Be careful, it has limited functionality and can be unstable to missteps and misspells! In this guide some of the parameters are set to their recommended values. Therefore, if you want to adjust them, use corresponding functions 'manually'. If you do not understand the meaning of the requested parameters please refer to the reference manual. \n")

# STEP I
# 1. Read data
# 2. Plot data
  act <- TRUE
  while(act){
    cat("=============\n")  
    cat("1. ")  
    dat <- step1()
    step2(dat)
	cat("\n Do you want to return and re-do this step? Please type 'yes', 'no' or 'exit' to leave the guide	\n")
    answ <- readline("")
    if(answ=="exit"){
      cat("saving results... done!\n")
	  return(dat)	  
	}
    if(answ=="no")    
      act <- FALSE
  }
  cat("\n")
# STEP II  
# 3. trim data
# 4. plot data
  act <- TRUE
  while(act){  
    cat("=============\n")  
    cat("2. ")  
    dat2 <- step3(dat)
    step2(dat2)
	cat("\n Do you want to return and re-do this step? Please type 'yes', 'no' or 'exit' to leave the guide	\n")
    answ <- readline("")
    if(answ=="exit"){
      cat("saving results... done!\n")
	  return(dat)	  
	}
    if(answ=="no"){    
      act <- FALSE
	  dat <- dat2
	}
  }
  cat("\n")


# STEP III  
# 5. Coherent baseline
  act <- TRUE
  while(act){  
    cat("=============\n")  
    cat("3. ")  
    dat <- step5(dat)
	cat("\n Do you want to return and re-do this step? Please type 'yes', 'no' or 'exit' to leave the guide	\n")
    answ <- readline("")
    if(answ=="exit"){
      cat("saving results... done!\n")
	  return(dat)	  
	}
    if(answ=="no")    
      act <- FALSE
  }
  cat("\n")
  

# STEP IV  
# 6. Noise
# 7. Plot
  act <- TRUE
  while(act){  
    cat("=============\n")  
    cat("4. ")  
    dat <- step6(dat)
    step7(dat)
	cat("\n Do you want to return and re-do this step? Please type 'yes', 'no' or 'exit' to leave the guide	\n")
    answ <- readline("")
    if(answ=="exit"){
      cat("saving results... done!\n")
	  return(dat)	  
	}
    if(answ=="no")    
      act <- FALSE
  }
  cat("\n")
  
  
# STEP V  
# 8. Lambda
  act <- TRUE
  while(act){  
    cat("=============\n")  
    cat("5. ")  
    dat <- step8(dat)
	cat("\n Do you want to return and re-do this step? Please type 'yes', 'no' or 'exit' to leave the guide	\n")
    answ <- readline("")
    if(answ=="exit"){
      cat("saving results... done!\n")
	  return(dat)	  
	}
    if(answ=="no")    
      act <- FALSE
  } 
  cat("\n")
  
  
# STEP VI
# 9. Gr
  act <- TRUE
  while(act){  
    cat("=============\n")  
    cat("6. ")  
    dat <- step9(dat)
	cat("\n Do you want to return and re-do this step? Please type 'yes', 'no' or 'exit' to leave the guide	\n")
    answ <- readline("")
    if(answ=="exit"){
      cat("saving results... done!\n")
	  return(dat)	  
	}
    if(answ=="no")    
      act <- FALSE
  } 
  cat("\n")
  
  
# STEP VII
# 10. fit params
  act <- TRUE
  while(act){  
    cat("=============\n")  
    cat("7. ")  
    ctrl <- step10()
	cat("\n Do you want to return and re-do this step? Please type 'yes', 'no' or 'exit' to leave the guide	\n")
    answ <- readline("")
    if(answ=="exit"){
      cat("saving results... done!\n")
	  return(dat)	  
	}
    if(answ=="no")    
      act <- FALSE
  } 
  cat("\n")

  
  cat("That's it! Press enter to start fit. Be patient, it may take a while \n")
  tmp <- scan()
  
  fit.res <- do.fit(data=dat, bounds.lower=ctrl$bl, bounds.upper=ctrl$bu, 
      scale=ctrl$sc, knots.x=ctrl$kx, knots.n=ctrl$kn, p.bkg=.5, 
      stdev=TRUE, control=ctrl$ctrl, save.to=ctrl$st)  
  

# Plot fit results
  cat("\n Do you want to plot background estimation (to do this packages 'ggplot2' and 'gridExtra' have to be installed)? Please type 'yes' or 'no' \n")
  answ <- readline("")
  if(answ=="yes")
    mPlot.results(fit.res)  

# Calculate and plot GR
  gr<-NA
  cat("\n Do you want to calculate and plot corrected PDF (to do this package 'ggplot2' have to be installed)? Please type 'yes' or 'no' \n")
  answ <- readline("")
  if(answ=="yes"){
    if(is.null(dat$Gr$rho.0)){
      cat("Please provide the atomic number density: \n")
      rho.0 <- readline("")
	  rho.0 <- as.numeric(unlist(strsplit(rho.0, ",")))[1]   
    }
    else
  	  rho.0 <- dat$Gr$rho.0
	  
   gr <- calc.Gr(fit.res, rho.0, r.max=10)
  }
# Save results
  cat("\n Do you want to save fit results to a text file? If yes, please type file name. If no, press enter \n")
  answ <- readline("")
  if(answ!="")
    write.fit.results(fit.res, file = answ) 


# Last note...
  cat("Finishing the guide... \n If this guide was called as 'myVar <- guide()' then: \n")
  cat("myVar$data contains experimental data, estimated noise, lambda and baseline. See set.data in reference manual for details \n")
  cat("myVar$fit.res contains results of the fit. See do.fit \n")
  cat("myVar$Gr contains the calculated corrected PDF. See calc.Gr \n")
  if(answ!="")
    cat("file ", answ, "contains fit results in a text format \n")  
  if(ctrl$st!="")
    cat("file ", ctrl$st, "contains fit results in a R format \n")
  
  return(list(fit.res=fit.res, data=dat, gr=gr))
}

#########################################################
#########################################################

step1 <- function(){
  cat("First step is to read data. Do you want to read data from \n 1. text file \n or \n 2. .sqb-file, returned by PDFgetN? \n Type 1 or 2 and press Enter. \n")
  step1 <- readline("")
  step1 <- as.numeric(unlist(strsplit(step1, ",")))   
  if(step1==1){
    cat("Provide file name (): \n")
    step1.name <- readline("")
    dat <- read.data(file=step1.name)
  }
  else if(step1==2){
    cat("Provide file name: \n")
    step1.name <- readline("")
    dat <- read.sqb(file=step1.name)
  }
  else{
    stop("Please type 1 or 2.")
  }
  return(dat)
}

#########################################################
step2 <- function(dat){
  cat("Do you want to plot data? Please type 'yes' or 'no' (no quotes) \n")
  step2 <- readline("")
  if(step2=="yes"){
    cat("Provide plotting region on this step. \n Enter x coordinate minimum \n") 
    step2.x1 <- readline("")
    step2.x1 <- as.numeric(unlist(strsplit(step2.x1, ",")))   
	
    cat("Enter x coordinate maximum \n") 
    step2.x2 <- readline("")
    step2.x2 <- as.numeric(unlist(strsplit(step2.x2, ",")))  
	
    cat("Enter y coordinate minimum \n") 
    step2.y1 <- readline("")
    step2.y1 <- as.numeric(unlist(strsplit(step2.y1, ",")))   
	
    cat("Enter y coordinate maximum \n") 
    step2.y2 <- readline("")
    step2.y2 <- as.numeric(unlist(strsplit(step2.y2, ",")))  
		 
	plot(dat$x, dat$y, t="l", xlab="x", ylab="y", 
	     xlim=c(step2.x1, step2.x2), ylim=c(step2.y1, step2.y2))
  }
}



#########################################################
step3 <- function(dat){
  cat("Do you want to truncate data (x-region)? Please type 'yes' or 'no'\n")
  step3 <- readline("")
  if(step3=="yes"){  
    cat("Enter new x minimum value \n")
    step3.min <- readline("")
    step3.min <- as.numeric(unlist(strsplit(step3.min, ",")))   

    cat("Enter new x maximum value \n")
    step3.max <- readline("")
    step3.max <- as.numeric(unlist(strsplit(step3.max, ",")))   
    dat <- trim.data(dat, x.min=step3.min, x.max=step3.max)
  }
  return(dat)
}
	
#########################################################
step5 <- function(dat){
  cat("Do your data contain smooth baseline that shouldn't be subtracted (for example elastic coherent baseline in neutron scattering) and wasn't specified on step 1? Please type 'yes' or 'no'")
  step5 <- readline("")
  if(step5=="yes"){  
    cat("For neutron total scattering experiment we can calculate smooth baseline as a sum of elastic coherent scattering and Laue diffuse scattering. If you have other form of baseline please provide it on step 1 (include it in your text file as a third column with header 'SB') or specify text file that contains it. For neutron scattering experiment, do you \n 1. know APD values \n 2. want to fit ADP values \n 3. it isn't a neutron scattering experiment (baseline will be set to 0)? \n Please type 1, 2 or 3 \n")
    step5.choice <- readline("")
    step5.choice <- as.numeric(unlist(strsplit(step5.choice, ",")))  
 
	if(step5.choice==3){
      if(any(is.na(dat$SB)))
        dat$SB <- rep(0, length(dat$x))
	}
    if((step5.choice==1) || (step5.choice==2)){
	  cat("Type number of atoms of each type per unit cell divided by space and press double enter (e.g. for  NaCl you should type '4 4') \n")
	  n.atoms <- scan()
	  cat("Type neutron scattering length for atoms of each type divided by space and press double enter (e.g. for  NaCl you should type '3.63 9.58') \n")	  
	  scatter.length <- scan()
	  cat("Do you want to use single value for the ADP for all atom types? Please type 'yes' or 'no' \n")
      step5.oneADP <- readline("")
      if(step5.oneADP=="yes")
	    oneADP=TRUE
	  else
	    oneADP=FALSE
	  if(step5.choice==1){
        cat("Please provide the ADP(s) and press double enter \n") 
        step5.ADP <- scan()  
		dat <- set.SB(dat, fit=FALSE, oneADP=oneADP, n.atoms=n.atoms, 
          scatter.length=scatter.length, ADP=step5.ADP)
      }
      if(step5.choice==2){
        cat("Please provide the limits for the ADP fit (upper and lower bounds, divided by space) \n") 
        step5.ADP.lim <- scan()  
		dat <- set.SB(dat, fit=TRUE, oneADP=oneADP, n.atoms=n.atoms, 
          scatter.length=scatter.length, ADP.lim=step5.ADP.lim)	  
	  }
	}
  }
  return(dat)
}




#########################################################
step6 <- function(dat){
  cat("Although noise in diffraction experiments is per se Poisson, various corrections can destroy its structure. We suggest considering the experimental uncertainty as having  Gaussian distribution with x-dependent amplitude. Splitting x-region into N segments and estimating Gaussian standard deviation over these segments allows us to approximate the true noise-distribution. The other way to approximate noise level is to consider it uniform. In that case the best approximation can be obtained on a signal-free region, i.e. on a region that contains only background. \n Please type integer number N if you want to divide x-range into N segments for independent noise-level estimation or type bounds for a signal-free region (two numbers divided by space), and press double enter \n")
  step6 <- scan("")
  cat("\n Thanks!\n")
  if(length(step6)==1)
    dat <- set.sigma(dat, n.regions=step6)
  
  if(length(step6)==2)
    dat <- set.sigma(dat, x.bkg.only=step6)
  return(dat)
}




#########################################################
step7 <- function(dat){
  cat("\n Do you want to plot data +/- 2 sd for estimated noise level? Please type 'yes' or 'no'\n")
  step7 <- readline("")
  if(step7=="yes"){
    cat("Provide plotting region on this step. \n Enter x coordinate minimum \n") 
    step7.x1 <- readline("")
    step7.x1 <- as.numeric(unlist(strsplit(step7.x1, ",")))   
	
    cat("Enter x coordinate maximum \n") 
    step7.x2 <- readline("")
    step7.x2 <- as.numeric(unlist(strsplit(step7.x2, ",")))  
	
    cat("Enter y coordinate minimum \n") 
    step7.y1 <- readline("")
    step7.y1 <- as.numeric(unlist(strsplit(step7.y1, ",")))   
	
    cat("Enter y coordinate maximum \n") 
    step7.y2 <- readline("")
    step7.y2 <- as.numeric(unlist(strsplit(step7.y2, ",")))  
		 
	plot(dat$x, dat$y, t="l", xlab="x", ylab="y", 
	     xlim=c(step7.x1, step7.x2), ylim=c(step7.y1, step7.y2))
		 
    lines(dat$x, dat$smoothed, col=2)
    lines(dat$x, dat$smoothed + 2* dat$sigma, col=4)
    lines(dat$x, dat$smoothed - 2* dat$sigma, col=4)
  }
}


#########################################################
step8 <- function(dat){
  cat("On this step we estimate the mean signal magnitude (lambda). lambda is calculated as a linear piecewise function which is equal to lambda_0 outside the [x.min, x.max] region. Inside this region lambda is approximated by a line connecting points (x_1;lambda_1) and (x_2;lambda_2). Estimate these parameters to obtain a line that connects centres of the two most distant peaks. High accuracy isn't required on this step. \n Type lambda_0, lambda_1, lambda_2, x_1, and x_2 divided by spaces and press double enter")
  
  step8 <- scan("")
  dat <- set.lambda(dat, lambda_0=step8[1], lambda_1=step8[2], 
    lambda_2=step8[3], x_1=step8[4], x_2=step8[5])
  return(dat)
}


#########################################################
step9 <- function(dat){
  cat("Do you want to include information on low-r G(r) behaviour in Bayesian model (see reference manual for details)? Please type 'yes' or 'no'\n")
  step9 <- readline("")
  if(step9=="yes"){
    cat("Provide the atomic number density: \n")
    step9.rho <- readline("")
	step9.rho <- as.numeric(unlist(strsplit(step9.rho, ",")))   
 
    cat("Indicate bounds for a low-r peak-free region (normally 0..1-2 Angstrom). Type two numbers divided by space and press double enter: \n")
    step9.r <- scan()
   	
    dat <- set.Gr(dat, r1=seq(step9.r[1], step9.r[2], 0.005), rho.0=step9.rho, 
	  type1="gaussianNoise")  
  
  }
  return(dat)
}

#########################################################
step10 <- function(){
  cat("Now the data are ready. Lets specify fit parameters. \n")
  
  cat("Enter bounds for background estimation (you want lower and upper bounds to be some smaller and bigger than real minimum and maximum background values, respectively). Put two numbers divided by space and press double enter \n")
  step10.bd <- scan("")
  
  cat("Enter bounds for normalization parameter. If you don't want normalization parameter to be fitted type '1 1' (no quotes) \n")
  step10.sc <- scan("")
  
  cat("Enter spline knot positions. Put more knots in the region where background demonstrates less smooth behaviour. To use N equidistant knots simply type integer number N \n")
  step10.kn <- scan("")
  if(length(step10.kn)==1){
    knots.x <- NA
	knots.n <- step10.kn
  }
  else{
    knots.x <- step10.kn
	knots.n <- length(knots.x)
  }
  cat("Enter the maximum iteration (population generation) allowed and number of population members (NP). For the most tasks it is best to set NP to be at least 10-15 times the length of the parameter vector, which includes spline knot positions, and, optionally, normalization and ADP parameters. Type two numbers divided by space and press double enter \n")
  step10.pp <- scan("")
  step10.cl <- set.control(NP=step10.pp[2], itermax=step10.pp[1], parallelType=1)
  
  cat("Enter name of the file where the results will be saved. The usual extension is '.RData'. If you don't want to save results to file simply press enter \n") 
  step10.fi <- readline("")
	
  ctrl <- list(bl=step10.bd[1], bu=step10.bd[2], scale=step10.sc, 
              kx=knots.x, kn=knots.n, ctrl=step10.cl, st=step10.fi)
  return(ctrl)
}
 

