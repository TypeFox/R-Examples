list2ascii <- function(x,file = paste(deparse(substitute(x)), ".txt", sep = "")) {
	## see in http://tolstoy.newcastle.edu.au/R/help/05/12/17151.html
  ## by Michael Prager
  # MHP July 7, 2004
  # R or S function to write an R list to an ASCII file.
  # This can be used to create files for those who want to use
  # a spreadsheet or other program on the data.
  #
  tmp.wid <- getOption("width") # save current width
  options(width = 10000)        # increase output width
  sink(file)                    # redirect output to file
  print(x)                      # print the object
  sink()                        # cancel redirection
  options(width = tmp.wid)      # restore linewidth
  return(invisible(NULL))       # return (nothing) from function
}


################################################################################################################################
############################### internal functions for chi2Corr() and chi2CorrGUI() ############################################
################################################################################################################################

obsdata_chi2corr <- function(formula, data, name1, name2) {
  ## observed contingency table
  tabobs <- table(data[[name1]], data[[name2]])
  
  ## calculates the theoretical frequencies according to the risk factors
  formula1 <- paste(name1, "~", formula, sep = "")
  formula2 <- paste(name2, "~", formula, sep = "")
  glm1 <- glm(formula1, binomial, data, na.action = na.omit) 
  glm2 <- glm(formula2, binomial, data, na.action = na.omit)
  pvir1 <- glm1$fitted.values 
  pvir2 <- glm2$fitted.values
  tabth <- cbind(c((sum((1 - pvir1) * (1 - pvir2))), (sum(pvir1 * (1 - pvir2)))), c((sum((1 - pvir1) * pvir2)), (sum(pvir1 * pvir2))))
  
  ## Observed corrected chi-square:
  chi2corrobs <- round(sum((tabobs - tabth) ^ 2 / tabth), digits = 5)
  tabth <- round(tabth, digits = 2)
  colnames(tabth) <- colnames(tabobs)
  rownames(tabth) <- rownames(tabobs)
  if(!is.numeric(chi2corrobs))
    stop("fatal error in calculation of chi-square between the real data and the theoretical infectious status of each real individual")
  
  return(list(chi2corrobs = chi2corrobs, tabobs = tabobs, tabth = tabth, pvir1 = pvir1, pvir2 = pvir2))
}


## Calculates the bootstrapped corrected chi-square using a risk factor matrix (data) and 2 vectors of serological status (sero1, sero2):
chi2corrboot <- function(data, formula, sero1, sero2) {
	df <- as.data.frame(cbind(data, sero1 = sero1, sero2 = sero2))
	f1 <- paste("sero1~", formula, sep = "")
	f2 <- paste("sero2~", formula, sep = "")
	glm1 <- glm(f1, binomial, df) 
	glm2 <- glm(f2, binomial, df)
	p.glm1 <- glm1$fitted.values 
	p.glm2 <- glm2$fitted.values
  
	## theoretical and observed contingency tables:
	tabth <- cbind(c((sum((1 - p.glm1) * (1 - p.glm2))), (sum(p.glm1 * (1 - p.glm2)))), c((sum((1 - p.glm1) * p.glm2)), (sum(p.glm1 * p.glm2)))) 
	tabobs <- table(sero1 = sero1, sero2 = sero2) 
  
	## corrected chi-square:
	return(sum((tabobs - tabth) ^ 2 / tabth)) 
}


simudata_chi2corr <- function(formula, data, name1, name2, nbsimu, pvir1, pvir2, chi2corrobs) {
  ## Dataframe with only the risk factors (without the serological status columns):
  data.fac <- subset(data, select = c(- which(colnames(data) == name1), - which(colnames(data) == name2)))
  
  chi2corrsim <- numeric(length = nbsimu) ## Vector of the bootstraped chi-square 
  ## Creates the serological status vectors for parasites 1 and 2:
  for(simu in 1:nbsimu) {
    nn <- NROW(data)
    sero1 <- numeric(nn)
    sero2 <- numeric(nn)	
    sero1[runif(nn) <= pvir1] <- 1
    sero2[runif(nn) <= pvir2] <- 1	
  	chi2corrsim[simu] <- chi2corrboot(data = data.fac, formula = formula, sero1 = sero1, sero2 = sero2) ## Vector of bootstrapped corrected chi-square
  }
  
  dispcoeff <- round(sum(chi2corrsim) / nbsimu, digits = 5) ## Dispersion coefficient
  pval1 <- round(1 - pchisq(chi2corrobs / dispcoeff, 1), digits = 5) ## P-value considering that the corrected chi2 is proportional to a chi2 with 1 df
  pval2 <- round((1 + sum(chi2corrsim > chi2corrobs)) / (1 + nbsimu), digits = 5) ## P-value estimated from the bootstrapped values of the corrected chi2
  chi2corrsim <- unlist(chi2corrsim)
  
  return(list(chi2corrsim = chi2corrsim, pval1 = pval1, pval2 = pval2, dispcoeff = dispcoeff))
}



################################################################################################################################
########################### internal functions for chi2CorrAge() and chi2CorrAgeGUI() ##########################################
################################################################################################################################

## return the vector of transmission rates for one parasite and the one of probabilities to be sensible
SensTransMatrix <- function(para, listmodel, rate, agenum, a) {
  n <- nrow(listmodel[[1]])
  nparam <- ncol(listmodel[[1]])
  agemax <- length(listmodel)
  
  lambda <- matrix(NA, nrow = n, ncol = agemax)
  propS <- matrix(NA, nrow = n, ncol = agemax)
  propS[, 1] <- rep(1, n)
  
  for(ageindex in 1:agemax) {
    lambda[, ageindex] <- exp(listmodel[[ageindex]] %*% as.numeric(para))
    if(ageindex < agemax) {
      propS[, ageindex + 1] <- (propS[, ageindex] - rate / (lambda[, ageindex] + rate)) * 
        exp(- (lambda[, ageindex] + rate) * (a[ageindex + 1] - a[ageindex])) + rate / (lambda[, ageindex] + rate)
    }
  }
  
  lambda_age <- sapply(1:n, function(x) {lambda[x, col(lambda)[x, agenum[x]]]})
  propS_age <- sapply(1:n, function(x) {propS[x, col(propS)[x, agenum[x]]]})
  return(list(lambda_age, propS_age))
}


## parameters estimation 
EstimParam <- function(paranum, rate, listmodel, agenum, v0, tol = 0.00000001, maxit = 50000, a, mort) {
  estimparam <- v0
  
  ## the function to be minimized
  minuslogL <- function(v) {
    senstransmatrix <- SensTransMatrix(v, listmodel, rate, agenum, a)
    lambda_age <- senstransmatrix[[1]] 
    propS_age <- senstransmatrix[[2]]
    proba <- (propS_age - rate / (lambda_age + rate)) * exp((lambda_age + rate) * a[agenum]) * 
      (exp(-(lambda_age + rate + mort[agenum]) * a[agenum]) - exp(-(lambda_age + rate + mort[agenum]) * a[agenum + 1])) / 
      (exp(- mort[agenum] * a[agenum]) - exp(- mort[agenum] * a[agenum + 1])) * 
      mort[agenum] / (mort[agenum] + lambda_age + rate) + rate / (rate + lambda_age)
    out <- sum(- paranum * log(1 - proba) - (1 - paranum) * log(proba))
    return(out)
  }
  
  ## do the optimisation until the convergence
  conv <- 1
  while((conv != 0)) { ## conv == 0 indicates successful completion
    o <- optim(estimparam, minuslogL, method = "Nelder-Mead", control = list(maxit = maxit, reltol = tol)) ## Nelder-Mead method is based on a simplex algorithm
    conv <- o$convergence
    estimparam <- o$par
  }
  
  estimparam <- as.data.frame(t(estimparam))
  colnames(estimparam) <- colnames(listmodel[[1]])
  return(estimparam)
}


## model construction with each age class 
ModelClass <- function(para, formula, data, agemax, nameage) {
  list_model_age <- list()
  data_age <- data
  for(ageindex in 1:agemax) {
    data_age[[nameage]] <- factor(ageindex, levels = levels(data[[nameage]])) ## replace age by ageindex  
    list_model_age[[ageindex]] <- model.matrix(object = as.formula(paste(para, "~", formula)), data = data_age)
  }
  return(list_model_age)
}


calcInfectProba <- function(data, formula, namepara1, namepara2, nameage, w1, w2, mort, a, v0para1, v0para2) {
  
  ## change age variable to become factor
  agenum <- as.numeric(data[[nameage]])
  agemax <- max(agenum)
  data[[nameage]] <- as.factor(data[[nameage]])
  levels(data[[nameage]]) <- as.character(1:agemax)
  
  para1num <- as.numeric(as.vector(data[[namepara1]]))
  para2num <- as.numeric(as.vector(data[[namepara2]]))
  testpara1num <- sum(names(table(para1num)) %in% c("0", "1")) == length(table(para1num))
  testpara2num <- sum(names(table(para2num)) %in% c("0", "1")) == length(table(para2num))
  if((testpara1num == FALSE) | (testpara1num == FALSE))
    stop("para1num and para2num must contain only '0' and '1' values")
  
  ## complete model construction with each age class and for each parasite
  lmodelpara1 <- ModelClass(para = namepara1, formula = formula, data = data, agemax = agemax, nameage = nameage)
  lmodelpara2 <- ModelClass(para = namepara2, formula = formula, data = data, agemax = agemax, nameage = nameage)
  
  ## partial model construction for unknown start point 
  if((sum(is.na(v0para1)) != 0) & (sum(is.na(v0para2)) != 0)) {
    if(length(levels(data[[nameage]])) == 1) {  ## only one age class, unknown age
      v0para1 <- rep(0, ncol(lmodelpara1[[1]]))
      v0para2 <- rep(0, ncol(lmodelpara1[[1]]))
    } else {
      partiallmodelpara1 <- ModelClass(para = namepara1, formula = nameage, data = data, agemax = agemax, nameage = nameage)
      partiallmodelpara2 <- ModelClass(para = namepara2, formula = nameage, data = data, agemax = agemax, nameage = nameage)
      
      ## parameters estimation for each parasite
      v0para1 <- rep(0, ncol(partiallmodelpara1[[1]]))
      v0para2 <- v0para1
      
      estimparam1 <- EstimParam(paranum = para1num, rate = w1, listmodel = partiallmodelpara1, agenum = agenum, v0 = v0para1, a = a, mort = mort)
      estimparam2 <- EstimParam(paranum = para2num, rate = w2, listmodel = partiallmodelpara2, agenum = agenum, v0 = v0para2, a = a, mort = mort)
      id <- sapply(colnames(estimparam1), function(X) {which(colnames(lmodelpara1[[1]]) == X)})
      
      v0para1 <- rep(0, ncol(lmodelpara1[[1]]))
      v0para2 <- rep(0, ncol(lmodelpara1[[1]]))
      v0para1[id] <- as.numeric(estimparam1)
      v0para2[id] <- as.numeric(estimparam2)
    }
  }
  
  ## parameters estimation for each parasite
  estimparam1 <- EstimParam(paranum = para1num, rate = w1, listmodel = lmodelpara1, agenum = agenum, v0 = v0para1, a = a, mort = mort)
  estimparam2 <- EstimParam(paranum = para2num, rate = w2, listmodel = lmodelpara2, agenum = agenum, v0 = v0para2, a = a, mort = mort)
  
  ## coinfection. calcul de la theoritical contengency table
  senstransmatrix_age1 <- SensTransMatrix(para = estimparam1, listmodel = lmodelpara1, rate = w1, agenum = agenum, a) ## list of two vectors with length equal to number of rows of data
  A0_age1 <- w1 / (senstransmatrix_age1[[1]] + w1) ## vector with length equal to number of rows of data
  Q0_age1 <- senstransmatrix_age1[[2]] - A0_age1 ## vector with length equal to number of rows of data
  
  senstransmatrix_age2 <- SensTransMatrix(para = estimparam2, listmodel = lmodelpara2, rate = w2, agenum = agenum, a) ## list of two vectors with length equal to number of rows of data
  A0_age2 <- w2 / (senstransmatrix_age2[[1]] + w2) ## vector with length equal to number of rows of data
  Q0_age2 <- senstransmatrix_age2[[2]] - A0_age2 ## vector with length equal to number of rows of data
  
  ## matproba include for each individual and about two virus the probability to be S(Susceptible)-S, S-I(Infected), I-S and I-I 
  ## the sum by row is equal to 1
  matprobainfect <- matrix(NA, nrow(data), 4)
  
  for(para1num in 1:2) {
    for(para2num in 1:2) {
      
      if(para1num == 1) {
        A_age1 <- A0_age1 
        Q_age1 <- Q0_age1
      } else {
        A_age1 <- 1 - A0_age1 
        Q_age1 <- - Q0_age1
      }
      
      if(para2num == 1){
        A_age2 <- A0_age2 
        Q_age2 <- Q0_age2
      } else {
        A_age2 <- 1 - A0_age2 
        Q_age2 <- - Q0_age2
      }
      
      ## c12, c1, c2 and c are coefficients for coinfection calculation
      c12 <- senstransmatrix_age1[[1]] + senstransmatrix_age2[[1]] + w1 + w2 + mort[agenum]
      c1 <- senstransmatrix_age1[[1]] + w1 + mort[agenum]
      c2 <- senstransmatrix_age2[[1]] + w2 + mort[agenum]
      
      delta_age <- a[2:length(a)] - a[1:(length(a) - 1)]
      c <- ((Q_age1 * Q_age2 * mort[agenum] / c12 * (1 - exp(- c12 * delta_age[agenum]))) +
	   (Q_age1 * A_age2 * mort[agenum] / c1 * (1 - exp(- c1 * delta_age[agenum]))) + 
           (A_age1 * Q_age2 * mort[agenum] / c2 * (1 - exp(- c2 * delta_age[agenum])))) / (1 - exp(- mort[agenum] * delta_age[agenum])) +
 	   A_age1 * A_age2
        
      matprobainfect[, (para1num - 1) * 2 + para2num] <- c
    }
  }
  
  if(!all.equal(rowSums(matprobainfect), rep(1, nrow(matprobainfect))))
    stop("fatal error in initial calculation of infectious probabilities")
  
  return(list(estimparam1 = estimparam1, estimparam2 = estimparam2, matprobainfect = matprobainfect))
}


obsdata_chi2corrage <- function(formula, data, name1, name2, nameage, w1, w2, mort, a, v0para1, v0para2) {

  ## v0para1 and v0para2 are the starting vector for optimisation step. on the begening, they are unknown
  
  ## observed contingency table
  tabobs <- table(data[[name1]], data[[name2]])
  
  ## infectprobaR yieds the matrix of theoretical infection probabilities for each individual ($matprobainfect)
  ## and the estimated parameters for the model with each parasite studied ($estimparam1 and $estimparam2)
  infectproba <- calcInfectProba(data = data, formula = formula, namepara1 = name1, namepara2 = name2, nameage = nameage, w1 = w1, w2 = w2, mort = mort, a = a, v0para1 = v0para1, v0para2 = v0para2)
  
  tabth <- colSums(infectproba$matprobainfect) ## tabth is the theoretical table based on real data 
  tabth <- rbind(tabth[1:2], tabth[3:4])
  
  ## Observed corrected chi-square:
  chi2corrobs <- round(sum((tabobs - tabth) ^ 2 / tabth), digits = 5)
  tabth <- round(tabth, digits = 2)
  colnames(tabth) <- colnames(tabobs)
  rownames(tabth) <- rownames(tabobs) 
  if(!is.numeric(chi2corrobs))
    stop("fatal error in calculation of chi-square between the real data and the theoretical infectious status of each real individual")
  
  return(list(infectproba = infectproba, chi2corrobs = chi2corrobs, tabth = tabth, tabobs = tabobs))
}



simudata_chi2corrage <- function(formula, data, name1, name2, nameage, w1, w2, mort, a, v0para1, v0para2, matprobainfect) {
  ## simulation of infectious status for two parasites
  multinomS <- apply(matprobainfect, 1, function(X) {which(rmultinom(1, 1, X) == 1)})
  sero2 <- (multinomS + 1) %% 2
  sero1 <- (multinomS - sero2 - 1) / 2
  
  ## construction of simulated data set
  sero1 <- as.factor(sero1)
  sero2 <- as.factor(sero2)
  levels(sero1) <- c("0", "1")
  levels(sero2) <- c("0", "1")
  dataS <- data ## dataS is Simulated data
  dataS[[name1]] <- sero1
  dataS[[name2]] <- sero2
  
  infectprobaS <- calcInfectProba(data = dataS, formula = formula, namepara1 = name1, namepara2 = name2, nameage = nameage, w1 = w1, w2 = w2, mort = mort, a = a, v0para1 = v0para1, v0para2 = v0para2) 
  tabth <- colSums(infectprobaS$matprobainfect) ## Theoretical table based on Simulated data
  tabobs <- as.vector(t(table(sero1, sero2))) ## Observed table based on Simulated data
  return(sum((tabth - tabobs) ^ 2 / (tabth)))
}


