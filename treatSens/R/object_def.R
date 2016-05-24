############
#combine.sensitivity
#combines the results in two sensitivity objects into single object for plotting
############
combine.sensitivity <- function(x1, x2){
  if(!class(x1) == "sensitivity" | !class(x2) == "sensitivity")
    stop("Objects must be of class \"sensitivity\"")
  if(!all.equal(c(x1$Y, x1$Z), c(x2$Y, x2$Z)))
    stop("Combining objects only allowed for sensitivity analyses on same data set")
  if(!(x1$sensParam == x2$sensParam))
    stop("Combining objects only allowed for sensitivity analyses using same type of sensitivity parameter")
  if(!(x1$model.type == x2$model.type))
    stop("Combining objects only allowed for sensitivity analyses using same model")

  grid1 = expand.grid(as.numeric(dimnames(x1$sp.y)[[1]]), as.numeric(dimnames(x1$sp.z)[[2]]))
  gridpts1 = matrix(c(rep(grid1[,1], dim(x1$tau)[3]), rep(grid1[,2], dim(x1$tau)[3])), ncol = 2)
  tauLong1 = c(x1$tau)
  seTauLong1 = c(x1$se.tau)
  sigRespLong1 <- c(x1$sig2.resp)
  sigTrtLong1  <- c(x1$sig2.trt)

  grid2 = expand.grid(as.numeric(dimnames(x2$sp.y)[[1]]), as.numeric(dimnames(x2$sp.z)[[2]]))
  gridpts2 = matrix(c(rep(grid2[,1], dim(x2$tau)[3]), rep(grid2[,2], dim(x2$tau)[3])), ncol = 2)
  tauLong2 = c(x2$tau)
  seTauLong2 = c(x2$se.tau)
  sigRespLong2 <- c(x2$sig2.resp)
  sigTrtLong2  <- c(x2$sig2.trt)

  result <- list(model.type = x1$model.type, sensParam = x1$sensParam, tau = c(tauLong1, tauLong2), se.tau = c(seTauLong1, seTauLong2), 
                 sp.z = c(gridpts1[,2], gridpts2[,2]), sp.y = c(gridpts1[,1], gridpts2[,1]), 
                 Y = x1$Y, Z = x1$Z, sig2.resp = c(sigRespLong1, sigRespLong2), sig2.trt = c(sigTrtLong1, sigTrtLong2),
            #should any of the following be averaged across objects?
                 tau0 = x1$tau0, se.tau0 = x1$se.tau0,
                 Xcoef = x1$Xcoef, Xcoef.plot = x1$Xcoef.plot,
                 varnames = x1$varnames,var_ytilde = x1$var_ytilde,var_ztilde = x1$var_ztilde, XpartCor = x1$XpartCor)
  class(result) <- "sensitivityCombo"
  return(result)
}

############
#summary.sensitivity
#Prints average value of tau (trt effect) for each cell in grid
#Default to dimensions labeled with target sensitivity parameters
#commented lines will relabel with average realized s.p.
############
summary.sensitivity <- function(object, ...){
  summary.sensitivity.default(object, ...)
}

summary.sensitivity.default <- function(object, digits = 3, signif.level = 0.05,...){
  Tau <- object$tau
  table <- round(apply(Tau, c(1,2), mean),digits)
  taus <- apply(object$tau, c(1,2), mean)
  part.cors = object$sensParam == "cor"
  
  if(part.cors){
	trt.coef <- as.numeric(dimnames(taus)[[2]])
	resp.coef <- as.numeric(dimnames(taus)[[1]])
	print.text <- "Partial correlations with U"
  }else {
  	resp.coef <- as.numeric(dimnames(taus)[[1]])
  	trt.coef <- as.numeric(dimnames(taus)[[2]])
	print.text <- "Coefficients on U" 
  }
  
  zeroCoords = contourLines(trt.coef, resp.coef, taus, levels = 0)
  if(class(unlist(zeroCoords))=="NULL") {
    cat("Sensitivity parameters where tau = 0 could not be calculated.\n\n")
  } else {    
    zeroCors = round(cbind(zeroCoords[[1]]$y, zeroCoords[[1]]$x),digits)
    colnames(zeroCors) <- c("Y", "Z")
    rownames(zeroCors) <- rep("", dim(zeroCors)[1])
    cat(print.text, "where tau = 0:\n")
    print(zeroCors)
    cat("\n\n")
  } 
  
  K = dim(object$se.tau)[3]
  W = apply(object$se.tau^2, c(1,2), mean, na.rm = T)
  B = apply(object$tau, c(1,2), sd, na.rm = T)^2
  se.est = t(sqrt(W+(1+1/K)*B))
  
  noSigCoords = contourLines(trt.coef, resp.coef, taus/t(se.est), levels = -sign(object$tau0)*qnorm(signif.level/2))
  if(class(unlist(noSigCoords))=="NULL") {
    cat("Sensitivity parameters where significance level", signif.level, "is lost could not be calculated.\n\n")
  } else {
    noSigCors = round(cbind(noSigCoords[[1]]$y, noSigCoords[[1]]$x),digits)
    colnames(noSigCors) <- c("Y", "Z")
    rownames(noSigCors) <- rep("", dim(noSigCors)[1])
    cat(print.text, "where significance level", signif.level, "is lost:\n")
    print(noSigCors)
    cat("\n\n")
  }
  
  colnames(table) <- trt.coef
  rownames(table) <- resp.coef
  cat("Estimated treatment effects\n")
  print(table)
}


############
#summary.sensitivityCombo
#Prints average value of tau (trt effect) for each cell in grid
#Default to dimensions labeled with target sensitivity parameters
#commented lines will relabel with average realized s.p.
############
summary.sensitivityCombo <- function(object, ...){
  summary.sensitivityCombo.default(object, ...)
}

summary.sensitivityCombo.default <- function(object, digits = 3, signif.level = 0.05,...){
  Zcors = object$sp.z #horizontal grids of U
  Ycors = object$sp.y #vertical grids of U
  taus = object$tau
  taus = tapply(taus, list(Zcors, Ycors), mean, na.rm = T)
  table <- round(taus,digits)
  part.cors = object$sensParam == "cor"
  
  if(part.cors){
    trt.coef <- as.numeric(dimnames(taus)[[2]])
    resp.coef <- as.numeric(dimnames(taus)[[1]])
    print.text <- "Partial correlations with U"
  }else {
    resp.coef <- as.numeric(dimnames(taus)[[1]])
    trt.coef <- as.numeric(dimnames(taus)[[2]])
    print.text <- "Coefficients on U" 
  }
  
  zeroCoords = contourLines(trt.coef, resp.coef, taus, levels = 0)
  if(class(unlist(zeroCoords))=="NULL") {
    cat("Sensitivity parameters where tau = 0 could not be calculated.\n\n")
  } else {    
    zeroCors = round(cbind(zeroCoords[[1]]$y, zeroCoords[[1]]$x),digits)
    colnames(zeroCors) <- c("Y", "Z")
    rownames(zeroCors) <- rep("", dim(zeroCors)[1])
    cat(print.text, "where tau = 0:\n")
    print(zeroCors)
    cat("\n\n")
  } 
  
  W = tapply(object$se.tau^2, list(Zcors, Ycors), mean, na.rm = T)
  K = length(W)/length(Ycors)
  B = tapply(object$tau, list(Zcors, Ycors), sd, na.rm = T)^2
  se.est = sqrt(W+(1+1/K)*B)
    
  noSigCoords = contourLines(trt.coef, resp.coef, taus/se.est, levels = -sign(object$tau0)*qnorm(signif.level/2))
  if(class(unlist(noSigCoords))=="NULL") {
    cat("Sensitivity parameters where significance level", signif.level, "is lost could not be calculated.\n\n")
  } else {
    noSigCors = round(cbind(noSigCoords[[1]]$y, noSigCoords[[1]]$x),digits)
    colnames(noSigCors) <- c("Y", "Z")
    rownames(noSigCors) <- rep("", dim(noSigCors)[1])
    cat(print.text, "where significance level", signif.level, "is lost:\n")
    print(noSigCors)
    cat("\n\n")
  }
  
  colnames(table) <- trt.coef
  rownames(table) <- resp.coef
  cat("Estimated treatment effects\n")
  print(table)
}



##############
#print.sensitivity 
#prints average values for each cell in grid of:
#tau, SE of tau, (realized sens params,) coefficients and their se's
##############
print.sensitivity = function(x, digits=3, part.cors = F, ...){
  if(part.cors){
	trt.coef <- as.numeric(dimnames(x$tau)[[2]])/sqrt(x$var_ztilde)
	resp.coef <- as.numeric(dimnames(x$tau)[[1]])/sqrt(x$var_ytilde) *(1-trt.coef^2)
	print.text <- "Partial correlations with U"
  }else {
  	resp.coef <- as.numeric(dimnames(x$tau)[[1]])
  	trt.coef <- as.numeric(dimnames(x$tau)[[2]])
	print.text <- "Coefficients on U" 
  }

  
  table <- round(apply(x$tau, c(1,2), mean),digits)
  colnames(table) <- trt.coef
  rownames(table) <- resp.coef
  cat("Estimated treatment effects\n")
  print(table, ...)
  
  K = dim(x$se.tau)[3]
  W = apply(x$se.tau^2, c(1,2), mean)
  B = apply(x$tau, c(1,2), sd)^2
  table <- round(sqrt(W+(1+1/K)*B),digits)
  #table <- round(apply(x$se.tau, c(1,2), mean),digits)
  colnames(table) <- trt.coef
  rownames(table) <- resp.coef
  cat("Standard error of estimated treatment effects\n")
  print(table, ...)
  
  if(!part.cors){
  	table <- round(apply(x$sp.y, c(1,2), mean),digits)
  	colnames(table) <- trt.coef
  	rownames(table) <- resp.coef
  	cat("Estimated zeta.y - coefficient of U in response model\n")
  	print(table, ...)
  
  	K = dim(x$se.spy)[3]
  	W = apply(x$se.spy^2, c(1,2), mean)
  	B = apply(x$sp.y, c(1,2), sd)^2
  	table <- round(sqrt(W+(1+1/K)*B),digits)
  	#table <- round(apply(x$se.zy, c(1,2), mean),digits)
  	colnames(table) <- trt.coef
  	rownames(table) <- resp.coef
  	cat("Standard error of zeta.y\n")
  	print(table, ...)
  
  	table <- round(apply(x$sp.z, c(1,2), mean),digits)
  	colnames(table) <- trt.coef
  	rownames(table) <- resp.coef
  	cat("Estimated zeta.z - coefficient of U in treatment model\n")
  	print(table, ...)
  
  	K = dim(x$se.spz)[3]
  	W = apply(x$se.spz^2, c(1,2), mean)
  	B = apply(x$sp.z, c(1,2), sd)^2
  	table <- round(sqrt(W+(1+1/K)*B),digits)
  	#table <- round(apply(x$se.zz, c(1,2), mean),digits)
  	colnames(table) <- trt.coef
  	rownames(table) <- resp.coef
  	cat("Standard error of zeta.z\n")
  	print(table, ...)
  }
  invisible(x)
}


#############
#plot.sensitivity
#wrapper function - plot function in separate file
#############

#setMethod("plot", c("sensitivity", "missing"),
#  definition = function(x,y,...){
plot.sensitivity = function(x,y,...){
    sensPlot(x,...)
  }
plot.sensitivityCombo = function(x,y,...){
  sensPlot(x,...)
}  
#}
#)
