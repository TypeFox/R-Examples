eegsim <-
  function(channel,time,coefs=rep(1,5),tshift=rep(0,5)){
    ###### Simulate event-related potential EEG data
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: February 16, 2015
    
    if(length(time)!=length(channel)){stop("Inputs channel and time should be same length.")}
    if(length(coefs)!=5L){stop("Incorrect number of input coefficients.")}
    if(length(tshift)!=5L){stop("Incorrect number of input time shifts.")}
    cbind(p1s(channel)*p1t(time-tshift[1]),n1s(channel)*n1t(time-tshift[2]),p2s(channel)*p2t(time-tshift[3]),n2s(channel)*n2t(time-tshift[4]),p3s(channel)*p3t(time-tshift[5]))%*%coefs
    
  }

# temporal functions
p1t <- function(x){exp(-600*((x-0.10)^2))}
n1t <- function(x){-exp(-600*((x-0.17)^2))}
p2t <- function(x){exp(-450*((x-0.23)^2))}
n2t <- function(x){-exp(-450*((x-0.28)^2))}
p3t <- function(x){exp(-300*((x-0.35)^2))}

# spatial functions
p1s <- function(x){
  volts <- c(0.65, 0.42, 1.18, 1.15, 0.85, 0.74, 0.31, 1.15, 0.05, -0.5, 2, 1.69, 3.03, 2.4, 3.38, 2.84, 2.95, 2.24, 1.9, 0.96, 2.07, 1.38, 2.20, 1.8,  2.4, 2.22, 0.26, 1.37, 0.73, 0.84, 0.365, 0.420, 1.475, 1.120, 1.200, 1.110, 0.575, 0.025, -0.250)
  vnames <- c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CPZ", "O1", "O2", "OZ", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "PO1", "PO2", "PO3", "PO4", "PO5", "PO6", "PO7", "PO8", "POZ", "PZ", "TP7", "TP8", "TP9", "TP10", "P9", "P10", "PO9", "PO10", "I1", "I2", "IZ")
  vidx <- match(x,vnames)
  funval <- volts[vidx]
  funval[is.na(funval)] <- 0
  funval
}

n1s <- function(x){
  volts <- c(1.58, 2.04, 2.74, 2.57, 3.69, 3.13, 1.09, 8.68, 8.45, 8.44, 4.79, 5.04, 5.86, 4.73, 6.7, 5.94, 7.06, 5.54, 7.88, 7.67, 8.09, 7.69, 8.31, 7.72, 8.52, 7.75, 8, 3.42, 3.56, 2.52, 1.78, 1.26, 3.53, 2.77, 4.26, 3.875, 4.34, 4.225, 4.22)
  vnames <- c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CPZ", "O1", "O2", "OZ", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "PO1", "PO2", "PO3", "PO4", "PO5", "PO6", "PO7", "PO8", "POZ", "PZ", "TP7", "TP8", "TP9", "TP10", "P9", "P10", "PO9", "PO10", "I1", "I2", "IZ")
  vidx <- match(x,vnames)
  funval <- volts[vidx]
  funval[is.na(funval)] <- 0
  funval
}

p2s <- function(x){
  volts <- c(0.51, 0.47, 0.52, 1.13, -0.29, -0.05, 0.27, 0.91, 0.57, 0.7, 2.3, 2.35, 2.39, 2.24, 1.74, 1.68, -0.22, 0.03, 3.01, 3.02, 2.93, 2.75, 2.84, 2.48, 2.76, 2.21, 2.49, 2.26, -2.29, -2.07, -1.145, -1.035, -0.11, 0.015, 1.38, 1.105, 0.455, 0.285, 0.35)
  vnames <- c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CPZ", "O1", "O2", "OZ", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "PO1", "PO2", "PO3", "PO4", "PO5", "PO6","PO7", "PO8", "POZ", "PZ", "TP7", "TP8", "TP9", "TP10", "P9", "P10", "PO9", "PO10", "I1", "I2", "IZ")
  vidx <- match(x,vnames)
  funval <- volts[vidx]
  funval[is.na(funval)] <- 0
  funval
}

n2s <- function(x){
  volts <- c(-0.68, -0.45, -0.85, -0.99, -0.53, -1.16, -0.37, 1.08, 1.24, 1.54, -0.61, -0.11, -0.96, -1.11, -1.28, -1.44, -0.07, -1.24, 0.3, 0.04, 0.40, -0.04, 0.49, -0.11, 0.59, -0.19, 0.99, -0.05, 0.7, -0.7, 0.35, -0.35, -0.035, -0.62, 0.295, -0.095, 0.54, 0.62, 0.77)
  vnames <- c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CPZ", "O1", "O2", "OZ", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "PO1", "PO2", "PO3", "PO4", "PO5", "PO6","PO7", "PO8", "POZ", "PZ", "TP7", "TP8", "TP9", "TP10", "P9", "P10", "PO9", "PO10", "I1", "I2", "IZ")
  vidx <- match(x,vnames)
  funval <- volts[vidx]
  funval[is.na(funval)] <- 0
  funval
}

p3s <- function(x){
  volts <- c(2.08, 2.13, 2.62, 3.61, 2.57, 3.33, 1.32, 2.21, 2.04, 1.58, 3.66, 3.5, 4.53, 4.9, 4.79, 5.1, 3.96, 4.52, 3.68, 3.7, 3.81, 3.85, 3.94, 3.99, 4.07, 4.14, 2.93, 3.71, 2.28, 3.94, 1.14, 1.97, 1.98, 2.26, 2.035, 2.07, 1.105, 1.02, 0.79)
  vnames <- c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CPZ", "O1", "O2", "OZ", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "PO1", "PO2", "PO3", "PO4", "PO5", "PO6","PO7", "PO8", "POZ", "PZ", "TP7", "TP8", "TP9", "TP10", "P9", "P10", "PO9", "PO10", "I1", "I2", "IZ")
  vidx <- match(x,vnames)
  funval <- volts[vidx]
  funval[is.na(funval)] <- 0
  funval
}
