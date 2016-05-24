# gld package for R

# Fit the FKML gld to data via:
# 1. Maximum Likelihood (ML)
# 2. Maximum Spacings Product (MSP) (MPS)
# 3. Titterington's Method (TM)
# 4. Starship Method (SM)
# 5. Method of TL-moments (TL)
# 6. Distributional Least Absolute (DLA)

# Written by Benjamin Dean (24/09/2010) 
# Changes 2013,2014 Robert King
# Example use:
# x <- rlogis(100)
# fit.fkml(x, method="ML")

########################################################################
#                                 Fit data                             #
########################################################################

fit.fkml <- function(x, method="ML", t1=0, t2=0, 
  l3.grid=c(-0.9,-0.5,-0.1,0,0.1,0.2,0.4,0.8,1,1.5), l4.grid=l3.grid,
  record.cpu.time = TRUE,optim.method="Nelder-Mead",
  inverse.eps = .Machine$double.eps, optim.control=list(maxit=10000), 
  optim.penalty=1e20, return.data=FALSE) {

# Start timing
if(record.cpu.time) { time.1 <- as.numeric(proc.time()[3]) 
} else { time.1 <- "not timing"
}
# Sample size
n <- length(x)
# Sort the data
x <- sort(x)

# Penalty value for optim
penalty.value <- optim.penalty

# Aliases
if (method == "ML")  {method.id <- 1; method.name="Maximum Likelihood"}
if (method == "MSP") {method.id <- 2; method.name="Maximum Spacings Product"}
if (method == "MPS") {method.id <- 2; method.name="Maximum Product of Spacings"}
if (method == "TM")  {method.id <- 3; method.name="Titterington's"}
if (method == "SM")  {method.id <- 4; method.name="Starship"}
if (method == "TL")  {method.id <- 5; method.name="Trimmed L-Moments"}
if (method == "DLA")  {method.id <- 6; method.name="Distributional Least Absolutes"}

# Perform the fitting process: 

if (method.id != 5) { # ML, MSP, TM, SM or DLA (not TL)

  # Starting values from grid search
  grid.results <- grid.search(l3.grid,l4.grid,method.id,x,n,
                              penalty.value, inverse.eps)

  pars.start <- grid.results$lambda
  # Optimise
  optim.results <- optim(par=pars.start,fn=obj.fnc,method=optim.method,
    control=optim.control,method.id=method.id,x=x,n=n,
    penalty.value=optim.penalty,inverse.eps=inverse.eps)
  # Add the optimisation method to the list
  optim.results$optim.method = optim.method
    
  # Optim results
  pars.est <- optim.results$par; conv <- optim.results$convergence
  its <- optim.results$counts[["function"]]; obj <- optim.results$value
  valid <- as.numeric(obj < penalty.value)
      
  # Store results
  if (record.cpu.time) {time.2 <- as.numeric(proc.time()[3]); runtime <- round(time.2-time.1,2) } else {runtime=NA}
  #
  old.style.results <- data.frame(l1=pars.est[1],l2=pars.est[2],
    l3=pars.est[3],l4=pars.est[4],valid=valid,conv=conv,
    its=its,runtime=runtime,
    obj=obj,l1.start=pars.start[1],l2.start=pars.start[2],
    l3.start=pars.start[3],l4.start=pars.start[4])
  result <- list(lambda=unclass(pars.est[1:4]),grid.results=grid.results,
      optim.results=optim.results,param="fkml",
      method.code=method,method.name=method.name,
      fkml.oldstyle.results=old.style.results)
  if (return.data) {result$data = x}
  class(result) <- "starship"
  names(result$lambda) <- paste("lambda",1:length(result$lambda),sep="")
} else { # TL

  # Starting values from grid search
  grid.results <- grid.search.tl(l3.grid,l4.grid,method.id,x,n,
    penalty.value,t1,t2,inverse.eps)
    # note that grid.min$lambda has only two elements, l3 & l4

  pars.start <- grid.results$lambda
  # Optimise
  optim.results <- optim(par=pars.start,fn=obj.fnc.tl,
        method=optim.method,control=list(maxit=10000),
        method.id=method.id,x=x,n=n,penalty.value=penalty.value,
        t1=t1,t2=t2,inverse.eps=inverse.eps)
     
  # Optim results
  l3.est <- optim.results$par[1]; l4.est <- optim.results$par[2]
  conv <- optim.results$convergence; its <- optim.results$counts[["function"]]
  obj <- optim.results$value
  # Add the optimisation method to the list
  optim.results$optim.method = optim.method 
  
  # Rerun C code at l3.est and l4.est to get l1.est, l2.est and valid  
  tl.final <- .C("fit_fkml",as.integer(5),as.double(0),as.double(1),
    as.double(l3.est),as.double(l4.est),as.double(x),as.integer(n),
    as.double(penalty.value),as.double(t1),as.double(t2),
    l1.tl=as.double(0),l2.tl=as.double(1),valid.tl=as.integer(0),
    as.double(0),as.double(inverse.eps)) 
  l1.est <- tl.final[["l1.tl"]]; l2.est <- tl.final[["l2.tl"]]
  valid <- (obj < penalty.value) * tl.final[["valid.tl"]]
    
  # Store results
if (record.cpu.time) {time.2 <- as.numeric(proc.time()[3]); runtime <- round(time.2-time.1,2)} else {runtime <- NA}
  old.style.results <- data.frame(l1=l1.est,l2=l2.est,l3=l3.est,
    l4=l4.est,valid=valid,conv=conv,its=its,runtime=runtime,obj=obj,
    l3.start=pars.start[1],l4.start=pars.start[2])
  result <- list(lambda=unclass(c(l1.est,l2.est,l3.est,l4.est)),
    grid.results=grid.results,optim.results=optim.results,param="fkml",
    method.code=method,method.name=method.name,
    fkml.oldstyle.results=old.style.results)
  if (return.data) {result$data = x}
  class(result) <- "starship"
  names(result$lambda) <- paste("lambda",1:length(result$lambda),sep="")
}

# Return result
return(result)
}

################################################################################
#                             Objective function                               #
################################################################################

obj.fnc <- function(par,method.id,x,n,penalty.value,inverse.eps) {

l1 <- par[1]
l2 <- par[2]
l3 <- par[3]
l4 <- par[4]

obj.result <- .C("fit_fkml",as.integer(method.id),as.double(l1),
  as.double(l2),as.double(l3),as.double(l4),as.double(x),
  as.integer(n),as.double(penalty.value),as.double(0),as.double(0),
  as.double(0),as.double(1),as.integer(0),result = as.double(0),
  as.double(inverse.eps))

return(obj.result[["result"]])
}

################################################################################
#                  Objective function (method of TL-moments)                   #
################################################################################

# In the method of TL-moments, the objective function only depends on l3 and l4.
# This is different to the other methods (which depend on all four parameters),
# so the objective function must be constructed differently.

obj.fnc.tl <- function(par,method.id,x,n,penalty.value,t1,t2,inverse.eps) {

l3 <- par[1]
l4 <- par[2]

obj.result.tl <- .C("fit_fkml",as.integer(method.id),as.double(0),
  as.double(1),as.double(l3),as.double(l4),as.double(x),
  as.integer(n),as.double(penalty.value),as.double(t1),as.double(t2),
  as.double(0),as.double(1),as.integer(0),result = as.double(0),
  as.double(inverse.eps))

return(obj.result.tl[["result"]])
}

################################################################################
#                      Grid search for starting values                         #
################################################################################

grid.search <- function(l3.grid,l4.grid,method.id,x,n,penalty.value,inverse.eps) {

# Initialise variables
l1.s <- 0; l2.s <- 1; l3.s <- -0.1; l4.s <- -0.1; obj.min <- Inf

# Loop through grid values and record obj
for (i in 1:length(l3.grid)) {
  for (j in 1:length(l4.grid)) {
    l3 <- l3.grid[i]; l4 <- l4.grid[j]
    l2 <- (qgl(0.75,0,1,l3,l4) - qgl(0.25,0,1,l3,l4))/IQR(x)
    l1 <- median(x) - qgl(0.5,0,l2,l3,l4)
    obj <- obj.fnc(c(l1,l2,l3,l4),method.id,x,n,penalty.value,inverse.eps)      
    if (obj < obj.min) {
      l1.s <- l1; l2.s <- l2; l3.s <- l3; l4.s <- l4; obj.min <- obj
    }
  }
}

return(list(response=obj.min,lambda=c(l1.s,l2.s,l3.s,l4.s)))
}

################################################################################
#            Grid search for starting values (method of TL-moments)            #
################################################################################

# A different grid search is also needed for the method of TL-moments.

grid.search.tl <- function(l3.grid,l4.grid,method.id,x,n,penalty.value,t1,t2,inverse.eps) {

# Initialise variables
l3.s <- -0.1; l4.s <- -0.1; obj.min <- Inf

# Loop through grid values and record obj
for (i in 1:length(l3.grid)) {
  for (j in 1:length(l4.grid)) {
    l3 <- l3.grid[i]; l4 <- l4.grid[j]
    obj <- obj.fnc.tl(c(l3,l4),method.id,x,n,penalty.value,t1,t2,inverse.eps)      
    if (obj < obj.min) {
      l3.s <- l3; l4.s <- l4; obj.min <- obj
    }
  }
}

return(list(response=obj.min,lambda=c(l3.s,l4.s)))
}
