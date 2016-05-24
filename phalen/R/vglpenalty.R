vglpenalty <-
function(X,x0 = NA,x1 = NA,plotpenalty = TRUE,allowed.error = 0.005,invert = FALSE) {
  # check for errors in passed values
  if (length(x0)!=length(x1)) {
    stop("x0 and x1 must be equal lengths") 
  } else if (length(x0)>2) {
    stop("vector x1 must be contain only 1 or 2 numbers")
  } else if (length(x0)==2) {
    if (all(c(x0[1],x1[1],x1[2],x0[2]) == sort(c(x0[1],x1[1],x1[2],x0[2])))==FALSE) {
      stop("x thresholds not in order, reorder so that x0[1] < x1[1] < x1[2] < x0[2]")
    } else if (x0[1]==x1[1] || x0[2]==x1[2]) {
      stop("x0 and x1 cannot be the same value")
    }
  } else if (x0[1] == x1[1]) {
    stop("x0 and x1 cannot be the same value")
  } else if (allowed.error <= 0 || allowed.error >= 1) {
    stop("allowed error must be greater than 0 and less than 1")
  }
  
  if (length(x0) == 1) {
    M = mean(c(x0,x1))  
    if (x0 < x1) {
      b = -(log((1/allowed.error)-1)/(x0-M))
      vglpenalty = 1/(1+exp(-b*(X-M)))
    } else {
      b = -(log((1/allowed.error)-1)/(x1-M))
      vglpenalty = 1-(1/(1+exp(-b*(X-M))))
    }
  }
  
  if (length(x0)==2) {
    M1 = mean(c(x0[1],x1[1]))
    M2 = mean(c(x0[2],x1[2]))
    b1 = -(log((1/allowed.error)-1)/(x0[1]-M1))
    b2 = -(log((1/allowed.error)-1)/(x1[2]-M2))
    brk = mean(c(x1[1],x1[2]))
    
    vglpenalty = ifelse(X<brk,1/(1+exp(-b1*(X-M1))),1 - (1/(1+exp(-b2*(X-M2)))))
  }
  
  if (invert==TRUE) {
    vglpenalty = 1-vglpenalty
  }
  
  if (plotpenalty == TRUE) {
    plot(vglpenalty,type='l',lwd = 2,col = "#4E4F52",ylim = c(0,1),ylab="glpenalty")
  }
  
  vglpenalty
}
