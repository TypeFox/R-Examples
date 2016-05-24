# private method
# Runge-Kutta method for solving differential equations. It is used to generate both Lorenz and
# Rossler systems.
rungeKutta=function(func,initial.condition,time,params){
  n.samples = length(time)
  h = time[[2]]-time[[1]]
  y = matrix(ncol=length(initial.condition), nrow=n.samples)
  y[1,] = initial.condition
  for (i in 2:n.samples){
    k1 = h*func(y[i-1,],time[[i-1]],params)
    k2 = h*func(y[i-1,] + k1/2 , time[[i-1]] + h/2, params)
    k3 = h*func(y[i-1,] + k2/2 , time[[i-1]] + h/2, params)
    k4 = h*func(y[i-1,] + k3 , time[[i-1]] + h, params)

    y[i,] = y[i-1,] + (k1 + 2*k2 + 2*k3 + k4)/6
  }
  
  return (y)
}

# private method
# Trapezoidal rule for numerical integration 
trapezoidalRule = function(x, integrand ){
  index = 2:length(x)
  par(mfrow=c(1,1))
  plot(x,integrand,'l')
  return (as.double( (x[index] - x[index-1]) %*% (integrand[index] + integrand[index-1])) / 2)
}
  

# private method implementing (y(x+h)-y(x-h))/2h 
differentiate = function(h,y){
  len = length(y)
  if (len >= 3){
    derivative = (y[3:len]-y[1:(len-2)])/(2*h)  
  }else{
    # if not possible... use (y(x+h)-y(x))/h
    derivative = diff(y)/(h)  
  }
  
  return(derivative)
}
  
differentiateAxis = function(x){
  len = length(x)
  if (len >= 3){
    # We have used the (y(x+h)-y(x-h))/2h  rule ...
    # Eliminate first and last 
    axis = x[-c(1,len)]
  }else{
    # We have used the (y(x+h)-y(x))/h rule ...
    # Eliminate last
    axis = x[-c(len)]
  }
  
  return(axis)
}