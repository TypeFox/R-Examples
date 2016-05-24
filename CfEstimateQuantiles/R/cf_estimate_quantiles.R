cf_estimate_quantiles <-
function(cf_expansion_order, probabilities, cumulants) {
    
  cf<-function(count, x, a) {
    
    get<-function(poly, h) {
      result = 0;
      for (i in 1:length(poly)) {
        result = result + h[i]*poly[i]
      }
      return(result)
    }
    
    prod<-function(vec1, vec2) {
      n = length(vec1)
      res = vec1*0
      for (i1 in 1:n) {
        for (i2 in 1:n) {
          if (i1 + i2 <= n) {
            res[i1+i2] = res[i1+i2] + vec1[i1]*vec2[i2]      
          }
        }    
      }  
      return(res)
    }
    
    xi_h = (1:(count*count*3)) * 0
    dim(xi_h) = c(count, count*3)
    xi = (1:count) * 0
    h = 1:(count*3)
    h[1] = x;
    h[2] = x*x - 1;
    
    poly_h = c(1:((3*count)^2)) * 0;
    dim(poly_h) = c(3*count, 3*count)
    
    #Hermite
    for (i in 3:(count*3)) {
      h[i] = +(x * h[i-1] - (i-1) * h[i-2])
    }
    
    for (i in 1:(count*3)) {
      poly_h[i,i] = 1;
    }
    
    result =x
    
    for (k in 1:count) {
      xi_h[k,] = a[k]*poly_h[k+1,]
      if (k-1 != 0) 
        for (j in 1:(k-1)) {
          q = 1:(count*3) * 0
          
          q11 = prod(xi_h[k-j,]*xi[j], poly_h[1,])
          q12 = prod(xi_h[k-j,]*a[j], poly_h[j+2,])
          q21 = xi[k-j]*xi[j]*poly_h[1,]
          q22 = xi[k-j] * a[j] * poly_h[j+2,];
          
          q = q11 - q12 - q21 + q22;
          xi_h[k,] = xi_h[k,] - (j/k) * (q)
        }
      xi[k] = get(xi_h[k,], h)
      result = result + xi[k]
    }
    return(result)
  }  
  
  ## main body 
  
  cf_expansion_order = as.numeric(cf_expansion_order)
  probabilities = as.numeric(probabilities)
  cumulants = as.numeric(cumulants)
  m = cumulants[1];
  s2 = cumulants[2];
  a = rep(0, cf_expansion_order-2)
  for (i in 3:length(cumulants)) {
    a[i-2] = cumulants[i] / factorial(i) / (sqrt(s2) ^ (i))
  }
  
  
  result = rep(0, length(probabilities))
  
  for (i in 1:length(probabilities))
  {
    if (cf_expansion_order == 2)
      result[i] =  (m + 1/2 * (s2 + 1)*qnorm(probabilities[i]))    
    else     
      result[i] = ((m +  cf(cf_expansion_order - 2, qnorm(probabilities[i]), a)*sqrt(s2)))    
  }  
  return(result)
}

