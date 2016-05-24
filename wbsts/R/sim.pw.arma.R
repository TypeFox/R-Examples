sim.pw.arma <-
function(N,sd_u,b.slope,b.slope2, mac,br.loc) {
  num.breaks = length(br.loc)
  num.coeff = length(b.slope)
  if (length(sd_u)==1) sd_u=rep(sd_u,num.coeff)
  phi = c()
  phi2 = c()
  ma = c()
  p.sd = c()
  y = rep(0,N)
  for (i in 1:num.coeff) {
    
    if (i == 1)   {      
      phi=c(phi,rep(b.slope[1],br.loc[1]))
      phi2=c(phi2,rep(b.slope2[1],br.loc[1]))
      ma=c(ma,rep(mac[1],br.loc[1]))
      p.sd = c(p.sd,rep(sd_u[1],br.loc[1]))
      next
    }
    if (i == length(b.slope)) {
      phi=c(phi,rep(tail(b.slope,1),N - br.loc[length(b.slope)-1]))
      phi2=c(phi2,rep(tail(b.slope2,1),N - br.loc[length(b.slope2)-1]))
      ma=c(ma,rep(tail(mac,1),N - br.loc[length(mac)-1]))
      p.sd=c(p.sd,rep(tail(sd_u,1),N - br.loc[length(sd_u)-1]))        
      break
    } else {
      phi=c(phi,rep(b.slope[i],br.loc[i] - br.loc[i-1]))
      phi2=c(phi2,rep(b.slope2[i],br.loc[i] - br.loc[i-1]))
      ma=c(ma,rep(mac[i],br.loc[i] - br.loc[i-1]))
      p.sd=c(p.sd,rep(sd_u[i],br.loc[i] - br.loc[i-1])) 
    }
  }
  y[1] =  rnorm(1,0,sd_u[1])
  y[2] =  rnorm(1,0,sd_u[1])
  error_lag.t = rnorm(1,0,sd_u[1])
  for (i in 3:N) {
    error_t = rnorm(1,0,p.sd[i])
    y[i] = phi[i] * y[i-1] + phi2[i] * y[i-2] + ma[i]*error_lag.t +  error_t
    error_lag.t=error_t
  }
  output = list(NULL)
  output[[1]] = phi
  output[[3]] = phi2
  output[[2]] = y
  output[[4]] = br.loc
  output[[5]] = ma
  return(output)
}
