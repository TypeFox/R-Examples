sim.pw.ar2 <-
function(N,sd_u,b.slope,b.slope2,br.loc) {
  num.breaks = length(br.loc)
  phi = c()
  phi2 = c()
  y = rep(0,N)
  for (i in 1:length(b.slope)) {
    if (i == 1)   {
      phi=c(phi,rep(b.slope[1],br.loc[1]))
      phi2=c(phi2,rep(b.slope2[1],br.loc[1]))
      next
    }
    if (i == length(b.slope)) {
      phi=c(phi,rep(tail(b.slope,1),N - br.loc[length(b.slope)-1]))
      phi2=c(phi2,rep(tail(b.slope2,1),N - br.loc[length(b.slope2)-1]))
      break
    }
    else {
      phi=c(phi,rep(b.slope[i],br.loc[i] - br.loc[i-1]))
      phi2=c(phi2,rep(b.slope2[i],br.loc[i] - br.loc[i-1]))
    }
  }
  y[1] =  rnorm(1,0,sd_u)
  y[2] =  rnorm(1,0,sd_u)
  for (i in 3:N) {
    y[i] = phi[i] * y[i-1] + phi2[i] * y[i-2] + rnorm(1,0,sd_u)
  }
  output = list(NULL)
  output[[1]] = phi
  output[[3]] = phi2
  output[[2]] = y
  output[[4]] = br.loc
  return(output)
}
