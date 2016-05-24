`fitFun` <-
function(par, t=29, sigma_t= 0.3 * t, V=1:100)  {
  s<-vector()
  for(i in 1:length(V))
    s <- append(s,volEq7(t=t, sigma_t=sigma_t, r=par$r,
                         sigma_r=par$sigma_r, V=V[i]))
  ## constrain expected freq. sum to 1 
  s/sum(s) 
}

