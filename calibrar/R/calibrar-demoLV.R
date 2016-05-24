# T = 100
# par = list(r=0.5, l=0.2, K=100, alpha=0.1, gamma=0.1, initial=list(N=10, P=1))


.generatePredatorPreyModel = function(path, r=0.5, l=0.2, alpha=0.1, gamma=0.1, K=100, T=100, 
                                      N0=10, P0=1, ...) {

  # 'real' parameters
  par_real = list(r=r, l=l, K=K, alpha=alpha, gamma=gamma, initial=list(N=N0, P=P0))
  
  pop = .PredatorPreyModel(par=par_real, T=T)
  
  # observed abundances
  n = rapply(pop, f=jitter, how = "list") 
  
  main.folder   = file.path(path, "PredatorPreyDemo")
  data.folder   = file.path(main.folder, "data")
  
  if(!file.exists(data.folder)) dir.create(data.folder, recursive=TRUE)
    
  for(i in c("prey", "predator")) {
    ifile = paste0(i, ".csv")
    dat = matrix(n[[i]], ncol=1)
    colnames(dat) = i
    write.csv(dat, file.path(data.folder, ifile))
  }
  
  # parInfo.csv
  
  parInfo = list()
  parInfo$guess = list(r=0.1, l=0.1, K=1.1*max(n$prey), alpha=0.05, gamma=0.1, initial=list(N=n$prey[1], P=n$predator[1]))
  parInfo$lower = list(r=0, l=0, K=0.25*max(n$prey), alpha=0, gamma=0, initial=list(N=0.5*n$prey[1], P=0.5*n$predator[1]))
  parInfo$upper = list(r=2, l=2, K=5*max(n$prey), alpha=1, gamma=1, initial=list(N=1.5*n$prey[1], P=1.5*n$predator[1]))
  parInfo$phase = list(r=1, l=1, K=1, alpha=1, gamma=1, initial=list(N=NA, P=NA))
  
  # calibrationInfo.csv
  
  calibrationInfo = list()
  calibrationInfo$variable  = c("prey", "predator")
  calibrationInfo$type      = "lnorm2"
  calibrationInfo$calibrate = TRUE
  calibrationInfo$weights   = 1
  calibrationInfo$useData   = TRUE
  
  calibrationInfo = as.data.frame(calibrationInfo)
  
  write.csv(calibrationInfo, file.path(main.folder, "calibrationInfo.csv"), row.names=FALSE)
  
  constants = list(T=T)
  
  output = c(list(path=main.folder, par=par_real), constants, parInfo)
  
  return(output)
  
}


.PredatorPreyModel = function(par, T) {
  if(!requireNamespace("deSolve", quietly = TRUE)) 
    stop("You need to install the 'deSolve' package.")# check on hadley
  # par is a list with 'alpha', 'beta' 'gamma', 'sd' and 'mu_ini'.
  LV = function(t, y, parms, ...) {
    r = parms$r
    l = parms$l
    alpha = parms$alpha
    gamma = parms$gamma
    K = parms$K
    dN = r*y[1]*(1-(y[1]/K)) - alpha*y[1]*y[2]
    dP = -l*y[2] + gamma*alpha*y[1]*y[2]
    return(list(c(dN, dP)))
  }
  times = seq(0, T)
  y0 = c(par$initial$N, par$initial$P)
  sol = deSolve::ode(y=y0, times=times, func=LV, parms=par, method="ode45")
  out = as.list(as.data.frame(sol[,-1]))
  names(out) = c("prey", "predator")
  out$prey[is.na(out$prey)] = 0
  out$predator[is.na(out$predator)] = 0
  return(out)
}

