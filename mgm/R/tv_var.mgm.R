

tv_var.mgm <- function(data, # input tv mgm
                      type, 
                      lev,
                      lags = 1, 
                      tsteps = 'default', 
                      bandwidth,
                      gam = .25,                       
                      d = 2, 
                      rule.reg = 'AND', 
                      pbar = TRUE,
                      method ='glm',
                      missings = 'error',
                      ret.warn = FALSE) {
  
  n <- nrow(data)
  nNode <- ncol(data)
  
  # create time sequence
  tpoints <- seq(0,1, length=(n-1)) # n sequence before and after for folded gaussian weight

  if(tsteps=='default') {tsteps <- round(2/bandwidth)} # reasonable coveres the data for any choice of bandwidth
  est_steps <- seq(0, 1, length = tsteps) # times at which estimation takes places
  
  # storage
  a_wadj <- array(NA, dim=c(nNode,nNode,tsteps))
  l_mpar.matrix <- list()
  
  #progress bar
  if(pbar==TRUE) {pb <- txtProgressBar(min = 0, max=tsteps, initial=0, char="-", style = 3) }
  
  # fit weighted mgm model
  tv_list <- list() # storage
  for(t in 1:tsteps) {
    
    kernel_mean <- est_steps[t]
    
    # define weights with folded gaussian kernel    
    weights <- dnorm(tpoints, mean = kernel_mean, sd=bandwidth) # we use unit time interval
    weights <- weights / max(weights) 
    
    # estimate using mgm
    tv_list[[t]] <- mgmfit_core(data, type, lev, lambda.sel="EBIC", folds = 10, gam = gam,  # force EBIC, CV makes no sense in nonstationary data
                                d=d, rule.reg = rule.reg, pbar = FALSE, method = method, missings = missings, 
                                weights=weights, ret.warn=FALSE, VAR=TRUE)
    
    a_wadj[,,t] <- tv_list[[t]]$wadj
    l_mpar.matrix[[t]] <- tv_list[[t]]$mpar.matrix
    
    if(pbar==TRUE) { setTxtProgressBar(pb, t) } # update progress bar
    
  }
  
  # write model parameter matrix in array
  a_mpar.matrix <- array(NA, dim=c(dim(l_mpar.matrix[[t]])[1], dim(l_mpar.matrix[[t]])[2],tsteps))
  for(t in 1:tsteps) a_mpar.matrix[,,t] <- l_mpar.matrix[[t]]
  
  # get weighted N for each time step; this will later also be used as a missing data signal
  l_wN <- unlist(lapply(tv_list, function(x) sum(x$call$weights)))
  
  # outputlist
  l_call <- list('type'=type, 'lev'=lev, 'tsteps'=tsteps, 'bandwidth'=bandwidth, 'lambda.sel'="EBIC",
                 'gam'=gam, 'd'=d, 'rule.reg'=rule.reg, 'method'=method)
  outlist <- list('call'=l_call, 'wadj'=a_wadj, 'mpar.matrix'=a_mpar.matrix, 't.models'=tv_list, 'Nt'=l_wN)
  
  return(outlist) 
  
}
