#########################################################
### Cubic splines estimation method for 'couponbonds' ###
#########################################################

estim_cs <- function(bonddata, group, matrange="all",rse=TRUE) UseMethod("estim_cs")

### Cubic spline term structure estimation 
estim_cs.couponbonds <- function(bonddata, group, matrange="all",rse=TRUE) {

  ## data preprocessing 
  prepro <- prepro_bond(group=group,bonddata=bonddata,
           matrange=matrange)


  n_group=prepro$n_group
  sgroup=prepro$sgroup
  cf=prepro$cf
  cf_p=prepro$cf_p
  m=prepro$m
  m_p=prepro$m_p
  p=prepro$p
  ac=prepro$ac
  y=prepro$y
  duration=prepro$duration

    
  # Choosing knot points (McCulloch)
  # number of bonds in each group 
  K <- mapply(function(k) ncol(m[[k]]),sgroup,SIMPLIFY=FALSE)  
  # number of basis functions
  s <-  mapply(function(k) round(sqrt(K[[k]])),sgroup,SIMPLIFY=FALSE)
  
  # only perfom spline estimation if number of bonds per group >= 9
  if(sum(s>=3) != length(s))  stop(cat("Estimation aborted:
    For cubic splines estimation more than 9 observations per group are required","\n",
    "Check group(s):", group[which((s>=3)==FALSE)]),
    "\n" )
  
  # only used for knot point finding
  i <- mapply(function(k) 2:(max(2,(s[[k]]-2))),sgroup,SIMPLIFY=FALSE)  
  
  h <-  mapply(function(k) trunc(((i[[k]]-1)*K[[k]])/(s[[k]]-2)),sgroup,SIMPLIFY=FALSE)
             
  theta <- mapply(function(k)((i[[k]]-1)*K[[k]])/(s[[k]]-2)-h[[k]],sgroup,SIMPLIFY=FALSE)
  
  # calculate knot points
  T <- mapply(function(k) if(s[[k]]>3) c(floor(min(y[[k]][,1])),
       apply(as.matrix(m[[k]][,h[[k]]]),2,max)
       + theta[[k]]*(apply(as.matrix(m[[k]][,h[[k]]+1]),2,max)-apply(as.matrix(m[[k]][,h[[k]]]),2,max)),
       max(m[[k]][,ncol(m[[k]])])) else c(floor(min(y[[k]][,1])),max(m[[k]][,ncol(m[[k]])])),sgroup,SIMPLIFY=FALSE)
 
  # parameter estimation with OLS
  # dependent variable
  Y <- mapply(function(k) apply(cf_p[[k]],2,sum),sgroup,SIMPLIFY=FALSE)
  

  # k ... group index
  # sidx ... index for spline function  
  # independetn variable	
  X <- mapply(function(k) mapply(function(sidx)  apply(cf[[k]]*mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]])),2,sum), 1:s[[k]] ),sgroup,SIMPLIFY=FALSE)

  # OLS parameter estimation
  regout <- mapply(function(k) lm(-Y[[k]]~X[[k]]-1),sgroup,SIMPLIFY=FALSE) # parameter vector
  
  # estimated paramters  
  alpha <- lapply(regout, coef)
  
  # calculate discount factor matrix 
  dt <- list()
   for (k in sgroup){
   dt[[k]] <- matrix(1,nrow(m[[k]]),ncol(m[[k]]))
   for(sidx in 1:s[[k]]){
    dt[[k]] <- dt[[k]] + alpha[[k]][sidx]* mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]]))
   }
  }  

  # calculate estimated prices 
  phat <- mapply(function(k) apply(cf[[k]]*dt[[k]],2,sum),sgroup,SIMPLIFY=FALSE)
  
  # price errors
  perrors <- mapply(function(k) cbind(y[[k]][,1],phat[[k]] - p[[k]]),sgroup,SIMPLIFY=FALSE)     
  for (k in sgroup) class(perrors[[k]]) <- "error"
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(-phat[[k]],cf[[k]]),m_p[[k]]),sgroup,SIMPLIFY=FALSE)
  
  # yield errors
  yerrors <- mapply(function(k) cbind(y[[k]][,1], yhat[[k]][,2] - y[[k]][,2]),sgroup,SIMPLIFY=FALSE)
  for (k in sgroup) class(yerrors[[k]]) <- "error"
  
  # maturity interval
  t <- mapply(function(k) seq(max(round(min(T[[k]]),2),0.01),max(T[[k]]),0.01), sgroup,SIMPLIFY=FALSE)  
 
  # calculate mean and variance of the distribution of the discount function 
  mean_d <- mapply(function(k) apply(mapply(function(sidx) alpha[[k]][sidx]*
  			gi(t[[k]],T[[k]],sidx,s[[k]]),1:s[[k]]),1,sum) +1, sgroup, SIMPLIFY=FALSE)

  #  variance covariance matrix for estimated parameters
  if(rse) Sigma <- lapply(regout,vcovHAC.default) else Sigma <- lapply(regout,vcov) 
  
  var_d <- mapply(function(k) apply(mapply(function(sidx) gi(t[[k]],T[[k]],
  			sidx,s[[k]]),1:s[[k]]),1,function(x) t(x)%*%Sigma[[k]]%*%x), sgroup, SIMPLIFY=FALSE) 
  
  # lower 95% confidence interval
  cl <- mapply(function(k) mean_d[[k]] + rep(qt(0.025,nrow(X[[k]])- ncol(X[[k]])),
                    length(mean_d[[k]]))*sqrt(var_d[[k]]),sgroup,SIMPLIFY=FALSE)	
  # upper 95 % confidence interval	
  cu <- mapply(function(k) mean_d[[k]] + rep(qt(0.975,nrow(X[[k]])- ncol(X[[k]])),
                    length(mean_d[[k]]))*sqrt(var_d[[k]]),sgroup,SIMPLIFY=FALSE) 
  
  # zero cupon yield curves for maturity interval t 
  zcy_curves <-  mapply(function(k)  cbind(t[[k]],-log(mean_d[[k]])/t[[k]],-log(cl[[k]])/t[[k]],
 				 -log(cu[[k]])/t[[k]]),sgroup, SIMPLIFY=FALSE )  
 
   
  for (k in sgroup) class(zcy_curves[[k]]) <- "ir_curve"
  class(zcy_curves) <- "spot_curves"
  
  # calculate spread curves           	    
   if(n_group != 1) {     
   srange <- seq(max(unlist(lapply(t,min))),min(unlist(lapply(t,max))),0.01)
   s_curves <- mapply(function(k) cbind(srange,zcy_curves[[k]][c(which(zcy_curves[[k]][,1]== srange[1]): which(zcy_curves[[k]][,1] == srange[length(srange)])),2] 
    - zcy_curves[[1]][c(which(zcy_curves[[1]][,1]== srange[1]): which(zcy_curves[[1]][,1]== srange[length(srange)])),2]),sgroup, SIMPLIFY=FALSE) 
   
   } else s_curves = "none"
   for (k in sgroup) class(s_curves[[k]]) <- "ir_curve" 
   class(s_curves) <- "s_curves"
  
  # create discount factor curves 
  df_curves <- mapply(function(k) cbind(t[[k]],mean_d[[k]]),sgroup,SIMPLIFY=FALSE)
  
  for (k in sgroup) class(df_curves[[k]]) <- "ir_curve"
  class(df_curves) <- "df_curves"	
  
  # calculate forward rate curves
  
  fwr_curves <- mapply(function(k) cbind(t[[k]],impl_fwr(zcy_curves[[k]][,1],zcy_curves[[k]][,2])),sgroup,SIMPLIFY=FALSE)
  
   for (k in sgroup) class(fwr_curves[[k]]) <- "ir_curve"
  class(fwr_curves) <- "fwr_curves"	

 # return list of results
 result <- list(  group=group,          # e.g. countries, rating classes
                  matrange=matrange,    # maturity range of bonds
                  n_group=n_group,      # number of groups
                  knotpoints=T,         # knot points
                  rse=rse,              # robust standard errors
                  spot=zcy_curves, 	# zero coupon yield curves
                  spread=s_curves,      # spread curves
                  discount=df_curves,	# forward rate curves
                  forward=fwr_curves,	# discount factor curves
                  cf=cf,                # cashflow matrix
                  m=m,                  # maturity matrix
                  p=p,                  # dirty prices
                  phat=phat,            # estimated prices
                  perrors=perrors,     	# price errors
                  y=y,                  # maturities and yields
                  yhat=yhat,            # estimated yields
                  yerrors=yerrors,	# yield errors
                  alpha=alpha,          # cubic splines parameters                             
                  regout=regout         # OLS output
                  
                 )
                 
  # assign names to results list 
  for ( i in 6:length(result)) names(result[[i]]) <- group 
    
  class(result) <- "termstrc_cs"
  result
 
}

  




