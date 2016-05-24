ftest <-
  function(data.new, ne, nx, nz, nh, robust, clustervar){
  # get endogenous variable projections 
  cols.1 = (1:ncol(data.new))[colnames(data.new)!=clustervar]
  if(is.null(clustervar)){
    data.temp = data.new
  } else {
    data.temp = data.new[,(1:ncol(data.new))[colnames(data.new)!=clustervar]]
  }
       
  # get projections
  for(i in 1:ne){
    d.temp = as.matrix(cbind(data.temp[,(i+1)],data.temp[,(2+ne):ncol(data.temp)]))
    r.temp = lm(d.temp[,1] ~ -1 + d.temp[,-1])
    data.temp[[paste("py2",i,sep=".")]] = fitted.values(r.temp)
  }
  
  f.test.stats = rep(NA, ne)
  
  # get residuals and run ftests
  for(i in 1:ne){
    # make formula
    lhs.1 = paste("y2",i,sep=".")
    rhs.1 = ifelse(ne>1, paste("py2", (1:ne)[!(1:ne %in% i)], sep=".", collapse=" + "), "")
    rhs.2 = paste("x1", 1:nx, sep=".", collapse=" + ")
    rhs.3 = ifelse(ne>1, paste(rhs.1, rhs.2, sep = " + "), rhs.2)
    form.1 = as.formula(paste(lhs.1, rhs.3, sep=" ~ -1 + "))
    # reg endog regressor on exog
    data.temp[["residuals"]] = residuals(lm(form.1, data.temp))
    # take residuals and run on IV
    rhs.4 = ifelse(nz>0, paste("z1", 1:nz, sep = ".", collapse=" + "), "")
    rhs.5 = paste("het.z", 1:nh, sep = "", collapse=" + ")
    rhs.6 = ifelse(nz>0, paste(rhs.4, rhs.5, sep = " + "), rhs.5)
    form.2 = as.formula(paste("residuals", rhs.6, sep=" ~ -1 + "))
    
    m1 = lm(form.2, data.temp)
    m2 = lm(residuals ~ -1 , data.temp)
    
    if(is.null(clustervar) & robust == FALSE){f1 = waldtest(m1, m2)}
    if(is.null(clustervar) & robust == TRUE){f1 = waldtest(m1, m2, vcov = vcovHC(m1, type = "HC0"))}
    if(!is.null(clustervar)){
      data.temp[[clustervar]] = data.new[[clustervar]]
      f1 = waldtest(m1, m2, vcov = clusterVCV(data.temp, m1, cluster1=clustervar))
    }    
    f.test.stats[i] = f1$F[2]
  }
  
  return(f.test.stats)
}