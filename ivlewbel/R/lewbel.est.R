lewbel.est <-
function(formula, data, clustervar = NULL, robust = TRUE){
  # library(gmm)
  # remove missing data
  mns = unlist(strsplit(as.character(formula)[c(2,3)],"\\|"))  
  if(length(mns)<4){
    stop("Error in formula")
  }
  mns = gsub("factor\\(", "", mns)
  mns = gsub("\\)", "", mns)
  mns = gsub(" ", "", mns)
  mns = unlist(strsplit(mns, "\\+"))
  mns = unique(unlist(strsplit(mns,":")))
  mns = unique(unlist(strsplit(mns,"\\*")))
  mns = gsub("\\)", "", mns)
  mns = gsub("\\n", "", mns)
  mns = gsub(" ", "", mns)
  mns = unlist(strsplit(mns, "\\+"))
  
  mns = unique(mns[mns!="NULL"])
  if(!is.null(clustervar)){mns = c(mns, clustervar)}
  
  inc.obs = complete.cases(data[,mns])
  data = data[inc.obs, ]
  
  # what variables 
  outc = as.character(formula)[2]  
  regs = unlist(strsplit(as.character(formula)[3],"\\|"))  
  endg = regs[1]
  endg = gsub("\\n", "", endg)  
  exog = regs[2]
  exog = gsub("\\n", "", exog)
  hetv = regs[3]
  hetv = gsub("\\n", "", hetv)
  
  if(length(grep("NULL", as.character(regs[4])))!=0){
    regs = regs[-4]
  } 
  if(!is.na(regs[4])){
    inst = regs[4]
    inst = gsub("\\n", "", inst)
    } else {
      inst = NULL
    }
  
  # should we exclude any for singularities?
  ex.form1 = as.formula(paste(outc, paste(endg, exog, sep=" + "), sep=" ~ "))
  ex.var = NULL
  ex.var = coef(lm(ex.form1, data))
  ex.var = names(ex.var[is.na(ex.var)])
  
  # outcome
  Y2 = model.matrix(as.formula(paste("~ -1", outc, sep=" + ")), data = data)
  colnames(Y2) = "y1"
  # endog regressors
  Y1 = model.matrix(as.formula(paste("~ -1", endg, sep=" + ")), data = data)
  colnames(Y1) = paste("y2", 1:ncol(Y1), sep=".")
  # exog regressors
  X1 = model.matrix(as.formula(paste("~ ", exog, sep=" + ")), data = data)
  X1 = X1[,!colnames(X1) %in% ex.var]
  colnames(X1) = paste("x1", 1:ncol(X1), sep=".")
  # instruments (if necessary)
  if(!is.na(regs[4])){
    Z1 = model.matrix(as.formula(paste("~", inst, sep=" + ")), data = data)
    Z1 = matrix(Z1[,-which(colnames(Z1) %in% "(Intercept)")], nrow = nrow(Y2))
    colnames(Z1) = paste("z1", 1:ncol(Z1), sep=".")
  }  
  # heteroskedastic instruments  
  Z2 = model.matrix(as.formula(paste("~", hetv, sep=" + ")), data = data)
  Z2 = Z2[,!colnames(Z2) %in% ex.var]
  Z2 = matrix(Z2[,-which(colnames(Z2) %in% "(Intercept)")], nrow = nrow(data))
  Z2 = apply(Z2, 2, function(x){x-mean(x)} )
  colnames(Z2) = paste("z2", 1:ncol(Z2), sep=".")
  
  # get residual based IVs
  XF = X1
  if(!is.na(regs[4]) | length(grep("NULL",regs[4]))!=0){XF = cbind(X1, Z1)}
  E1 = matrix(residuals(lm(Y1 ~ XF - 1)), nrow = dim(data)[1])
  Z3 = matrix(apply(E1, 2, function(x){x*Z2}), nrow = dim(data)[1])
  colnames(Z3) = paste("het.z", 1:ncol(Z3), sep="")
  data.new = data.frame(Y2, Y1, XF, Z3)
  if(!is.null(clustervar)){data.new[[clustervar]] = data[,clustervar]}
  
  frm1 = paste(colnames(Y2), paste(colnames(Y1), collapse=" + "), sep=" ~ ")
  frm1 = paste(frm1, paste(colnames(X1), collapse=" + "), sep=" + ")
  frm2 = paste("", paste(colnames(X1), collapse=" + "), sep=" ~ ")
  if(!is.na(regs[4])){
    frm2 = paste(frm2, paste(colnames(Z1), collapse=" + "), sep=" + ")
  }
  frm2  = paste(frm2, paste(colnames(Z3), collapse=" + "), sep=" + ") 
  
  frm1 = paste(frm1, " -1" , sep="")
  frm2 = paste(frm2, " -1" , sep="")
  frm1 = as.formula(frm1)  
  frm2 = as.formula(frm2)  
    
  # efficient gmm coefficient estimates
  if(is.null(clustervar) & robust == FALSE){
    gmm1 = gmm(frm1, frm2, data = data.new, vcov="iid")
  }
  if(is.null(clustervar) & robust == TRUE){
    data.new$fakeid = as.character(1:nrow(data.new))
    gmm1 = suppressMessages(gmmcl(frm1, frm2, data = data.new, cluster = "fakeid"))
  }
  if(!is.null(clustervar)){
    gmm1 = suppressMessages(gmmcl(frm1, frm2, data = data.new, cluster = clustervar))
  }
  
  coef.names = c(colnames(model.matrix(as.formula(paste("~ -1", endg, sep=" + ")), data = data)),
                 colnames(model.matrix(as.formula(paste("~ ", exog, sep=" + ")), data = data)))
  # excluded
  coef.names = coef.names[!coef.names %in% ex.var]
  results = summary(gmm1)$coefficients[1:length(coef.names),]
  rownames(results) = coef.names
  
  # add excluded
  if(length(ex.var)>0){
    excl.results = t(as.matrix(results[1:length(ex.var),]))
    excl.results[,] = NA
    row.names(excl.results) = ex.var
    results = rbind(results, excl.results)
  }
  
  # now get f-tests
  ne = ncol(Y1)
  nx = ncol(X1)
  if(!exists("Z1")){
    Z1 = NULL
  }
  nz = ifelse(is.null(ncol(Z1)), 0, ncol(Z1))
  nh = ncol(Z2)
  
  data.new$fakeid = NULL
  f.test.stats = ftest(data.new, ne=ne, nx=nx, nz=nz, nh=nh, clustervar=clustervar, robust=robust)
  names(f.test.stats) = colnames(model.matrix(as.formula(paste("~ -1", endg, sep=" + ")), data = data)) 
  
  results = list(coefs = results, 
                 j.test = summary(gmm1)$stest, 
                 num.obs = nrow(data.new),
                 f.test.stats = f.test.stats)
  return(results)
}
