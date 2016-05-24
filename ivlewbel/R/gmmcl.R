gmmcl <-
function(formula1, formula2, data, cluster){
  # library(plyr) 
  # library(gmm)
  # create data.frame
  data$id1 = 1:dim(data)[1]
  formula3 = paste(as.character(formula1)[3],"id1", sep=" + ")
  formula4 = paste(as.character(formula1)[2], formula3, sep=" ~ ")
  formula4 = as.formula(formula4)
  formula5 = paste(as.character(formula2)[2],"id1", sep=" + ")
  formula6 = paste(" ~ ", formula5, sep=" ")
  formula6 = as.formula(formula6)
  frame1 = model.frame(formula4, data)
  frame2 = model.frame(formula6, data)
  dat1 = join(data, frame1, type="inner", match="first")
  dat2 = join(dat1, frame2, type="inner", match="first")
  
  # matrix of instruments
  Z1 = model.matrix(formula2, dat2)
  
  # step 1
  gmm1 = gmm(formula1, formula2, data = dat2, vcov="iid")
  
  # clustering weight matrix
  cluster = factor(dat2[,cluster])
  u = residuals(gmm1)
  estfun = sweep(Z1, MARGIN=1, u,'*')
  u = apply(estfun, 2, function(x) tapply(x, cluster, sum))  
  S = 1/(length(residuals(gmm1)))*crossprod(u)
  
  # step 2
  gmm2 = gmm(formula1, formula2, data=dat2, 
             vcov="TrueFixed", weightsMatrix = chol2inv(chol(S)))
  return(gmm2)
}
