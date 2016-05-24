methmage <-
function(X, Y){
  ## CONSTRUCT THE GEE MODEL AND REPORT TEST STATISTIC AND PVALUE FROM THE GEEGLM TEST
  
  NofSites <- ncol(X);
  n1 <- nrow(X);
  n2 <- nrow(Y);
  
  Data <- list();
  Data$X <- matrix(rbind(X, Y),byrow=FALSE);
  Data$Y <- factor(c(rep(0,n1*NofSites),rep(1,n2*NofSites)));
  A <- rep(seq(1,n1+n2),rep(NofSites,n1+n2));

  ## Calculate the three test statistics - PA, CQ and SD
  Model <- geeglm(X ~ Y, data = Data, id = A, corstr = "ar1",control=geese.control(maxit=100,epsilon=1e-8));
  R <- c(summary(Model)$coefficients[3][[1]][2], summary(Model)$coefficients[4][[1]][2]);

  return(R)
}
