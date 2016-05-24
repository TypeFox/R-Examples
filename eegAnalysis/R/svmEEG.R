svmEEG <-
function(x, method="C-classification",scale = TRUE, kernel =
      "radial", degree = 3, gamma = if (is.vector(x)) 1 else 1 / ncol(x),
    coef0 = 0, cost = 1, nu = 0.5,
    class.weights = NULL, cachesize = 40, tolerance = 0.001, epsilon = 0.1,
    shrinking = TRUE, cross = 0, probability = TRUE, fitted = TRUE, seed = 1L, subset, na.action = na.omit)
{
  dados<-x$FinalFea
  if (is.null(nrow(dados))) dados<-as.matrix(dados)
  nfea<-nrow(dados)
  gamma<-1/nrow(dados)
  type = NULL

  
  dados<-as.data.frame(t(dados))
  y<-x$label
  dados<-cbind(dados,y)
  pars<-tune.svm(y~.,data=dados, gamma = (gamma*c(0.5,1,2,10)), cost = c(0.1,1,5,10))
  bestcost<-pars$best.parameters$cost
  bestgamma<-pars$best.parameters$gamma

  
  model<- svm(t(x$FinalFea), x$label ,method= method,kernel=kernel,cost=bestcost,
  gdamma=bestgamma,coef0=coef0,degree=degree,probability=probability,type=type,scale=scale,
  nu=nu, class.weights=class.weights, cachesize=cachesize, tolerance=tolerance,
  epsilon=epsilon, shrinking=shrinking, cross=cross, fitted=fitted, seed=seed, subset=subset,
  na.action=na.action)
  
  result<-list(model=model, ncomps=x$ncomps, W=x$W, featype=x$featype,
  win=x$win,stat=x$stat,power=x$power,abs=x$abs,log=x$log,mintomax=x$mintomax,nfea=nfea,nch=x$nch,
  which.classes=x$which.classes,L=x$L,wavelet=x$wavelet,variance=x$variance)
  
  class(result)<-"svmEEG"
  
  return(result)
}



print.svmEEG <-
  function(x,...){
    print(x$model)
  }


summary.svmEEG <-
  function(object,...){
    x<-object
    print(x$model)
  }
