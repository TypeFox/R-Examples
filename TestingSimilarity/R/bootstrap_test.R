################################################################################
#' Bootstrap test for testing dose response curves for similarity concerning the maximum absolute deviation
#'
#'Function for testing whether two dose response curves can be assumed as equal concerning the
#'hypotheses \deqn{H_0: \max_{x\in\mathcal{X}} |m_1(d,\theta_1)-m_2(d,\theta_2)|\geq \epsilon\ vs.\ 
#'H_1: \max_{x\in\mathcal{X}} |m_1(d,\theta_1)-m_2(d,\theta_2)|< \epsilon.}
#'See \url{http://arxiv.org/pdf/1505.05266.pdf} for details.
#' 
#' @name bootstrap_test
#' @export
#' @importFrom lattice xyplot
#' @importFrom DoseFinding fitMod
#' @importFrom DoseFinding defBnds
#' @importFrom alabama auglag
#' @importFrom stats coef 
#' @importFrom stats ecdf 
#' @importFrom stats rnorm
#' @param data1,data2 data frame for each of the two groups containing the variables referenced in dose and resp  
#' @param m1,m2 model types. Built-in models are "linlog",  "linear",  "quadratic",  "emax",  "exponential",  "sigEmax",  "betaMod" and "logistic" 
#' @param epsilon positive argument specifying the hypotheses of the test
#' @param B number of bootstrap replications. If missing, default value of B is 5000
#' @param bnds1,bnds2 bounds for the non-linear model parameters. If not specified, they will be generated automatically 
#' @param plot if TRUE, a plot of the absolute difference curve of the two estimated models will be given 
#' @param scal,off fixed dose scaling/offset parameter for the Beta/ Linear in log model. If not specified, they are 1.2*max(dose) and 1 respectively
#' @return A list containing the p.value, the maximum absolute difference of the models, the estimated model parameters and the number of bootstrap replications. Furthermore plots of the two models are given.
#' @examples
#' library("DoseFinding")
#' library("alabama")
#' data(IBScovars)
#' male<-IBScovars[1:118,]
#' female<-IBScovars[119:369,]
#' bootstrap_test(male,female,"linear","emax",epsilon=0.35,B=300) 
#' @references \url{http://arxiv.org/pdf/1505.05266.pdf}
################################################################################  
library('DoseFinding')
library('alabama')

bootstrap_test <- function(data1,data2,m1,m2,epsilon,B=2000,bnds1=NULL,bnds2=NULL,plot=FALSE,scal=NULL,off=NULL) {
  
  #wrong specified data
  if (is.null(data1$dose)){
    stop("The data is not referenced in dose and resp")
  }
  
  x1 <- data1$dose;
  x2 <- data2$dose;
  y1 <- data1$resp;
  y2 <- data2$resp;
  
  y <- seq(min(x1,x2),max(x1,x2),length=500); #creating a grid for searching for the maximum
  
  #missing models or incorrect epsilon
  if (missing(m1) | missing(m2)) {
    stop("Need to specify the models that should be fitted")
  };
  if (epsilon<0) {
    stop("epsilon has to be positive")
  };
  
  
  #fixing the scaling parameter for the beta model
  if (m1=="betaMod" & is.null(scal)) {
    scal=1.2*max(x1)
  }
  if (m2=="betaMod" & is.null(scal)) {
    scal=1.2*max(x2)
  };
  
  #fixing the offset parameter for the linear in log model
  if ((m1=="linlog" | m2=="linlog") & is.null(off)) {
    off=1
  };
  
  builtIn <- c("linlog", "linear", "quadratic", "emax", 
               "exponential", "logistic", "betaMod", "sigEmax");
  mods <- c(linlog=function(d,e){linlog(d,e,off=off)},linear,quadratic,emax,exponential,logistic,betaMod=function(d, e){betaMod(d, e, scal=scal)},sigEmax);
  modelNum1 <- match(m1, builtIn);
  modelNum2 <- match(m2, builtIn);
  
  #wrong specification of the model
  if (is.na(modelNum1) | is.na(modelNum2)) {
    stop("Invalid dose-response model specified")
  };
  
 
  #fix the model type for calculating the test statistic
  m1f <- mods[[modelNum1]]; 
  m2f <- mods[[modelNum2]];
  
  #placeholders for the functions defining parameter bounds
  hin1 <- function(x){1};
  hin2 <- function(x){1}; 

  
  #fitting of the first model
  if (is.null(bnds1)) {
    bnds1=defBnds(max(x1))[[m1]]
  }
  mod1 <- fitMod (x1,y1, model=m1, bnds=bnds1);
  nop1 <- length(coef(mod1)) #number of parameters of the first model  
  
  #fitting of the second model
  if (is.null(bnds2)) {
    bnds2=defBnds(max(x2))[[m2]]
  }
  mod2 <- fitMod (x2,y2, model=m2, bnds=bnds2);
  nop2 <- length(coef(mod2)) #number of parameters of the second model

  #constructing a function defining the bounds for the model parameters for the optimization under nullhypothesis 
  if (m1=="emax" | m1=="exponential"){
    hin1 <- function(x) {
      h <- rep(NA, 1)
      h[1] <- x[3]-bnds1[1]
      h[2] <- -x[3]+bnds1[2]
      h
    }
  }else if (m1=="logistic" | m1=="sigEmax" | m1=="betaMod"){
    hin1 <- function(x) {
      h <- rep(NA, 1)
      h[1] <- x[3]-bnds1[1,1]
      h[2] <- -x[3]+bnds1[1,2]
      h[3] <- x[4]-bnds1[2,1]
      h[4] <- -x[4]+bnds1[2,2]
      h
    }
  };
  
  #same for the second model
  if (m2=="emax" | m2=="exponential"){
    hin2 <- function(x) {
      h <- rep(NA, 1)
      h[1] <- x[3]-bnds2[1]
      h[2] <- -x[3]+bnds2[2]
      h
    }
  }else if (m2=="logistic" | m2=="sigEmax" | m2=="betaMod"){
    hin2 <- function(x) {
      h <- rep(NA, 1)
      h[1] <- x[3]-bnds2[1,1]
      h[2] <- -x[3]+bnds2[1,2]
      h[3] <- x[4]-bnds2[2,1]
      h[4] <- -x[4]+bnds2[2,2]
      h
    }
  };
  
  coeff <- c(coef(mod1),coef(mod2)); #saving the coefficients of the models
  tstat <- max(dff(y,coef(mod1),coef(mod2),m1f,m2f),na.rm=TRUE);#calculation of the test statistic 

  if (tstat>=epsilon){
    minimum <- coeff} else {
      loglikelihood <- function(v){
        alpha <- v[1:nop1]; beta <- v[(nop1+1):length(v)];
        sum((y1-m1f(x1,alpha))^2)+sum((y2-m2f(x2,beta))^2)
      };
      softmax <- function(v){
        alpha <- v[1:nop1]; beta <- v[(nop1+1):length(v)];
        sum(dff(y,alpha,beta,m1f,m2f)*exp(50*dff(y,alpha,beta,m1f,m2f)))/(sum(exp(50*dff(y,alpha,beta,m1f,m2f))))-epsilon
      };
      hin <-function(v){
        alpha <- v[1:nop1]; beta <- v[(nop1+1):length(v)];
        h <- rep(NA,1)
        h[1:length(hin1(c(rep(1,4))))] <- hin1(alpha)
        h[(length(hin1(c(rep(1,4))))+1):(length(hin1(c(rep(1,4))))+length(hin2(c(rep(1,4)))))] <- hin2(beta)
        h
      }
      minimum <- auglag(par=coeff,loglikelihood,hin=hin,heq=softmax,control.outer=list(trace=0))$par
    }
  
  tstern <- rep(NA,1);#vector for the bootstrap test statistics

  #bootstrap 
  for (i in 1:B)
  {
    y1stern <- m1f(x1,minimum[1:nop1])+rnorm(length(x1),sd=sqrt(mod1$RSS/mod1$df));
    y2stern <- m2f(x2,minimum[(nop1+1):length(minimum)])+rnorm(length(x2),sd=sqrt(mod2$RSS/mod2$df));#Erzeugung der Bootstrap-Daten
    #fitting the first bootstrap model
    mod1stern <- fitMod (x1,y1stern, model=m1, bnds=bnds1);
    #fitting the second bootstrap model
    mod2stern <- fitMod (x2,y2stern, model=m2, bnds=bnds2)
    tstern[i] <- max(dff(y,coef(mod1stern),coef(mod2stern),m1f,m2f));#calculating the bootstrap test statistic
  }

  fn <- ecdf(tstern);#calculating the p-value by using the ecdf of the bootstrap distribution

  #small B
  if (B<200) {
    warning("Warning: A larger B should be choosen for higher accuracy of the test.")
  }
  
  #print plots of the two estimated models
  #return p-value and the test statistic (the maximum absolute deviation of the curves)
  #print plot, if desired
  plot1<-plot(mod1,ylab="Response",xlab="Dose",main="Model 1");
  plot2<-plot(mod2,ylab="Response",xlab="Dose",main="Model 2");
  if(plot==TRUE) {
    plot3<-xyplot(dff(y,coef(mod1),coef(mod2),m1f,m2f)~y,type=c("l","g"),xlab="Dose",ylab="Response",main="Absolute difference curve of the two fitted models",col="black")
    print(plot2, position=c(0.5,0.5, 1, 1), more=TRUE)
    print(plot1, position=c(0, 0.5, 0.5, 1),more=TRUE)
    print(plot3, position=c(0.1,0,0.9,0.5))
  }
  else {
    print(plot2, position=c(0.5,0.1, 1, 0.9), more=TRUE)
    print(plot1, position=c(0, 0.1, 0.5, 0.9))
  }
  return(list(p.value=fn(tstat),max.abs.difference=tstat,bootstrap.replications=B,estimated.model.param.model1=coeff[1:nop1],estimated.model.param.model2=coeff[(nop1+1):(nop1+nop2)]))
  
}