
#############################################
#Modifying the original cross validation function to fit our purpose#
#The original function forces R to display CV results from each step along with the graphic#
cvlm <-function (form.lm, data, m=10, seed = 178) {
  
  df <- data.frame(data)
  
  if (class(form.lm) == "formula") {
    form <- form.lm
  } else if (class(form.lm) %in% c("call", "lm")) {
    form <- formula(form.lm) 
  } else { 
    stop("form.lm must be formula or call or lm object") 
  }
  formtxt <- deparse(form)
  mf <- model.frame(form, data = df)
  ynam <- attr(mf, "names")[attr(attr(mf, "terms"), "response")]
  df.lm <- lm(mf)   
  tm <- terms(mf)
  xcolumns <- labels(tm)
  n <- nrow(df)
  df[, ynam] <- model.response(mf)
  df[, "Predicted"] <- predict(df.lm)
  df[, "cvpred"] <- numeric(n)
  yval <- mf[, ynam]
  if (!is.null(seed)) 
    set.seed(seed) 
  n <- dim(df)[1]
  rand <- sample(n)%%m + 1
  foldnum <- sort(unique(rand))
  
  sumss <- 0
  for (i in foldnum) {
    rows.in <- rand != i
    rows.out <- rand == i
    subs.lm <- lm(form, data = df[rows.in, ])
    df[rows.out, "cvpred"] <- suppressWarnings(predict(subs.lm, newdata = df[rows.out,]))
    resid <- df[rows.out, ynam] - df[rows.out, "cvpred"]
    ss <- sum(resid^2)
    sumss <- sumss + ss 
  }  
  sumdf <- sum(!is.na(df[, "Predicted"])) - length(xcolumns)
  sumres <- sumss/sumdf
  df <- sumdf
  outcv <- list(sumres = sumres, df = df, m = m)
}
#End of the modification### 

#Defining the function to generate the power series#
Generate.Powers <- function(X, lambda){
  
  FX <- data.frame(X)
  
  Full.Data<-do.call(polym,c(as.list(FX),degree=lambda,raw=TRUE))
  tab <- cor(Full.Data) == 1
  tab[upper.tri(tab, diag=1)] <- NA
  ind <- which(tab==1, arr.ind=TRUE)[,1] 
  if ( length(ind) > 0 ) {  Full.Data <- Full.Data[,-ind] }
  
  K.list<-rep(0,dim(Full.Data)[2])
  for(i in 1:length(K.list)){
    K.list[i]<-sum(as.numeric(unlist(strsplit(colnames(Full.Data)[i],split="[.]"))))
  }
  
  ranks<-rank(K.list,ties.method="first") 
  Sorted.Cov<-data.frame(Full.Data[,order(ranks)])
  return(Sorted.Cov)
  
}

### The main power series function
npse <- function(formula, order = 3, m = 10, seed = 178){
  
  ## extract response, covaraites, treatment, and instrument
  mf <- model.frame(formula = formula)
  Z <- model.response(mf)
  X <- model.matrix(attr(mf, "terms"), data = mf)[,-1]

  #Default maximum lambda is set to be 3#
  CV.Residuals<-rep(0,order)  

  for (i in 1:order){

    Sorted.Full.Data <- cbind(Z, Generate.Powers(X, i))   
    
    #Starting cross validation#
    Power.Regression <- lm(Z ~ .,data=as.data.frame(Sorted.Full.Data))    
    Power.CV <- cvlm(form.lm=Power.Regression, data=Sorted.Full.Data, m=m)
  
    #Storing the residuals# 
    CV.Residuals[i] <- Power.CV$sumres
    
  }
  
  #Choose the optimum lambda and estimate the tao#
  Lambda.Optimum<-which.min(CV.Residuals)
  Data.Optimum<-Generate.Powers(X, Lambda.Optimum)
  TaoHat <- lm(Z ~ .,data=as.data.frame(Data.Optimum))$fitted.values
  
  cvout <- list(fitted = TaoHat, Lambda = Lambda.Optimum, Data.opt = Data.Optimum, CV.Res = CV.Residuals)
  return(cvout)
}
