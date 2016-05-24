#' Linear regression using CorReg's method, with variable selection.
#' @import Rcpp
#' @import Rmixmod
#' @import lars
#' @import glmnet
#' @import elasticnet
#' @import Matrix
#' @import mclust 
#' @import glmnet
#' @import MASS
#' @import methods
# ' @import parcor
# '  @import  clere
# ' @import spikeslab
#' @import  rpart
#' @import corrplot
#' @import mvtnorm
# ' @import rtkpp
#' @importFrom grDevices col2rgb gray rgb colors
#' @importFrom utils install.packages
#' @importFrom graphics abline arrows boxplot legend matplot par points rect text title
#' @importFrom stats AIC BIC aov as.formula chisq.test coef confint.default cor dnorm lm pf predict qnorm qt rbinom rgamma rmultinom rnorm rpois rstudent runif sd var
#' @useDynLib CorReg
#' @export
#' @param B The (d+1)xd matrix associated to Z and that contains the parameters of the sub-regressions
#' @param lambda (optional) parameter for elasticnet or ridge (quadratic penalty) if select="elasticnet" or "ridge".
#' @param X The data matrix (covariates) without the intercept
#' @param Y The response variable vector
#' @param Z The structure (adjacency matrix) between the covariates
#' @param compl (boolean) to decide if the complete modele is computed
#' @param expl (boolean) to decide if the explicative model is in the output
#' @param pred (boolean) to decide if the predictive model is computed
#' @param select selection method in ("lar","lasso","forward.stagewise","stepwise", "elasticnet", "NULL","ridge","adalasso","clere","spikeslab")
#' @param criterion the criterion used to compare the models
#' @param K the number of clusters for cross-validation
#' @param groupe a vector of integer to define the groups used for cross-validation (to obtain a reproductible result)
#' @param Amax the maximum number of non-zero coefficients in the final model
# ' @param returning boolean: second predictive step (selection on I1 knowing I2 coefficients)
#' @param X_test validation sample
#' @param Y_test response for the validation sample
#' @param intercept boolean. If FALSE intercept will be set to 0 in each model.
#' @param alpha Coefficients of the explicative model to coerce the predictive step. if not NULL explicative step is not computed.
#' @param g (optional) number of group of variables for clere if select="clere"
# ' @param compl2 boolean to compute regression (OLS only) upon [X_f,epsilon] instead of [X_f,X_r]
# ' @param explnew alternative estimation
#' @description Computes three regression models: Complete (regression on the wole dataset X), marginal (regression using only independant covariates: \code{X[,colSums(Z)==0]}) and plug-in (sequential regression based on the marginal model and then use redundant covariates by plug-in, 
#' with a regression on the residuals of the marginal model by the residuals of the sub-regressions). Each regression can be computed with variable selection (for example the lasso).
#' @return a list that contains:
#' \item{compl}{Results associated to the regression on X}
#' \item{expl}{Results associated to the marginal regression on explicative covariates (defined by colSums(Z)==0)}
#' \item{pred}{Results associated to the plug-in regression model.}
#' \item{compl$A}{Vector of the regression coefficients (the first is the intercept).(also have expl$A and pred$A) }
#' \item{compl$BIC}{BIC criterion associated to the model (also have expl$A and pred$A) }
#' \item{compl$AIC}{AIC criterion associated to the model (also have expl$A) }
#' \item{compl$CVMSE}{Cross-validated MSE associated to the model (also have expl$A) }
#Attention cette fonction degage une puissance phenomenale (it's over 9000!)
# correg<-function (X = NULL, Y = NULL, Z = NULL, B = NULL, compl = TRUE, expl = FALSE, pred = FALSE,
#                 select = "lar",
#                 criterion = c("MSE", "BIC"),
#                 X_test = NULL, Y_test = NULL, intercept = TRUE, 
#                 K = 10, groupe = NULL, Amax = NULL, lambda = 1,returning=FALSE,
#                 alpha=NULL,g=5,compl2=FALSE,explnew=FALSE) 
#' @examples
#' \dontrun{
#' require(CorReg)
#'    #dataset generation
#'    base=mixture_generator(n=15,p=10,ratio=0.4,tp1=1,tp2=1,tp3=1,positive=0.5,
#'                           R2Y=0.8,R2=0.9,scale=TRUE,max_compl=3,lambda=1)

#'    X_appr=base$X_appr #learning sample
#'    Y_appr=base$Y_appr #response variable for the learning sample
#'    Y_test=base$Y_test #responsee variable for the validation sample
#'    X_test=base$X_test #validation sample
#'    TrueZ=base$Z#True generative structure (binary adjacency matrix)
#'    
#'    #Regression coefficients estimation
#'     select="lar"#variable selection with lasso (using lar algorithm)
#'    resY=correg(X=X_appr,Y=Y_appr,Z=TrueZ,compl=TRUE,expl=TRUE,pred=TRUE,
#'                select=select,K=10)
#'    
#'    #MSE computation
#'    MSE_complete=MSE_loc(Y=Y_test,X=X_test,A=resY$compl$A)#classical model on X
#'    MSE_marginal=MSE_loc(Y=Y_test,X=X_test,A=resY$expl$A)#reduced model without correlations
#'    MSE_plugin=MSE_loc(Y=Y_test,X=X_test,A=resY$pred$A)#plug-in model
#'    MSE_true=MSE_loc(Y=Y_test,X=X_test,A=base$A)# True model
#'    
#'    
#'    #MSE comparison
#'    MSE=data.frame(MSE_complete,MSE_marginal,MSE_plugin,MSE_true)
#'    MSE#estimated structure
#'    compZ$true_left;compZ$false_left
#'    
#'   barplot(as.matrix(MSE),main="MSE on validation dataset", sub=paste("select=",select))
#'   abline(h=MSE_complete,col="red")
#'    }



   correg<-function (X = NULL, Y = NULL, Z = NULL, B = NULL, compl = TRUE, expl = FALSE, pred = FALSE,
                     select = "lar",
                     criterion = c("MSE", "BIC"),
                     X_test = NULL, Y_test = NULL, intercept = TRUE, 
                     K = 10, groupe = NULL, Amax = NULL, lambda = 1,
                     alpha=NULL,  g=5) 
      
   
   {
      compl2=FALSE
      returning=FALSE
      explnew=FALSE
    
  if(is.null(X)){
     dat<- data.frame(t=seq(0, 2*pi, by=0.1) )
     xhrt <- function(t) 16*sin(t)^3
     yhrt <- function(t) 13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t)
     dat$y=yhrt(dat$t)
     dat$x=xhrt(dat$t)
     with(dat, plot(x,y, type="l"))  
     with(dat, polygon(x,y, col="red"))   
     points(c(10,-10, -15, 15), c(-10, -10, 10, 10), pch=169, font=5)
     title(main="I Love CorReg !")
     return("I Love CorReg !")
  }
  OLS=OLS
   res = list()
  X = 1*as.matrix(X)
  K = abs(K)
  K = min(K, nrow(X))
  Y=1*as.matrix(Y)
  if (is.null(groupe)) {
    groupe = rep(0:(K - 1), length.out = nrow(as.matrix(X)))
    groupe = sample(groupe)
  }
  select = select[1]
  if(select=="NULL"){
    returning=FALSE
  }
  if(select=="adalasso"){ requireNamespace("parcor")}
  criterion = criterion[1]
  if (is.null(Amax)) {
    Amax = ncol(X) + 1
  }
  if (sum(Z) == 0) {
    compl = TRUE
  }
  if(!is.null(alpha)){
  }
  if (compl) {
    if (select == "NULL") {
      res$compl$A = c(OLS(X = X, Y = Y, intercept = intercept)$beta)
    }else if (select != "elasticnet" & select != "ridge" & select != "adalasso" & select!="clere" & select!="spikeslab") {
      lars_compl = lars(x = X, y = Y, type = select, intercept = intercept)
      res$compl = meilleur_lars(lars = lars_compl, X = X, 
                                Y = Y, mode = criterion, intercept = intercept, 
                                K = K, groupe = groupe, Amax = Amax)
    }else if (select=="elasticnet"){
      lars_compl = renet(x = X, y = Y, intercept = intercept, 
                        lambda = lambda)
      names(lars_compl)[4] = "coefficients"
      res$compl = meilleur_lars(lars = lars_compl, X = X, 
                                Y = Y, mode = criterion, intercept = intercept, 
                                K = K, groupe = groupe, Amax = Amax)
    }else if(select=="adalasso"){
       resada=parcor::adalasso(X=X,y=Y,k=K)    
       if(intercept){
          if(is.null(resada$intercept.adalasso)){
             resada$intercept.adalasso=0
          }
          res$compl$A=c(resada$intercept.adalasso,resada$coefficients.adalasso)
       }else{
          res$compl$A=c(resada$coefficients.adalasso)
       }
       Xloc=X[,resada$coefficients.adalasso!=0]
       res$compl$A[res$compl$A!=0]=c(OLS(X=Xloc,Y=Y,intercept=intercept)$beta)
       res$compl$A=c(res$compl$A)
    }else if(select=="clere"){
       requireNamespace("clere")
       res$compl$A=A_clere(y=as.numeric(Y),x=X,g=g)
    }else if(select=="spikeslab"){
       respike=spikeslab::spikeslab(x=X,y=Y)
       res$compl$A=rep(0,times=ncol(X)+intercept)
       if(intercept){
          res$compl$A[c(1,1+which(respike$gnet.scale!=0))]=OLS(X=X[,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
       }else{
          res$compl$A[respike$gnet.scale!=0]=OLS(X=X[,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
       }
   }else{#ridge
      #res_ridge = linearRidge(Y~.,data=data.frame(X))
      res_ridge=glmnet(as.matrix(X),Y,alpha=0)
      res$compl$A=c(matrix(coef(res_ridge,lambda)))
    }
    res$compl$BIC = BicTheta(X = X, Y = Y, intercept = intercept, 
                             beta = res$compl$A)
    res$compl$AIC=mon_AIC(theta=res$compl$A,Y=Y,X=X,intercept=intercept) 
   res$compl$CVMSE = CVMSE(X = X, Y = Y, intercept = intercept, K = K, groupe = groupe)
  }
#explicatif
  if (sum(Z) != 0 & (expl | pred)) {
    qui = WhoIs(Z = Z)
    I1 = qui$I1
    I2 = qui$I2
    if(!is.null(alpha)){
       res$expl$A=alpha
    }else if (select == "NULL") {
      res$expl$A = OLS(X = as.matrix(X[, I1]), Y = Y, intercept = intercept)$beta
    }else if (select != "elasticnet" & select != "ridge" & select != "adalasso"  & select!="clere" & select!="spikeslab") {
      lars_expl = lars(x = as.matrix(X[, I1]), y = Y, type = select, 
                       intercept = intercept)
      res$expl = meilleur_lars(lars = lars_expl, X = as.matrix(X[, I1]), Y = Y, 
                               mode = criterion, intercept = intercept, 
                               K = K, groupe = groupe, Amax = Amax)
    }else if (select=="elasticnet"){
      lars_expl = renet(x = as.matrix(X[, I1]), y = Y, intercept = intercept, 
                       lambda = lambda)
      names(lars_expl)[4] = "coefficients"
      res$expl = meilleur_lars(lars = lars_expl, X = as.matrix(X[,I1]), Y = Y, 
                               mode = criterion, intercept = intercept, 
                               K = K, groupe = groupe, Amax = Amax)
    }else if(select=="adalasso"){
       resada=parcor::adalasso(X=X[,I1],y=Y,k=K)    
       if(intercept){
          if(is.null(resada$intercept.adalasso)){
             resada$intercept.adalasso=0
          }
          res$expl$A=c(resada$intercept.adalasso,resada$coefficients.adalasso)
       }else{
          res$expl$A=c(resada$coefficients.adalasso)
       }
       Xloc=X[,I1][,resada$coefficients.adalasso!=0]
       res$expl$A[res$expl$A!=0]=c(OLS(X=Xloc,Y=Y,intercept=intercept)$beta)       
    }else if(select=="clere"){
       res$expl$A=A_clere(y=as.numeric(Y),x=X[,I1],g=g)
    }else if(select=="spikeslab"){
       respike=spikeslab::spikeslab(x=X[,I1],y=Y)
       res$expl$A=rep(0,times=ncol(X[,I1])+intercept)
       if(intercept){
          res$expl$A[c(1,1+which(respike$gnet.scale!=0))]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
       }else{
          res$expl$A[respike$gnet.scale!=0]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
       }
    }else{#ridge
      #lars_expl = linearRidge(Y~.,data=data.frame(X[,I1]))
      res_ridge=glmnet(as.matrix(X[,I1]),Y,alpha=0)
      res$expl$A=c(matrix(coef(res_ridge,lambda)))
      #res$expl$A=coef(lars_expl)
    }
    if(is.null(alpha)){
       A_expl = rep(0, times = ncol(X) + intercept)
       if(intercept){
          A_expl[c(intercept, I1 + intercept)] = res$expl$A
       }else{
          A_expl[I1] = res$expl$A
       }
       res$expl$A = A_expl
    }else{
       A_expl=alpha
    }
    res$expl$BIC = BicTheta(X = X, Y = Y, intercept = intercept, beta = A_expl)
    res$expl$AIC=mon_AIC(theta=res$expl$A,Y=Y,X=X,intercept=intercept) 
    res$expl$CVMSE = CVMSE(X = X[,res$expl$A[-intercept]!=0], Y = Y, intercept = intercept, K = K, groupe = groupe)
    

#predictif
    if (pred & length(I2)>0) {     
      if (is.null(B)) {
        B = hatB(Z = Z, X = X)
      }
      Xtilde = X[, I2] - cbind(rep(1, times = nrow(X)), X[, I1]) %*% B[c(1, I1 + 1), I2]
      Xtilde = as.matrix(Xtilde)
      if(intercept){
        Ytilde = Y - as.matrix(X[, I1]) %*% A_expl[-1][I1] - A_expl[1]
      }else{
        Ytilde = Y - as.matrix(X[, I1]) %*% A_expl[I1]
      }
      if (select == "NULL") {
        A_inj = OLS(X = Xtilde, Y = Ytilde, intercept = F)$beta
      }else if (select != "elasticnet"  & select != "ridge" & select != "adalasso"  & select!="clere" & select!="spikeslab") {
        lars_inj = lars(x = Xtilde, y = Ytilde, type = select, 
                        intercept = F)
        if(max(lars_inj$R2)==0){
           A_inj=rep(0,times=ncol(Xtilde))
        }else{
         A_inj = meilleur_lars(lars = lars_inj, X = Xtilde, 
                              Y = Ytilde, mode = criterion, intercept = F, 
                              K = K, groupe = groupe)$A
        }
      }else if (select=="elasticnet") {
        lars_inj = renet(x = Xtilde, y = Ytilde, intercept = F, 
                        lambda = lambda)
        names(lars_inj)[4] = "coefficients"
        A_inj = meilleur_lars(lars = lars_inj, X = Xtilde, 
                              Y = Ytilde, mode = criterion, intercept = F, 
                              K = K, groupe = groupe)$A
      }else if(select=="adalasso"){
         resada=parcor::adalasso(X=Xtilde,y=Ytilde,k=K)    
         A_inj=c(resada$coefficients.adalasso)
         if(length(A_inj[A_inj!=0])>0){
            Xloc=Xtilde[,A_inj!=0]
            A_inj[A_inj!=0]=c(OLS(X=Xloc,Y=Ytilde,intercept=F)$beta)
         }
      }else if(select=="clere"){
         A_inj=A_clere(y=as.numeric(Ytilde),x=Xtilde,g=g)
         A_inj=A_inj[-1]#vraiment pas propre
      }else if(select=="spikeslab"){
         respike=spikeslab::spikeslab(x=X[,I1],y=Ytilde)
         A_inj=rep(0,times=ncol(Xtilde))
         A_inj[which(respike$gnet.scale!=0)]=OLS(X=Xtilde[,respike$gnet.scale!=0],Y=as.numeric(Ytilde),intercept=intercept)$beta 
      }else{#ridge
         if(ncol(Xtilde)>1){
            res_ridge=glmnet(as.matrix(Xtilde),Y,alpha=0,intercept=FALSE)
            A_inj=c(matrix(coef(res_ridge,lambda)))
            #ridge_pred = linearRidge(Ytilde~0+.,data=data.frame(Xtilde))
            #A_inj=coef(ridge_pred)
         }else{
            ridge_pred = OLS(X = Xtilde,Y=Ytilde,intercept=FALSE)#ridge has no sens on only one covariate
            A_inj=ridge_pred$beta
         }
      }
      A_pred = rep(0, times = ncol(X) + intercept)
      A_pred[I2 + intercept] = A_inj
      if(returning){
        Ytildebis=Y-as.matrix(X[,I2])%*%A_pred[I2 + intercept]
        Ytildebis=as.matrix(Ytildebis)
        if (select != "elasticnet" & select != "ridge" & select != "adalasso"  & select!="clere" & select!="spikeslab") {
          lars_retour=lars(x = as.matrix(X[,I1]), y = Ytildebis, type = select, 
                         intercept = intercept)
          A_retour = meilleur_lars(lars = lars_retour, X = as.matrix(X[,I1]), 
                                   Y = Ytildebis, mode = criterion, intercept = intercept, 
                                   K = K, groupe = groupe)$A
        }else if (select=="elasticnet"){
          lars_retour= renet(x = as.matrix(X[,I1]), y = Ytildebis, intercept =intercept, 
                  lambda = lambda)
          names(lars_retour)[4] = "coefficients"
          A_retour = meilleur_lars(lars = lars_retour, X = as.matrix(X[,I1]), 
                                   Y = Ytildebis, mode = criterion, intercept = intercept, 
                                   K = K, groupe = groupe)$A
        }else if(select=="adalasso"){
           resada=parcor::adalasso(X=X[,I1],y=Ytildebis,k=K)    
           if(intercept){
              A_retour=c(resada$intercept.adalasso,resada$coefficients.adalasso)
           }else{
              A_retour=c(resada$coefficients.adalasso)
           }
           Xloc=X[,I1][,resada$coefficients.adalasso!=0]
           A_retour[A_retour!=0]=c(OLS(X=Xloc,Y=Y,intercept=intercept)$beta)   
         }else if(select=="clere"){
            A_retour=A_clere(y=as.numeric(Y),x=X[,I1],g=g)
         }else if(select=="spikeslab"){
            respike=spikeslab::spikeslab(x=X[,I1],y=Y)
            res$compl$A=rep(0,times=ncol(X[,I1])+intercept)
            if(intercept){
               res$compl$A[c(1,1+which(respike$gnet.scale!=0))]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
            }else{
               res$compl$A[respike$gnet.scale!=0]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
            }
         }else{#ridge
          #ridge_pred = linearRidge(Ytildebis~.,data=data.frame(X[,I1]))
          #A_retour=coef(ridge_pred)
          res_ridge=glmnet(as.matrix(X[,I1]),Ytildebis,alpha=0,intercept=FALSE)
          A_retour=c(matrix(coef(res_ridge,lambda)))
        }
        
        if(intercept){
          A_pred[c(intercept, I1 + intercept)] =A_retour
        }else{
          A_pred[ I1 ] = A_retour
        }
      }else{
        if(intercept){
          A_pred[c(intercept, I1 + intercept)] = A_expl[c(intercept,I1 + intercept)] - B[c(1, I1 + 1), I2] %*% as.matrix(A_inj)
        }else{
          A_pred[ I1 ] = A_expl[I1] - B[c(1, I1 + 1), I2] %*% as.matrix(A_inj)
        }
        
      }    
      res$pred$A = A_pred
    #  res$pred$CVMSE = CVMSE(X = X[, which(A_pred[-1] != 0)], Y = Y, intercept = intercept, K = K, groupe = groupe)
      res$pred$BIC = BicTheta(X = X, Y = Y, intercept = intercept, beta = A_pred)
      res$pred$AIC=mon_AIC(theta=res$pred$A,Y=Y,X=X,intercept=intercept) 
      if(compl2){
         Xloc=X
         Xloc[,I2]=Xtilde
         res$compl2$A=OLS(X =as.matrix(Xloc) ,Y=Y,intercept = intercept)$beta
         res$compl2$BIC = BicTheta(X = as.matrix(Xloc), Y = Y, intercept = intercept, beta = res$compl2$A)
         res$compl2$CVMSE = CVMSE(X = Xloc, Y = Y, intercept = intercept, K = K, groupe = groupe)
      }
    }
  }else if (sum(Z) == 0 & (expl | pred)) {
    res$expl$A = res$compl$A
    if (pred) {
      res$pred$A = res$compl$A
    }
    if (explnew){
       res$expl2$A = res$compl$A
    }
  }
  return(res)
}