#' CorReg: see www.correg.org for article and Phd Thesis about CorReg.
#' @name CorReg-package
#' @aliases CorReg-package
#' @docType package
#' @title Quick tutorial for CorReg package
#' @description Sequential linear regression based on a structural equation model(explicit correlations). 
#' It permits to face highly correlated datasets.
#' We first search for an explicit model of correlations within the covariates by linear regression, 
#' then this structure is interpreted and used to reduce dimension and correlations for the main regression on the response variable.

#' @author Maintainer: Clement THERY <clement.thery@@arcelormittal.com>

#' @references Model-based covariable decorrelation in linear regression (CorReg): application to missing data and to steel industry. C Thery - 2015. see \url{http://www.theses.fr/2015LIL10060} to read the associated PhD Thesis.
#' @keywords package


#' @examples
#'    \dontrun{
#' require(CorReg)
#'    #dataset generation
#'    base=mixture_generator(n=15,p=10,ratio=0.4,tp1=1,tp2=1,tp3=1,positive=0.5,
#'                           R2Y=0.8,R2=0.9,scale=TRUE,max_compl=3,lambda=1)

#'    X_appr=base$X_appr #learning sample
#'    Y_appr=base$Y_appr #response variable for the learning sample
#'    Y_test=base$Y_test #responsee variable for the validation sample
#'    X_test=base$X_test #validation sample
#'    
#'    TrueZ=base$Z#True generative structure (binary adjacency matrix)
#'    #Z_i,j=1 means that Xj linearly depends on Xi
#'    
#'    #density estimation for the MCMC (with Gaussian Mixtures)
#'    density=density_estimation(X=X_appr,nbclustmax=10,detailed=TRUE)
#'    Bic_null_vect=density$BIC_vect# vector of the BIC found (1 value per covariate)
#'    
#'    #MCMC to find the structure
#'    res=structureFinder(X=X_appr,verbose=0,reject=0,Maxiter=900,
#'                nbini=20,candidates=-1,Bic_null_vect=Bic_null_vect,star=TRUE,p1max=15,clean=TRUE)
#'    hatZ=res$Z_opt #found structure (adjacency matrix)
#'    hatBic=res$bic_opt #associated BIC
#'    
#'    #BIC comparison between true and found structure
#'    bicopt_vect=BicZ(X=X_appr,Z=hatZ,Bic_null_vect=Bic_null_vect)
#'    bicopt_vrai=BicZ(X=X_appr,Z=TrueZ,Bic_null_vect=Bic_null_vect)
#'    sum(bicopt_vect);sum(bicopt_vrai)
#'    
#'    #Structure comparison
#'    compZ=compare_struct(trueZ=TrueZ,Zalgo=hatZ)#qualitative comparison
#'    
#'    #interpretation of found and true structure ordered by increasing R2
#'    readZ(Z=hatZ,crit="R2",X=X_appr,output="all",order=1)# <NA>line : name of subregressed covariate
#'    readZ(Z=TrueZ,crit="R2",X=X_appr,output="all",order=1)# <NA>line : name of subregressed covariate
#'    
#'    #Regression coefficients estimation
#'     select="NULL"#without variable selection (otherwise, choose "lar" for example)
#'    resY=correg(X=X_appr,Y=Y_appr,Z=hatZ,compl=TRUE,expl=TRUE,pred=TRUE,
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

#'   barplot(as.matrix(MSE),main="MSE on validation dataset", sub=paste("select=",select))
#'   abline(h=MSE_complete,col="red")
#'    }
NULL
#la puissance CorReg!!!!!!!!!