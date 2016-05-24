# ' Selection method based on p-values (coefficients)
# ' @description This function computes linear regression using the lm function and then delete covariates 
# ' whose associated p-value is too small (compared to a user defined parameter \code{pvalmin}) and iterates still it finds a model with only high p-values (or if it finds nothing)
# ' @param Y the response variable
# ' @param X the dataset of the covariates (without the response)
# ' @param pvalmin the maximal bound for p-value associated to remaining coefficients
# ' @param bonferroni boolean defining wether a Bonferroni correction is applied or not
# ' @param A optional vector of coefficients (length ncol(X)+1) where 0 values will coerce zero values in the final vector of regression coefficients
# ' @export
# ' 
# '@examples
# '    \dontrun{

# '    require(CorReg)
# ' data=mixture_generator(n=15,p=5,valid=0,ratio=0,Amax=3)#dataset generation
# ' X=as.matrix(data$X_appr)
# ' Y=as.matrix(data$Y_appr)
# ' TrueA=data$A
# ' res=cleanYtest(Y = Y, X = X, pvalmin = 0.5)
# '     }
# ' 
cleanYtest<-function (Y = Y, X = X, pvalmin = 0.05, bonferroni=F,A=NULL) 
{
  p=ncol(X)
  qui=NULL
  if(!is.null(A)){
    qui=which(A[-1]!=0)
    X=X[,qui]
  }
  #on regarde chaque coef donc on boucle jusqu'a stabilite
  pvalminini=pvalmin
  change=T#changement potentiel
  quinonzero=1:(ncol(X))
  loc=length(quinonzero)
  A=rep.int(0, times=ncol(X)+1)
  while(change & ncol(X)>0){
    if(bonferroni){pvalmin=pvalminini/(ncol(X))}
      lmloc=lm(Y~.,data=data.frame(X))
      summar=summary(lmloc)
      coefs_pval=coef(summar)[,4]#p-values des coefficients
      quivarzero=which(coefs_pval[-1]>pvalmin)
    if(length(quivarzero)>0){#on elague juste les coefs pourris
      quinonzero=quinonzero[-quivarzero]
      X=as.matrix(X[,quinonzero])
      loc=length(quinonzero)
    }else{#on n'a rien change
      change=F
    }
  }
  #on regarde la constante et on l'enleve si besoin
  if(coefs_pval[1]>pvalmin){
    Aloc=lm(Y~0+.,data=data.frame(X))$coefficients
    quinonzero=quinonzero[-1] 
    A[quinonzero+1]=Aloc
  }else{
    Aloc=lmloc$coefficients
  }
  A[c(1,quinonzero)]=Aloc
  if(!is.null(qui)){#on avait un A plus grand
    Along=rep(0,times=(p+1))
    Along[c(1,qui+1)]=A
    return(Along)
  }else{
    return(A)
  }
}

