#' MCMC algorithm to find a structure between the covariates
#' @description This function computes a random walk based on a full generative model on the dataset. We optimize a BIC-like criterion to find a model of sub-regressions within the covariates. If marginal density are unknown, Gaussian Mixture models are used automatically. We obtain the best structure as an adjacency matrix (binary squared matrix) that corresponds to the Directed Acyclic Graph of dependencies within the covariates.
#' 
#' @details
#' At each step we compare several candidates that are the local structure modified at one place (one coefficient of the adjacency matrix). Knowing the local structure a candidate is then just defined by the index of the position we modify in this local structure. So strategies are just choices of list of indices.
#' 
#' To avoid local extrema we allow constraints relaxation. If a modification is not feasible (it generates cycles for example) then the candidate is not rejected but modified. In fact, if we want to modify \code{Z[i,j]} then we modify Z at each position that makes the modification of \code{Z[i,j]} not feasible. It is like allowing several steps in one, a kind of simulated anhealing but without parameter to tune.
#' 
#'@param X the matrix of the dataset containing p correlated covariates (with n individuals)
#'@param Z (optional) initial structure. Binary adjacency matrix of size p. if NULL zero matrix is used
#' @param Bic_null_vect p-sized vector of the BIC values of the null hypothesis (used for independent variables). If NULL then it would be computed based on Gausian Mitures hypothesis.
#'@param candidates strategy to define a neighourhood (list of candidates). Each new candidate is a modification of the current model as an adjacency matrix. So candidates are defined by the position that will be modified in Z. One modification for each candidate. We have then several neighbourhood to propose. 0:row and column (randomly chosen),-1:column only (randomly chosen), int>0:random int candidates at each step, -2 : all (but the diagonal) so \code{p^2-p} candidates, -3 : non-zeros (at each step we test all possible link removal). Each strategy gives a distinct number of candidates at each step.
#'@param reject 0: constraint relaxation (if a candidate is not feasible then we modify it to make it feasible by deleting not compatible links), 1: reject mode (if a candidate is not feasible, we don't look at it).
#' @param methode  parameter for OLS (matrix inversion in Ordinary Least Squares) 1:householderQr, 2:colPivHouseholderQr
#'@param p1max maximum complexity (number of explaining covariates) for a sub-regression (positive integer)
#'@param p2max maximum number of sub-regressions (positive integer)
#'@param Maxiter number of steps (positive integer)
#'@param plot (boolean) TRUE: returns for each step the type of move, complexity and BIC. If nbini>1 then it returns the values associated to the chain that found the best BIC.
#'@param best (boolean) TRUE: systematically jumps to the best BIC seen ever when seen (it is stored even if best=FALSE)
#'@param better (boolean) TRUE: systematically jumps to the best candidate if better than stationarity (random jump weighted by the BIC otherwise)
#'@param random (boolean) if FALSE:moves only to improve and only to the best. Otherwise random jump weighted by the BIC if no deterministic jump due to parameters \code{best} and/or \code{better}.
#'@param verbose level of printed informations during the walk. 0:none, 1:BIC,step and complexity when best BIC found 2:BIC, step, complexity, nb candidates and best candidate when best BIC found.
#'@param nb_opt_max stop criterion defining how many times the chain can walk (or stay) on the max found
#'@param exact (boolean) If exact sub-regression is found it gives its content (another verbose mode). One of the covariates can then be deleted manually by the user without loss of information.
#'@param nbini Number of initialisations (using initialisation based on correlation matrix if Z is NULL). if NULL and Z is NULL : only one chain starting with zero matrix (model without any sub-regression)
#'@param star (boolean) to compute BIC* instead of BIC (stronger penalization of the complexity based on a hierarchical uniform hypothesis on the probability of each structure). WARNING : star=TRUE implies p2max<=p/2.
#'@param clean (boolean) if TRUE then we add cleaning steps at the end of the walk (testing each remainging 1 for removal). So it is only few additional steps with \code{candidates=-3}
#'@param ... optional parameters to be passed (for initialization).
#'@return a list that contains:
#'\item{Z}{The local structure of the last step (adjacency matrix)}
#'\item{Z_opt}{The best structure seen during the walk in terms of the BIC-like criterion.}
#'\item{bic_opt}{Value of the global BIC-like criterion associated to \code{Z_opt}}
#'\item{step_opt}{The index of the step where \code{Z_opt} was found}
#'\item{Bic_null_vect}{p-sized vector of the BIC values associated to the model without sub-regressions. For use in a later search.}
#'\item{bic_step}{if \code{plot=TRUE}, vector of the BIC at each step}
#'\item{complexity_step}{if \code{plot=TRUE}, vector of the complexities at each step (\code{=sum(Z)})}
#'\item{step}{if \code{plot=TRUE}, vector of the type of modification at each step.  0:delete, 1: add, 2: stationarity}
#'
#'
#' @examples
#'\dontrun{
#'   rm(list=ls())#clean the workspace
#'   
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
#'    
#'    #density estimation for the MCMC (with Gaussian Mixtures)
#'    density=density_estimation(X=X_appr,nbclustmax=10,detailed=TRUE)
#'    Bic_null_vect=density$BIC_vect# vector of the BIC found (1 value per covariate)
#'    
#'    #MCMC to find the structure
#'    res=structureFinder(X=X_appr,verbose=0,reject=0,Maxiter=900,plot=TRUE,
#'                nbini=20,candidates=-1,Bic_null_vect=Bic_null_vect,star=TRUE,p1max=15,clean=TRUE)
#'    hatZ=res$Z_opt #found structure (adjacency matrix)
#'    hatBic=res$bic_opt #associated BIC
#'    
#'    #looking inside the walk
#'   
#'par(mar=c(5,4,4,5)+.1)
#'plot(res$bic_step,type="l",col="red",ylab="BIC",
#'     sub="blue: complexity, red: BIC", main="Evolution of BIC and complexity during the walk")
#'par(new=TRUE)
#'plot(res$complexity_step,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
#'axis(4)
#'mtext("Complexity",side=4,line=3)
#'#legend("topleft",col=c("red","blue"),lty=1,legend=c("BIC","Complexity"))
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
#'   } 
#'
#'
#'@export
structureFinder<-function(X=X,Z=NULL,Bic_null_vect=NULL,candidates=-1,reject=0,methode=1,p1max=5,p2max=NULL,Maxiter=1,plot=FALSE,best=TRUE,better=FALSE,random=TRUE,verbose=1,nb_opt_max=NULL,exact=TRUE,nbini=NULL,star=TRUE,clean=TRUE,...){
  params=match.call()
  Wini=FALSE
  X=1*as.matrix(X)
  if(is.null(p2max)){
    p2max=ncol(X)+1 
  }
  if(star){
     p2max=floor(min(p2max,ncol(X)*0.644))
     p1max=floor(min(p1max,ncol(X)*0.644))
  }
  if(is.null(nb_opt_max)){
    nb_opt_max=Maxiter
  }
  if(is.null(Bic_null_vect)){
     Bic_null_vect=density_estimation(X=X,nbclustmax=10,verbose=FALSE,detailed=FALSE)$BIC_vect
  }
  if(!is.null(nbini)){
     if (nbini<1){nbini=NULL}
  }
  if(is.null(nbini)){
     if(is.null(Z)){
        Z=matrix(0,ncol=ncol(X),nrow=ncol(X))
     }
     if(reject==0){#relax mode
        res=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        return(res)
     }else{# reject mode
        res=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        return(res)
     }
  }else{
     BICnull=sum(Bic_null_vect)
     if(!("W" %in% names(params))){
        W=cor(X)
     }
     if(is.null(Z)){
        Z=matrix(0,ncol=ncol(X),nrow=ncol(X))
        Wini=TRUE
     }
     res=list()
     if(nbini>1){#first try with provided Z matrix (or null if not provided)
        if(reject==0){#relax mode
           res=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }else{# reject mode
           res=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }
        nbini=nbini-1#to finally have the exact number of tries
        if(clean){
           resclean=cleanZ(X=X,Z=res$Z_opt,Bic_null_vect=Bic_null_vect,star=star,verbose=verbose)#nettoyage colonnes puis ponctuel
           res$Z_opt=resclean$Z_opt
           res$bic_opt=resclean$bic_opt
        }
     }
    
     for(i in 1:nbini){
        if(Wini){#only if no Z provided
           Z=Winitial(W=W,X=X,p1max=p1max,Bic_null_vect=Bic_null_vect,p2max=p2max)
        }
        if(reject==0){#relax mode
           resloc=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }else{# reject mode
           resloc=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }
        if(length(res)==0){res=resloc}
        if(resloc$bic_opt<=min(res$bic_opt,BICnull)){
           res=resloc
        }
     }
     if(clean){
        resclean=cleanZ(X=X,Z=res$Z_opt,Bic_null_vect=Bic_null_vect,star=star,verbose=verbose)#nettoyage colonnes puis ponctuel
        res$Z_opt=resclean$Z_opt
        res$bic_opt=resclean$bic_opt
     }
     res$Bic_null_vect=Bic_null_vect
     return(res)
  }

}