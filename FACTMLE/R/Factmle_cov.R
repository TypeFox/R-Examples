#' Calculates the Maximum likelihood Factor analysis with a covariance Matrix.
#' 
#' @param S         The Covariance Matrix. It is a p*p numeric matrix, where p is the number of variables.
#' @param rnk       Rank constraint for the Factor analysis problem. It must a positive integer less than the number of variables p
#' @param Psi_init  The initial value of Psi. It is a p*1 numeric vetor, where p is the number of variables.Default value is a vector of uniform random numbers.
#' @param lb        The lower bound on the Psi values. The default value is set to 0.05
#' @param index     This option is for modified version of factmle.The default value is a null vector. If assigned a zero vector, it will perform MLFA keeping some of the Psi 
#' values specified by the index at a specifed level *lb2*
#' @param lb2       This option of modified version of factmle algorithm. The default value is 0.001. The Psi values 
#' specified by the *index* is kept constant at *lb2* while doing MLFA.
#' @param tol       Precision parameter. Default is 10^-7
#' @param Max_iter  Maximum number of iterations. Default is 1000.
#' @return A list with the following components
#' \describe{
#' \item{Psi}{A vector containing the unique variances.}
#' \item{Lambda}{A p*rnk matrix containing the factor loadings in the columns.}
#' \item{Nll}{A vector containing the negative Log-likelihood values at every iteration.}
#' \item{Nllopt}{The value of the negative log-likelihood upon convergence.}
#' }
#' @seealso \code{eigs_sym}
#' @export
#' @importFrom rARPACK eigs_sym
#' @importFrom stats runif
#' @importFrom stats cov
#' @examples
#' 
#' 
#' library(MASS)
#' library(stats)
#' Psi=runif(15,min=0.2,max=1.3)
#' Lambda=mvrnorm(n=15,mu=rep(0,3),Sigma = diag(rep(1,3)))
#' data=mvrnorm(n=5000,mu=rep(0,15),Sigma = diag(Psi)+Lambda%*%t(Lambda))
#' S=cov(data)
#' x=Factmle_cov(S,3)
#' 
#' 
Factmle_cov<-function(S,rnk,Psi_init=c(),lb=.01,index=c(),lb2=.01,tol=10^-7,Max_iter=1000)
{
  l=length(Psi_init);
  dim=ncol(S);
    if (l==0)
  {
    Psi_init=runif(dim,min=lb,max=1);
  }
  
  f=rep(1,Max_iter);
  Psi=Psi_init;
  k=2;
  dim=length(Psi);
  check=1;
  while (check)
  {
    Old_Psi=Psi;
    x=1/Psi;
    x_half=x^(.5);
    y=rep(x_half,dim);
    x1=matrix(y,nrow=dim,ncol=dim,byrow = FALSE)
    x2=matrix(y,nrow = dim,ncol = dim,byrow = TRUE)
    s1=(x1*S)*x2;
    s1=0.5*(s1+t(s1))
    v1=eigs_sym(s1,rnk,which = "LM",opts = list(retvec = TRUE));
    v=v1$vectors;
    d=v1$values;
    diags=diag(S);
    f[k]=calc_fval(diags,x,d)
    diff_d=pmax(0,(1-1/d))
    
    x1=matrix(rep(x_half,dim),nrow = dim,ncol = dim,byrow = FALSE)
    x2=matrix(rep(1/x_half,dim),nrow = dim,ncol = dim,byrow = TRUE)
    A=(x1*S)*x2;
    A=t(v)%*%A;
    
    y=matrix(rep(diff_d,dim),nrow=dim,ncol = rnk,byrow = TRUE)
    B=v*y;
    diff_psi_0 = colSums(t(B)*A)
    Psi=pmax(diags-diff_psi_0,lb);
    Psi[index]=lb2;
    
    c=abs((f[k]-f[k-1])/f[k-1]);
    if ( (c < tol)||(k>Max_iter)  )
    {
      check=0;
    }
    k=k+1;
  }
  Lambda=calculate_Lambda(Psi,S,rnk)
  
  out=list(Psi=Psi,Lambda=Lambda,Nllopt=f[k-1],Nll=f[c(2:(k-1))])
  return(out)
}

calculate_Lambda<-function(Psi,S,rnk)
{
  
  dim=length(Psi);
  Psi_minus_half=1/(Psi^(.5) )
  
  y1=matrix(rep(Psi_minus_half,dim),nrow = dim,ncol =dim,byrow=FALSE)
  y2=matrix(rep(Psi_minus_half,dim),nrow=dim,ncol=dim,byrow=TRUE)
  
  Sstar=(y1*S)*y2;
  Sstar=0.5*(Sstar+t(Sstar));
  
  V1=eigs_sym(Sstar,rnk,which="LM",retvec=TRUE)
  V=V1$vectors;
  d=V1$values;
  modulus=sqrt(pmax(d-1,0));
  
  x1=matrix(rep(sqrt(Psi),rnk),nrow = dim,ncol = rnk,byrow = FALSE)
  x2=matrix(rep(modulus,dim),nrow=dim,ncol=rnk,byrow = TRUE)
  
  Lambda=(x1*V)*x2;
  
}

calc_fval <- function(diags,x,eig_val)
{
  fval= sum(-log(x) ) + t(diags)%*%x +sum(log(pmax(1,eig_val)) -pmax(1,eig_val) +1);
  return(fval)
}
