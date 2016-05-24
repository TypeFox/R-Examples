#' Compute generalized (ridge-adjusted) partial correlation coefficients from matrix R*.
#'
#' This function uses a generalized correlation matrix R* as input to compute
#' generalized partial correlations between \eqn{X_i} and one of the remaining
#' variables after removing the effect of all other variables in the matrix.
#'
#' @param gmc0 has a p by p matrix R* of generalized corr coefficients.
#' @param dig has the number of digits for reporting (=4, default)
#' @param idep is the column number of the first variable (=1, default)
#' @param verbo Make this TRUE for detailed printing of computational steps
#' @param incr {incremental constant for iteratively adjusting ridgek
#'        where ridgek is the constant times the identity matrix used to
#'    make sure that the gmc0 matrix is positive definite and iteratively
#'  increased till all partial correlations are within the [-1,1] interval.}
#' @return A five column `out' matrix containing partials. The first column=
#'   has the name of i=idep variable. The
#'    second column has the name of the j variable, while the third column has  r*(i,j|k).
#'   The 4-th column has  r*(j,i|k) (denoted partji), 5-th column has rijMrji
#'   that is difference in absolute values (abs(partij) - abs(partji)).
#'
#' @note The ridgek constant created by the function may not be large enough to make sure that
#'  that other pairs of r*(i,j|k) are within the [-1,1] interval. The user may have to choose
#'  a suitable input `incr' to get the partials in  the [-1,1] interval. This function
#'  needs cofactor, minor and other functions in the memory.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{parcor_ijk}}.
#' @keywords partial correlations
#' @references Vinod, H. D. "Generalized Correlations and Instantaneous
#'  Causality for Data Pairs Benchmark," (March 8, 2015)
#'  \url{http://ssrn.com/abstract=2574891}
#'
#' @references Vinod, H. D. "Matrix Algebra Topics in Statistics and Economics
#'  Using R", Chapter 4 in Handbook of Statistics: Computational Statistics
#'  with R, Vol.32, co-editors: M. B. Rao and C.R. Rao. New York:
#'  North Holland, Elsevier Science Publishers, 2014, pp. 143-176.
#' @examples
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c("V1", "v2", "V3")
#' gm1=gmcmtx0(x)
#' parcor_ridg(gm1, idep=1)
#' }
#'
#' @export

parcor_ridg <-
function(gmc0,dig=4,idep=1, verbo=FALSE,incr=3){
#input is asymmetric generalized correlation matrix
#dig=4 #digits of accuracy desired
#We want r* of idep, with respect to others removing all others as k
n=NROW(gmc0)
p=NCOL(gmc0)
if (n!=p) stop("gmc0 is not square matrix")
nam=colnames(gmc0) #R makes nam=NULL of lenghth 0 if gmc0 column names Missing
if (length(nam)==0) nam=paste("V",1:p,sep="")
print(c("We want partial Corr of",nam[idep],"w.r.t. others"))
j.other=setdiff(1:n, idep)
e0=eigen(gmc0)
sort.e0=sort(e0$values)
sort.abse0=sort(abs(e0$values))
min.e0=sort.e0[1];min.e0
diff.e0=incr*(sort.abse0[2]-sort.abse0[1])
print(c("increment=",diff.e0))
#diff.e0=round(diff.e0,3)
if (diff.e0==0) {print("difference between cosecutive eigenvalues is zero")
diff.e0=incr*(sort.abse0[3]-sort.abse0[1])}
if (diff.e0==0) {print("difference between cosecutive eigenvalues still zero")
diff.e0=incr*(sort.abse0[4]-sort.abse0[1])}

sgn.e0=sign(Re(min.e0));sgn.e0
ridgek=0 #initialize for positive definite case
if (sgn.e0==-1)ridgek=sgn.e0*round(Re(min.e0)+diff.e0,3)
print(c("ridgek=",ridgek))
#now we ridge adjust the generalized correlation matrix gmc0
#by adding ridgek times identity matrix

for (jridge in 1:10){
print(c("ridgek=",ridgek))

if(jridge>=2)ridgek=ridgek+incr*diff.e0
gmc1=gmc0+ridgek*diag(ncol(gmc0))
nam=colnames(gmc0)
if(jridge==1){
print(c("We want partial Corr of",nam[idep],"w.r.t. others"))
}#eindif printing
n=length(nam)
print(c("nam length=",n))
dig=4 #digits of accuracy desired
out1=matrix(NA, nrow=n-1, ncol=4)
partij=rep(NA,n-1) #place holders
partji=rep(NA,n-1)
ii=0
for(i in 2:n){
if (i==2){
print(c("partial r* between",nam[1],nam[i]))}
if (i==n){
print(c("partial r* between",nam[1],nam[i]))}
p1=parcor_ijk(gmc1,1,i)
ii=ii+1
partij[ii]=p1$ouij
partji[ii]=p1$ouji
}#end i lop
if (verbo){
print("partij")
print(partij)
print("partji")
print(partji)
} #endif verbo
#if ( (max(abs(partij))>1) | (max(abs(partji))>1)){
#print(c("iteration for ridgek=",ridgek, jridge))}

if ((max(abs(partij))<1) & (max(abs(partji))<1)) break
}#end jridge loop
rijMrji=(abs(partij)-abs(partji))
print(c("final ridgek=",ridgek))
cb1=cbind(partij,partji,rijMrji)
cb2=apply(cb1,2,round,dig)
if(verbo) print(cb2)
m=length(partij)
nam2=rep(nam[idep],m)
nam3=nam[j.other]
out=cbind(nam2,nam3,cb2)
if(verbo) print(out)
return(out)
}
