#' compute the matrix R* of generalized correlation coefficients.
#' 
#' function checks for missing data separately for each pair using
#' \code{kern} function to kernel regress x on y, and conversely y on x. It
#'  needs the library `np' which reports R-squares of each regression. This function
#' reports their square roots with the sign of the Pearson correlation coefficients.
#' Its appeal is that it is asymmetric yielding causal direction information,
#' by relaxing the assumption of linearity implicit in usual correlation 
#' coefficients.
#' 
#' @param mym {A matrix of data on variables in columns}
#' @param nam {Column names of the variables in the data matrix}
#' @return  non-symmetric R* matrix of generalized correlation coefficients
#' @importFrom stats cor
#' @note The user may want to modify this function without touching \code{gmcmtx0}.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
## @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Vinod, H. D."Generalized Correlation and Kernel Causality with 
#'  Applications in Development Economics" in Communications in 
#'  Statistics -Simulation and Computation, 2015, 
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048} 
#' @keywords kernel regression, asymmetric
#' @examples
#' 
#' \dontrun{
#' set.seed(34);x=matrix(sample(1:600)[1:99],ncol=3)
#' colnames(x)=c("V1", "v2", "V3")
#' gmcmtxZ(x)
#' }
#' 
#' @export

gmcmtxZ <-
function(mym, nam=colnames(mym)){ 
# mym is a data matrix with n rows and p columns
# some NAs may be present in the matrix
p=NCOL(mym)
#print(c("p=",p))
out1=matrix(1,p,p)# out1 stores asymmetric correlations
for (i in 1:p){
x=mym[,i]
ava.x=which(!is.na(x))
for (j in 1:p){
if (j>i){ y=mym[,j]
ava.y=which(!is.na(y))#ava means non-missing
ava.both=intersect(ava.x,ava.y)
newx=x[ava.both]
newy=y[ava.both]
c1=cor(newx,newy)
sig=sign(c1)
mod.1=kern(dep.y=x, reg.x=y)
mod.2=kern(dep.y=y, reg.x=x)
corxy= sqrt(mod.1$R2)*sig
coryx= sqrt(mod.2$R2)*sig
out1[i,j]=corxy
out1[j,i]=coryx
}#endif
}#end j loop
}#end i loop
colnames(out1)=nam
rownames(out1)=nam 
return(out1)     
}
