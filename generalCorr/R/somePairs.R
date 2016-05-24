#' Function reporting kernel causality results as a 7-column matrix.
#'
#' Three criteria are used: 
#' 
#' (typ=1) reports ("Y", "X", "Cause",
#' "SD1apd", "SD2apd", "SD3apd", "SD4apd") nameing variables identifying "cause"
#' and measures of stochastic dominance using absolute values of kernel
#' regression gradients comparing regresson of X on Y with that of Y on X.
#' 
#' {(typ=2)} 
#'  reports ("Y", "X", "Cause", "SD1res", "SD2res", "SD3res", "SD4res")
#'  and measures of stochastic dominance using absolute values of kernel
#'  regression residuals comparing regresson of X on Y with that of Y on X.
#' 
#' {(typ=3)} 
#' reports ("Y", "X", "Cause", "r*x|y", "r*y|x", "r", "p-val")
#' containing generalized correlation coefficients r*, 'r' refers to.
#' Pearson correlation coefficient p-val is the p-value for 
#' testing the significance of Pearson 'r'
#' @param mtx {the data matrix in the first column is paired with all others.}
#' @param dig {=6 (default) number of digits for reporting.}
#' @param verbo {=FALSE (default). Make it TRUE for printing steps.}
#' @param typ {must be 1 (default), 2 or 3 for the three criteria.}
#' @param rnam {=FALSE (default). Make it TRUE if created rownames are wanted}
#' @return A matrix containing causal identification results.  The last matrix
#' outVote reports summary results for all three criteria in column `sum'.
#' @keywords  causal criteria, SD1, SD2, SD3, SD4, generalized correlations
#'
#' @note will not work without several functions and packages in memory
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{some0Pairs}}.
#' @references H. D. Vinod "Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics" in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#'
#'
#' @examples
#'
#' \dontrun{
#' data(mtcars)
#' somePairs(mtcars)
#' }
#' 
#' @export

somePairs <-
function(mtx,dig=6,verbo=FALSE,typ=1,rnam=FALSE){

n=NROW(mtx)
p=NCOL(mtx)
rna=rep(NA,p)#for storing row names
print(c("n,p,digits",n,p,dig))
if (typ==1) print("absolute apds compared")
if (typ==2) print("absolute residuals compared")
if (typ==3) print("r* compared")
nam=colnames(mtx)
npair=p-1
outcause=matrix(NA, nrow=npair,ncol=7)
ii=0
#following loop is such that i<=j, which means
#[i,j] will have sup-diagonal
for (i in 2:p){
x0=mtx[,i]   #i has x   or all other columns
y0=mtx[,1]   #first col. has y  NOT x
na2=napair(x0,y0)
x=na2$newx
y=na2$newy
print(c("i,non-missing",i,length(x)))
ii=ii+1
print(c("i,ii",i,ii))
rna[i]=paste("1",i,sep=".")
if(typ==1)arxy=abs_stdapd(x,y)
if(typ==1)aryx=abs_stdapd(y,x)
if(typ==2)arxy=abs_stdres(x,y)
if(typ==2)aryx=abs_stdres(y,x)
if(typ<3){
crit4=comp_portfo2(arxy,aryx)
if(verbo) print(crit4)
round.crit4=round(crit4,dig)

outcause[ii,4:7]= round.crit4
outcause[ii,1]= nam[1]
outcause[ii,2]= nam[i]
outcause[ii,3]= nam[i]#i has x and x is the implicit cause
if (crit4[1]>0) outcause[ii,3]= nam[1]  #SD1>0 then cause=y
} #endif typ<3
if(typ==3){
rst=rstar(x,y)
rxy=rst$corxy
ryx=rst$coryx
#print(c(rxy,ryx,ii,i))
del=rxy^2-ryx^2
#del>0 means rxy>ryx or x on y good or cause=y
outcause[ii,4]=round(rst$corxy,dig)
outcause[ii,5]=round(rst$coryx,dig)
outcause[ii,6]=round(rst$pearson.r,dig)
outcause[ii,7]=round(rst$pv,dig)
outcause[ii,1]= nam[1]
outcause[ii,2]= nam[i]
outcause[ii,3]= nam[i]#cause= x
if (del>0) outcause[ii,3]= nam[1]#cause=y
} #endif typ==3
} #end of i loop
namout=c("Y", "X","Cause","SD1","SD2","SD3","SD4")
colnames(outcause)=namout
if(typ==1)  namout=c("Y", "X","Cause","SD1apd","SD2apd","SD3apd","SD4apd")
if(typ==2)  namout=c("Y", "X","Cause","SD1res","SD2res","SD3res","SD4res")
if(typ==3)  namout=c("Y", "X","Cause","r*x|y","r*y|x","r","p-val")
#rownames(outcause)=rownames(out2)
colnames(outcause)=namout
if(rnam) rownames(outcause)=rna[2:p] #first row name slot=NA
if(verbo)print(outcause)
if(verbo)print(xtable(outcause))
return(outcause)
}
