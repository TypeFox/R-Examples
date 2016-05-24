#' Function reporting kernel causality results as a 7-column matrix where
#' stochastic dominances orders are weighted by c(1.2,1.1, 1.05, 1) to
#' compute an overall result for all orders of stochastic dominance.
#'
#'
#' Three criteria are used and further summarized into one matrix called \code{outVote}: 
#'   
#' typ=1 reports ("Y", "X", "Cause",
#' "SD1apd", "SD2apd", "SD3apd", "SD4apd") nameing variables identifying "cause"
#' and measures of stochastic dominance using absolute values of kernel
#' regression gradients comparing regresson of X on Y with that of Y on X.
#' 
#' 
#' typ=2 reports ("Y", "X", "Cause", "SD1res", "SD2res", "SD3res", "SD4res")
#' and measures of stochastic dominance using absolute values of kernel
#' regression residuals comparing regresson of X on Y with that of Y on X.
#' 
#' 
#' typ=3 reports ("Y", "X", "Cause", "r*x|y", "r*y|x", "r", "p-val")
#' containing generalized correlation coefficients r*, 'r' refers to.
#' Pearson correlation coefficient p-val is the p-value for 
#' testing the significance of 'r'
#' @param mtx {the data matrix in the first column is paired with all others.}
#' @param dig {=6 (default) number of digits for reporting.}
#' @param verbo {=FALSE (default). Make it TRUE for printing steps.}
#' @param rnam {=FALSE (default). Make it TRUE if created rownames are wanted.}
#' @param wt {allows user to choose alternative weights for SD1 to SD4.}
#' @param sumwt {=4(default) sum of weights can be changed here.}
#' @return Prints three matrices detailing results for Cr1, Cr2 and Cr3.
#' It returns a matrix called `outVote' which summarizes the three criteria.
#' In general, a positive sign for weighted sum reported in the column `sum' means
#' that the first variable listed as the input to this function is the `kernel cause.'  
#' For example, crime `kernel causes' officer deployment (not vice versa) is indicated by 
#' the positive sign of `sum' (=3.175) reported for that example included in this package.
#' @importFrom xtable xtable
#' @note This will not work without several functions and packages in R memory.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{somePairs}}
#' @references "Generalized Correlation and Kernel Causality with
#'    Applications in Development Economics" in Communications in
#'    Statistics -Simulation and Computation, 2015,
#'    \url{http://dx.doi.org/10.1080/03610918.2015.1122048}
#' @keywords causal criteria, SD1, SD2, SD3, SD4, generalized correlations
#' @note The output matrix last column for `mtcars' example
#' has the sum of the scores by the three criteria
#' combined. If `sum' is positive, then variable X (mpg) is more likely to have been
#' engineerd to kernel cause the response variable Y, rather than vice versa.
#' @note The European Crime data has all three criteria correctly suggesting that
#' high crime rate kernel causes the deployment of a large number of police officers.
#' 
#' @examples
#'
#'
#' \dontrun{
#' data(mtcars) #brings the Motor Trend data into the memory of R
#' some0Pairs(mtcars) # first variable is mpg and effect on mpg is of interest
#' }
#' 
#' \dontrun{
#' data(EuroCrime)
#' attach(EuroCrime)
#' some0Pairs(cbind(crim,off))
#' }
#' 
#' 
#' @export

some0Pairs <-
function(mtx,dig=6,verbo=TRUE,
rnam=FALSE,wt=c(1.2,1.1,1.05,1),sumwt=4){
n=NROW(mtx)
p=NCOL(mtx)
npair=p-1  #number of pairs considered
cr1=rep(NA,npair)   #need to initialize to avoid R errors
cr2=rep(NA,npair)
cr3=rep(NA,npair)
crall=rep(NA,npair)

rna=rep(NA,p)#for storing row names
outVote=matrix(NA,nrow=npair,ncol=7)
colnames(outVote)=c("X", "Y","Cause","Cr1","Cr2","Cr3","sum")
#print(c("n,p",n,p,"digits=",dig))
nam=colnames(mtx) #R makes nam=NULL of lenghth 0 if mtx column names Missing
if (length(nam)==0) nam=paste("V",1:p,sep="")
for (typ in 1:3){
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
#if (i==2 & typ==1)   print(summary(cbind(x,y)))
#if (i==2)print(c("i=",i,"non-missing y=",length(y)))
if(verbo){
if (i>2)print(c("i=",i,"non-missing y=",length(y)))
} #end verbo if
if (length(x)<5) {print("available observations<5")
break}
ii=ii+1
if(verbo)print(c("i=",i,"ii=",ii))
rna[i]=paste("1",i,sep=".")

if (typ==1) outVote[ii,1]= nam[1]
if (typ==1) outVote[ii,2]= nam[i] #second col has i-th name
if(typ==1)arxy=abs_stdapd(x,y) #compare apd's
if(typ==1)aryx=abs_stdapd(y,x)
if(typ==2)arxy=abs_stdres(x,y) #compare residuals
if(typ==2)aryx=abs_stdres(y,x)
if(typ<3){
crit4=comp_portfo2(arxy,aryx)
if(verbo) print(crit4)
round.crit4=round(crit4,dig)
wtdsign= wt*sign(round.crit4)
av.crit4=round((sum(wtdsign,na.rm=TRUE)/sumwt),dig)
outcause[ii,4:7]= round.crit4
outcause[ii,1]= nam[1]
outcause[ii,2]= nam[i]
outcause[ii,3]= nam[i]#i has x and x is the implicit cause
if (av.crit4>0) outcause[ii,3]= nam[1]  #SD1>0 then cause=y
if (typ==1) {outVote[ii,4]= av.crit4
cr1[ii]=av.crit4}
if (typ==2) {outVote[ii,5]= av.crit4
cr2[ii]=av.crit4}
} #endif typ<3
if(typ==3){
rst=rstar(x,y)
rxy=rst$corxy
ryx=rst$coryx
#print(c(rxy,ryx,ii,i))
del=rxy^2-ryx^2
#del>0 means rxy>ryx or x on y good or cause=y
cr3[ii]=as.numeric(sign(del))
outcause[ii,4]=round(rst$corxy,dig)
outcause[ii,5]=round(rst$coryx,dig)
outcause[ii,6]=round(rst$pearson.r,dig)
outcause[ii,7]=round(rst$pv,dig)
outcause[ii,1]= nam[1]
outcause[ii,2]= nam[i]
outcause[ii,3]= nam[i]#cause= x
if (del>0) outcause[ii,3]= nam[1]#cause=y

outVote[ii,6]= sign(del)
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
#return(outcause)
if(typ==1) outCr1=outcause
if(typ==2) outCr2=outcause
if(typ==3) outCr3=outcause
}  #end for typ loop
#print(nam)
for (j in 1:npair){
cr13=c(cr1[j], cr2[j],cr3[j])
crall[j]= round(sum(cr13,na.rm=TRUE),dig)
if(verbo)print(c("cr=",cr1[j],cr2[j],cr3[j],crall[j]))
outVote[j,7]= crall[j]
if(!is.na(crall[j])) {
if(crall[j]>0)outVote[j,3]= nam[1]
if(crall[j]<=0)outVote[j,3]= nam[j+1] }
} #end j loop
#print("summary Cr1 to Cr3: positive values imply X is the kernel cause")
if(verbo)print(outVote)
if(verbo)print(xtable(outVote))
list(outCr1= outCr1,
outCr2= outCr2, outCr3= outCr3,
outVote=outVote)
}
