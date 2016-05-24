granova.2w <- function(data, formula = NULL, fit = "linear", ident = FALSE,  offset = NULL, ...){

# offset = NULL,
# Function depicts two-way ANOVA using data-based contrasts.
# Argument 'data.A.B' should be an n X 3 dataframe.  If the rows are named uniquely, then points
# will be able to be identified with those labels, otherwise the row number is used.
# The first column should be the response values, the second and third columns should
# be factors defining the levels.  The factors are ordered
# by the data (row/col means); they are are not taken to be ordered at the outset.
# Function 'scatter3d' be must be available in R, which is easiest to do by loading its source, Rcmdr (thanks,John Fox).
# The 'fit' is defaulted to 'linear' where interaction is not fit (i.e. flat
# surface based on predicting cell means using row and column marginal effects).
# Replace 'linear' with, say, quadratic to produce a curved surface.
# Note: right click on the scatterplot to terminate 'Identify' and return the output from the function.
# ... sends other optional commands to scatter3.
# If offset is NULL, then identif3d default is used for offset


require(rgl)
require(tcltk)
require(mgcv)
#car:::scatter3d
#car:::Identify3d

data.A.B<-data
mtx <- is.data.frame(data.A.B)
if(!mtx)data.A.B <- data.frame(data.A.B)
if(!is.factor(data.A.B[,2]))data.A.B[,2] <- as.factor(data.A.B[,2])
if(!is.factor(data.A.B[,3]))data.A.B[,3] <- as.factor(data.A.B[,3])

A.name<-names(data.A.B)[2]
B.name<-names(data.A.B)[3]
resp.name<-names(data.A.B)[1]

vA <- length(unique(data.A.B[,2]))
vB <- length(unique(data.A.B[,3]))
N  <- dim(data.A.B)[1]
rnd2 <- function(x)round(x,2)

A <-  data.A.B[,2]
B <-  data.A.B[,3]
yy <- data.A.B[,1]

#Means in each cell
mns <- tapply(yy, list(A,B), mean)
mns.vec <- as.vector(mns)
mns.matrx<-matrix(mns.vec,ncol=vB)

#Next two lines are of class 'array'; must be changed to class numeric
cell.mnsA <- tapply(yy, list(A), mean)
cell.mnsB <- tapply(yy, list(B), mean)

#For scat3d; points to group means as separate group and thus colors them differently.
group.factor <- factor(c(rep(0,N), rep(1, vA*vB)))

#mnsA,B are lengths of A,B vectors with appropriate means
mnsA <- as.numeric(cell.mnsA)[A]
mnsB <- as.numeric(cell.mnsB)[B]
ordA <- order(cell.mnsA)
ordB <- order(cell.mnsB)

mns.matrx <- mns.matrx[ordA,ordB]
cell.cnts <- table(A, B, dnn = colnames(data.A.B[2:3]) )[ordA, ordB]

#Line 55 is test.
dimnames(mns.matrx)<-dimnames(cell.cnts)

#grandmean
grndmean <- mean(yy)

#these two vectors are effects cum data-based contrasts for A & B factors
facA.mn.cntrst <- mnsA - grndmean
facB.mn.cntrst <- mnsB - grndmean

#Allow user defined model for basic two way analysis of variance summary.
if(is.null(formula)){formula<-as.formula(print(paste(resp.name,"~",A.name,"*",B.name),quote=F))}
aov.yy <- aov(formula=formula,data=data.A.B)

#Trying to put means in: something of a hack.  Adding points at the means for each cell,
#then using the group feature of scatter3d to give them a different color.   Scatter3d
#calculates the separate regression plane for the means, but it is the same as the plane
#calculated for the data.  

mnsAA <- rep(as.numeric(cell.mnsA), length(unique(B))) - grndmean
mnsBB <- rep(as.numeric(cell.mnsB), ea = length(unique(A))) - grndmean
facA.mn.cntrst <- c(facA.mn.cntrst, mnsAA)
facB.mn.cntrst <- c(facB.mn.cntrst, mnsBB)

yy <- c(yy,mns.vec)
aaa<- paste(1:length(unique(A)), rep(paste(1:length(unique(B)), "mean", sep=""), ea = length(unique(A))), sep = "")

out <- list(signif(sort(cell.mnsA-grndmean),3), signif(sort(cell.mnsB-grndmean),3), signif(cell.cnts,3), signif(mns.matrx,3), summary(aov.yy))
names(out)<-c(paste(colnames(data.A.B[2]),".effects",sep=""),paste(colnames(data.A.B[3]),".effects",sep=""), 
  "CellCounts.Reordered", "CellMeans.Reordered", "aov.summary")

if(is.null(rownames(data.A.B))){rownames(data.A.B) <- 1:N}
if(is.null(offset)){
offset<-(((100/length(mnsA))^(1/3)) * 0.02)
}


scatter3d(facA.mn.cntrst, yy, facB.mn.cntrst, xlab = colnames(data.A.B)[2], ylab = colnames(data.A.B[1]), 
    zlab = colnames(data.A.B)[3], group = group.factor, fogtype='exp2',fov=55, surface = TRUE, fit = fit, surface.col = c(4,8), ...)
if(ident){
	if(is.null(offset)){offset<-((100/length(facA.mn.cntrst))^(1/3)) * 0.02}
	Identify3d(facA.mn.cntrst, yy, facB.mn.cntrst, labels = c(rownames(data.A.B), aaa), offset = offset)}

return(out)
}
