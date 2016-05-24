scree.plot <-
function(namefile, title="Scree Plot", type="R", use="complete.obs", simu="F")
{

mat <- namefile
if (use=="complete.obs") mat <- na.omit(namefile)

if (type=="R") eigenval <- eigen(cor(mat,use="pairwise.complete.obs"), symmetric=TRUE)$values
if (type=="V") eigenval <- eigen(cov(mat,use="pairwise.complete.obs"), symmetric=TRUE)$values
if (type=="E") eigenval <- namefile
if (type=="M") eigenval <- eigen(namefile, symmetric=TRUE)$values

nev <- length(eigenval)
plot(eigenval, type = "b", pch = 16, bty = "o", main = title, xlab = "Dimension", ylab = "Eigenvalue")
lines(c(1,nev),c(1,1),lty=2)	

if (is.numeric(simu) && (type=="R"))
{
n <- dim(mat)[1]
p <- dim(mat)[2]

matsimu <- matrix(nrow=n,ncol=p)
int <- rep(1,n*p)
attr(int,"dim") <- c(n,p)
mat <- pmax(as.matrix(mat),int)

for(i in 1:simu)
	{
	matnorm <- rnorm(n*p)
	attr(matnorm,"dim") <- c(n,p)
	matsimu <- (mat/mat)*matnorm
	eigenval <- eigen(cor(matsimu,use="pairwise.complete.obs"))$values
	points(eigenval,type="l")
	}
}
}
