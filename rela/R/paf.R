paf <- function(object, eigcrit=1, convcrit=.001) {

if (is.null(row.names(t(object)))) {print("Dataset must be coerced into matrix from prior data frame.")}

else {

if (!is.numeric(object)) {print("Your dataset is not a numeric object.")} 

else {

object<-(na.omit(object))
yy<-sort(abs(as.numeric(var(object))))
yy<-yy[1]
yy<-ifelse(yy==0, NA, yy)

if (is.na(yy)) {print("One of your variables is a constant. Constants are disallowed as part of a scale.")}

else {

if ( dim(cbind(object))[2]< 2 ) {print("Your set contains only one variable. At least two are required.")}

else {

if ( dim(cbind(object))[1]< 2 ) {print("Your set contains only one case. You need at least two cases.")}

else {


same <- function(x) {
for (i in 1:ncol(x)) {
	for (j in 1:ncol(x)) {
		if (i!=j) {
			test <- x[,i]-x[,j]
			if (sum(abs(test))==0) {print(paste("WARNING: Items ",i," (",colnames(x)[i],") ","and ",j," (",colnames(x)[j],") ","are duplicates."))} 
			}
		}
	}
}

same(object)


x <- na.omit(object)
x <- cor(x)

N <- nrow(na.omit(object))

options(digits=5)

x <- as.matrix(x)
x1 <- x; rownames(x1)=colnames(x)

bt<- -((N-1)-((2*ncol(x)+5)/6))*log(det(x))

# Correlation matrix inverse
invx <- solve(x)

# S-squared matrix
ssqr <- matrix(0,nrow(x),ncol(x))
diag(ssqr) <- diag(invx)
ssqr <- solve(ssqr)

# Partial correlations (anti-image)
Q <- ssqr%*%invx%*%ssqr
colnames(Q) <- colnames(x)
rownames(Q) <- colnames(x)

# Anti-image correlation matrix
qdi <- diag(diag(Q))
sqdi <- solve(qdi)
ssqdi <- sqrt(sqdi)
Qr <- ssqdi%*%Q%*%ssqdi
colnames(Qr) <- colnames(x)
rownames(Qr) <- colnames(x)

# KMO
xnod <- x
diag(xnod) <- 0
Qrnod <- Qr
diag(Qrnod) <- 0
KMO <- sum(xnod^2)/(sum(xnod^2)+sum(Qrnod^2))

# MSA
MSA <- 0
for (i in 1:nrow(xnod)) (MSA<-c(MSA,sum(xnod[i,]^2)/(sum(xnod[i,]^2)+sum(Qrnod[i,]^2))))
MSA<- matrix(MSA[-1],,1)
rownames(MSA) <- colnames(x)
colnames(MSA) <- "MSA"

# First communialities
comm0 <- 1-diag(Q)

# Computing subsequent iterrations
comm1 <- comm0
allcomm <- 0
diffscores <- 0
iter <- 0
count <- 0
eigenval <- 0
x0 <- x

repeat {

allcomm <- cbind(allcomm,comm1)

eigs <- 0
diag(x) <- comm1

for (i in 1:length(eigen(x)$values))
if (eigen(x0)$values[i] > eigcrit) {eigs <- c(eigs,eigen(x)$values[i]) }

eigs <- eigs[-1]
eigenval <- cbind(eigenval, eigen(x)$values)

eigmat <- sqrt(diag(eigs, length(eigs), length(eigs)))

eigvecr <- matrix(eigen(x)$vector[,0:length(eigs)],,length(eigs))
one <- c((1:ncol(eigvecr))/(1:ncol(eigvecr)))

factload <- eigvecr%*%eigmat

comm2 <- t(rowsum(t(factload^2), one))
dif <- abs(comm1-comm2)

iter <- iter+1
count <- c(count,iter)
diffscores <- cbind(diffscores, dif)

comm1 <- comm2

endtest <- matrix(1, nrow(dif),1)
for (i in 1:nrow(dif)) if (dif[i,1]<convcrit) {endtest[i,1]=NA}
if (length(na.omit(endtest))==0) break }

# Preparing and labeling output
firstlast <- cbind(comm0,comm1); colnames(firstlast)=c("Initial Communalities", "Final Extraction")
allcomm <- cbind(allcomm, comm1)
allcomm <- allcomm[,-1]; colnames(allcomm)=count
diffscores <- diffscores[,-1]; colnames(diffscores)=c(1:iter)
eigenval <- cbind(eigen(x0)$values, eigenval[,-1]); colnames(eigenval)=c(0:iter)

rownames(factload) <- colnames(x)
facttest <- factload[,1]
for (i in 1:length(facttest)) if (facttest[i]<0) {facttest[i]<-NA}
if (length(na.omit(facttest))==0) (factload[,1] <- -factload[,1])

correp <- factload %*% t(factload)
residuals <- x-correp
colnames(correp) <- colnames(x)
rownames(correp) <- colnames(x)
colnames(residuals) <- colnames(x)
rownames(residuals) <- colnames(x)
diag(correp) <- 1
diag(residuals) <- 0
RMS <- sqrt(sum(residuals^2)/((nrow(x)^2)-(nrow(x))))

output <- list("PRINCIPAL AXIS FACTORING", Correlation=x1, Anti.Image.Cov=Q,
Anti.Image.Cor=Qr, KMO=KMO, MSA=MSA, 
Bartlett=bt, Communalities=firstlast, Iterations=iter, Eigenvalues=eigenval, 
Communality.Iterations=allcomm, Criterion.Differences=diffscores, Factor.Loadings=factload, Reproduced.Cor=correp,
Residuals=residuals, RMS=RMS)

ob <- as.character(match.call()[2])
cl <- call("paf","object"=ob,"eigcrit"=eigcrit, "convcrit"=convcrit)

output$call <- cl
output$items <- names(output)

class(output) <- c("paf", class(output))

output

}}}}}}
