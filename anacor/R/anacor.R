`anacor` <-
function(tab, ndim = 2, row.covariates, col.covariates, scaling = c("Benzecri","Benzecri"), 
ellipse = TRUE, eps = 1e-6) 
{
#tab ... 2-way frequency table
#ndim ... number of dimensions
#row.covariates ... matrix A of dimension n x a
#col.covariatess ... matrix B of dimension m x b 
#scaling ... either "standard", "centroid", "Benzecri", "Goodman". 


#----------- sanity checks ---------------
tab <- as.matrix(tab)
if (length(scaling) == 1) scaling <- rep(scaling, 2) 
scaling[1] <- match.arg(scaling[1],c("standard", "centroid", "Benzecri", "Goodman"))
scaling[2] <- match.arg(scaling[2],c("standard", "centroid", "Benzecri", "Goodman"))
if (missing(row.covariates)) {
  row.covariates <- NULL
} else {
  row.covariates <- as.matrix(row.covariates)
  if (nrow(row.covariates) != nrow(tab)) stop("Matrix with row covariates does not match table dimensions!")
}
if (missing(col.covariates)) {
  col.covariates <- NULL 
} else {
  col.covariates <- as.matrix(col.covariates)
  if (nrow(col.covariates) != ncol(tab)) stop("Matrix with column covariates does not match table dimensions!")
}
if (ndim > min(dim(tab)) - 1) stop("Too many dimensions!")
#------------------ end sanity checks -----------------



name<-deparse(substitute(tab)) 

if (any(is.na(tab))) tab <- reconstitute(tab, eps = eps)				#NORA reconstitution of order 0 

n <- dim(tab)[1]							#number of rows
m <- dim(tab)[2]							#number of columns
N <- sum(tab)							#grand total (sample size)
tab <- as.matrix(tab)	
prop <- as.vector(t(tab))/N #relative frequencies 					
r <- rowSums(tab)							#row margins
c <- colSums(tab)							#column margins
qdim <- ndim+1
r <- ifelse(r==0,1,r)						#check 0 row/column margins and set to 1
c <- ifelse(c==0,1,c)
ROW <- FALSE
COL <- FALSE

z <- tab/sqrt(outer(r,c))						#D^(-1/2)*F*E^(-1/2) 
chisq <- N*(sum(z^2)-1)

if (is.matrix(row.covariates)) {
  xr <- cbind(rep(1,n),row.covariates)                          #add 1-vector to rowres-matrix
  nr <- dim(xr)[2]                                      #a+1
	cx<-crossprod(xr,r*xr)                              #t(xr)%*%(r*xr)                    
  or<-symSqrt(cx,inv=TRUE)                            #needed for computation of Z
  ROW<-TRUE
	}
if (is.matrix(col.covariates)) {
  yr<-cbind(rep(1,m),col.covariates)                          #add 1-vector to colres-matrix
  mr<-dim(yr)[2];
	cy<-crossprod(yr,c*yr)  
  oc<-symSqrt(cy,inv=TRUE)                            #needed for computation of Z
  COL<-TRUE
	}

# prepare expressions for SVD
if (ROW) {
  z <- or%*%crossprod(xr,tab) 
} else { 
  z <- tab/sqrt(r)
}

if (COL) {
  z <- z%*%yr%*%oc 
} else {
  z <- z/outer(rep(1,dim(z)[1]),sqrt(c))                         	#if ROW and COL false, z as above
}			

sv<-svd(z,nu=qdim,nv=qdim)					#SVD of Z --> $u lsv (P), $v rsv (Q), $d singular values
								#number of lsv and rsv is qdim = ndim+1
sval<-N*((sv$d[-1])^2)						#singval^2 times sample size (sv$d[1] trivial solution)	
sigmavec <-(sv$d)[2:qdim] 						#singval

dimlab <- paste("D", 1:ndim, sep = "")
Tr <- Xstar <- isetCorRow <- NULL                               #init Ter Braak 
Tc <- Ystar <- isetCorCol <- NULL
if (ROW) {
  x<-(xr%*%or%*%(sv$u))[,-1]                                    #canonial scores
  U <- or%*%(sv$u)[,-1]                                         #row weights
  Tr <- U%*%(diag(sigmavec))                                    #canonical coefficients (rows)
  dimnames(Tr) <- list(paste("beta",1:(dim(Tr)[1]), sep = ""), dimlab)
  #Tr2 <- solve(t(xr)%*%diag(r)%*%xr)%*%t(xr)%*%diag(r)%*%(diag(1/r)%*%tab%*%y)  #double check Tr (alternative formulation) 
  if (COL)  y <- (yr%*%oc%*%(sv$v))[,-1] else y <- ((sv$v)/sqrt(c))[,-1]       
  Xstar <- (diag(1/r)%*%tab%*%y)                                #site scores (rows), unscaled
  dimnames(Xstar) <- list(rownames(tab), dimlab)
  isetCorRow <- cor(Xstar,row.covariates)                       #instraset correlations for the rows
}
else {
  x <- ((sv$u)/sqrt(r))[,-1]                                    #compute row scores (without trivial solution)
 
}

if (COL) {
  y <- (yr%*%oc%*%(sv$v))[,-1]
  V <- oc%*%(sv$v)[,-1]                                         #column weights
  Tc <- V%*%(diag(sigmavec))                                    #canonical coefficients
  dimnames(Tc) <- list(paste("beta",1:(dim(Tc)[1]), sep = ""), dimlab)
  Ystar <- (diag(1/c)%*%t(tab)%*%x)                             #site scores (columns), unscaled
  dimnames(Ystar) <- list(colnames(tab), dimlab)
  isetCorCol <- cor(Ystar, col.covariates)                      #instraset correlations for the rows
} else {
  y <- ((sv$v)/sqrt(c))[,-1]                                  	#compute column scores (without trivial solution)
}

#x... matrix of dimension n x ndim with row scores; X = D^(-1/2)*P 
#y... matrix of dimension m x ndim with column scores; Y = E^(-1/2)*Q

x <- x*sqrt(N)								#standard coordinates
y <- y*sqrt(N)								#for each column holds: x'Dx = N, y'Ey=N

#--------------- rescale scores ----------------
if (scaling[2] == "centroid") y<-y*outer(rep(1,m),sigmavec)						
if (scaling[1] == "centroid") x <-x*outer(rep(1,n),sigmavec)						
if (scaling[2] == "Goodman") 	y<-y*outer(rep(1,m),sqrt(sigmavec))
if (scaling[1] == "Goodman")	x<-x*outer(rep(1,n),sqrt(sigmavec))
if (scaling[2] == "Benzecri")	y<-y*outer(rep(1,m),sigmavec)
if (scaling[1] == "Benzecri")   x<-x*outer(rep(1,n),sigmavec)

#compute Benzecri distances
benzres <- benzdist(scaling, x, y, z, tab, n, m, r, c, row.covariates, col.covariates)           
  
prob <- sval/chisq; pcum<-cumsum(prob)
cs.mat <- cbind(sval, prob, pcum)

#------------------------- generalized SVD and derivatives ---------------------
if (ellipse)
{
  res.deriv <- gsvdDer(tab,ndim)                             #Compute gradients
  res.deriv.scaling <- gsvdScal(res.deriv,scaling)           #Rescale gradients
  res.acov <- acovuv(res.deriv.scaling, prop, n, m, N, ndim) #VC-matrices
  se.sigma <- sqrt(diag(res.acov$acovd))                     #standard errors for singular values
} else {
  res.acov <- NULL
  se.sigma <- NULL
}
#------------------------ end generalized SVD and derivatives ------------------

#---------------------- labels ---------------------
colnames(x) <- colnames(y) <- dimlab
rownames(x) <- rownames(tab)
rownames(y) <- colnames(tab)
rownames(cs.mat) <- paste("Component", 1:dim(cs.mat)[1])
colnames(cs.mat) <- c("Chisq","Proportion","Cumulative Proportion")

                                               
result <- list(datname = name, tab = tab, ndim = ndim, row.covariates = row.covariates, 
               col.covariates = col.covariates, row.scores = x, col.scores = y, 
               chisq.decomp = cs.mat, chisq = chisq, singular.values = sigmavec, se.singular.values = se.sigma,
               left.singvec = sv$u[,-1], right.singvec = sv$v[,-1],
               eigen.values = sigmavec^2, scaling = scaling,
               bdmat = benzres[1:4], rmse = benzres[5:6], 
               row.acov = res.acov$acovu, col.acov = res.acov$acovv, cancoef = list(rows = Tr, columns = Tc),
               sitescores = list(rows = Xstar, columns = Ystar), isetcor = list(rows = isetCorRow, columns = isetCorCol)) 
class(result) <- "anacor"
result

}


