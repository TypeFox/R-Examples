#' @title GWAS for Studies having Repeated Measurements on Related Individuals
#'
#' @description
#'  It is used to perform genome-wide association studies on individuals that are both related and have repeated measurements. 
#'  The function computes score statistic based p-values for a linear mixed model including random polygenic effects and 
#'  a random effect for repeated measurements. A p-value is computed for each marker and the null hypothesis tested is a
#'   zero additive marker effect. 
#'
#' @param formula.FixedEffects Formula including the response variable and cofactors as fixed effects.
#' @param genabel.data An GenABEL object including marker information. This object has one observtion per individuals. 
#' @param phenotype.data A data frame including the repeated observations and IDs. 
#' @param id.name The column name of the IDs in phen.data
#' @param GRM An optional genetic relationship matrix (GRM) can be included as input. Otherwise the GRM is computed within the function.
#' @param V An optional (co)variance matrix can be included as input. Otherwise it is computed using the hglm function.
#' @param memory Used to optimize computations. The maximum number of elements in a matrix that can be stored efficiently.
#' @details 
#' A generalized squares (GLS) is fitted for each marker given a (co)variance matrix V. 
#' The computations are made fast by transforming the GLS to 
#' an ordinary least-squares (OLS) problem using an eigen-decomposition of V. 
#' The OLS are computed using QR-factorization. If V is not specified then a model 
#' including random polygenic effects and permanent environmental effects is 
#' fitted (using the hglm package) to compute V. A GenABEL object (scan.gwaa class) 
#' is returned (including also the \code{hglm} results).
#' Let e.g. GWAS1 be an object returned by the \code{rGLS} function. 
#' Then a Manhattan plot can be produced by calling \code{plot(GWAS1)} and 
#' the top SNPs using \code{summary(GWAS1)}. Both of these functions are 
#' generic GenABEL functions. \cr
#' The results from the fitted linear mixed model without any SNP effect included 
#' are produced by calling \code{summary(GWAS1@@call$hglm)}.
#' 
#' @author Lars Ronnegard
#' 
#' 
#' @examples
#'  data(Phen.Data) #Phenotype data with repeated observations
#'  data(gen.data) #GenABEL object including IDs and marker genotypes
#'  GWAS1 <- rGLS(y ~ age + sex, genabel.data = gen.data, phenotype.data = Phen.Data) 
#'  plot(GWAS1, main="")
#'  summary(GWAS1)
#'  #Summary for variance component estimation without SNP effects
#'  summary(GWAS1@@call$hglm) 
#'
rGLS <-
function(formula.FixedEffects = y ~ 1, genabel.data, phenotype.data, id.name = "id", GRM = NULL, V = NULL, memory=1e8) {
	#Check input data
	if (class(genabel.data) != "gwaa.data") stop("The input of genabel.data is not a GenABEL object")
	if (is.null(genabel.data@phdata$id)) stop("IDs not given as id in the phdata list")
    if (!is.null(GRM)) { if(!isSymmetric(GRM)) warning("The given GRM must be a symmetric matrix!") }
    if (!is.null(V)) { if(!isSymmetric(V)) warning("The given V must be a symmetric matrix!") }
    V.input  <- V
    #require(hglm)
    #require(GenABEL)
    #Get trait name
	trait <- all.vars(formula.FixedEffects)[1]
    #Remove NAs from phenotypic data
    y.all <- phenotype.data[,names(phenotype.data)%in%trait]
    phenotype.data <- phenotype.data[!is.na(y.all),]
	#Connect IDs in GenABEL data set with IDs in the phenotype file
	id1 <- phenotype.data[,names(phenotype.data) %in% id.name] #ID for phenotype data
	id2 <- genabel.data@phdata$id #ID for genotype data
    test1 <- id1 %in% id2
    test2 <- id2 %in% id1
    genabel.data <- genabel.data[test2,] #Exclude individuals having no phenotype information
    phenotype.data <- phenotype.data[test1,] #Exclude individuals having no genotype information
    id1 <- phenotype.data[,names(phenotype.data) %in% id.name] #ID for phenotype data for cleaned data
	id2 <- genabel.data@phdata$id #ID for genotype data for cleaned data
	#####################
	#Construct incidence matrix for repeated observations
    N=length(id2)
    n=length(id1)
    indx <- numeric(n)
    for (i in 1:N) {
        indx <- indx + i * (id1 %in% id2[i])
    }
    Z.indx <- diag(N)[indx,]
    #Construct response and design matrix
	y <- phenotype.data[ , names(phenotype.data) %in% trait] #Create the response variable
    X <- model.matrix(formula.FixedEffects, data = phenotype.data) #Fixed effect design matrix
	#####################
 if (is.null(V)) {
    #Construct GRM
    if (is.null(GRM)) {
      autosomalMarkers <- which(chromosome(genabel.data)!= "X")
      GRM <- compute.GRM(genabel.data[ , snpnames(genabel.data)[autosomalMarkers]])
    }
    eig <- eigen(GRM)
    if (max(diag(GRM)) > 1.6) print("There seems to be highly inbred individuals in your data")
    if (min(eig$values < -0.5)) print("The genetic relationship matrix is far from positive definite")
    non_zero.eigenvalues <- eig$values>(1e-6) #Put numerically small eigenvalues to zero
    eig$values[ !non_zero.eigenvalues ] <- 0
    print("GRM ready")
    #####################
    #Fit hglm
    Z.GRM <- ( eig$vectors %*% diag(sqrt(eig$values)) )[indx, ]
    Z <- (cbind(Z.GRM, Z.indx))
    mod1 <- hglm(y=y, X=X, Z=Z, RandC = c(ncol(Z.GRM), ncol(Z.indx)), maxit = 200)
    if (mod1$Converge != "converged") stop("The variance component estimation did not converge in 200 iterations. Try to estimate them separately and provide the estimated (co)variance matrix V as input. \n\n")
    print("Variance component estimation ready")
    #####################
    #Construct rotation matrix
    ratio <- mod1$varRanef/mod1$varFix
    V <- constructV(Z=Z, RandC = c(ncol(Z.GRM), ncol(Z.indx)), ratio)
  }
	eig.V <- eigen(V)
	transf.matrix <- diag(1/sqrt(eig.V$values)) %*% t(eig.V$vectors)
	y.new <- transf.matrix %*% y
	X.new <- transf.matrix %*% X
	print("Rotation matrix ready")
	#####################
	#Fit a linear model for each SNP
	SNP.matrix <- as.double(genabel.data)
    if (sum(is.na(SNP.matrix)) > 0) {
        SNP.matrix <- SmoothSNPmatrix(SNP.matrix)
    }
	m <- ncol(SNP.matrix)
	p.val <- SNP.est <- rep(1, m)
	colnames(X.new) <- as.character(1:ncol(X.new)) #To avoid columns having strange names
	print("Rotate LMM started")
    #Fit using QR factorization
    #Null model
	qr0 <- qr(X.new)
	est0 <- qr.coef(qr0, y.new)
	res <- y.new - X.new %*% est0
    n <- length(y.new)
	RSS.0 <- sum(res^2)/n
	#Split computations into reasonably sized blocks
    if (memory < n) memory <- n
	step.size <- floor(memory/n)
    steps <- ceiling(m/step.size)
	jj=1
	kk=0
	for (step.i in 1:steps) {
		if (step.i == steps) kk <- kk + m %% step.size else kk <- kk + step.size
		markers.to.fit <- jj:kk
		snp.new <- transf.matrix %*% SNP.matrix[indx, markers.to.fit]
		mm <- 0
		for (j in markers.to.fit) {
			mm <- mm+1
			X1 <- cbind(snp.new[, mm], X.new)
			qr1 <- qr(X1)
			est1 <- qr.coef(qr1, y.new)
			res <- y.new - X1 %*% est1
			RSS.1 <- sum(res^2)/n
			SNP.est[j] <- est1[1]
			LRT <- -n * (log(RSS.1) - log(RSS.0))
			p.val[j] <- 1 - pchisq(LRT, df=1)
        }
        jj <- jj + step.size
	}
	print("Rotate LMM ready")
	#####################
	qt.results <- Create_gwaa_scan(genabel.data, p.val, SNP.est)
	if (is.null(V.input)) qt.results@call$hglm <- mod1
	return(qt.results)
}
