AminoAcids = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
AAbyGroup = c("aliphatic", "cysteine", "acidic", "acidic", "aromatic", "aliphatic", "basic", "aliphatic", "basic", "aliphatic", "aliphatic", "aminic", "proline", "aminic", "basic", "hydroxylated", "hydroxylated", "aliphatic", "aromatic", "aromatic")
AAGroups = c("acidic","aliphatic", "aminic",  "aromatic", "basic","cysteine","hydroxylated", "proline")
 small = c(1,2,3,6,12,13,16,17,18)
 polar = c(2,3,4,7,9,12,14,15,16,17,19,20)
 hydrophobic = c(1, 2,5,6,7,8,9,10,11,17,18,19,20)


#############

# Factor Scores as computed by SAS Factor Analysis of 54 Amino Acid Indices 
AAMetric.Atchley = matrix( 
	c(
 -0.59145974, -1.30209266, -0.7330651,  1.5703918, -0.14550842,
 -1.34267179,  0.46542300, -0.8620345, -1.0200786, -0.25516894,
  1.05015062,  0.30242411, -3.6559147, -0.2590236, -3.24176791,
  1.35733226, -1.45275578,  1.4766610,  0.1129444, -0.83715681,
 -1.00610084, -0.59046634,  1.8909687, -0.3966186,  0.41194139,
 -0.38387987,  1.65201497,  1.3301017,  1.0449765,  2.06385566,
  0.33616543, -0.41662780, -1.6733690, -1.4738898, -0.07772917,
 -1.23936304, -0.54652238,  2.1314349,  0.3931618,  0.81630366,
  1.83146558, -0.56109831,  0.5332237, -0.2771101,  1.64762794,
 -1.01895162, -0.98693471, -1.5046185,  1.2658296, -0.91181195,
 -0.66312569, -1.52353917,  2.2194787, -1.0047207,  1.21181214,
  0.94535614,  0.82846219,  1.2991286, -0.1688162,  0.93339498,
  0.18862522,  2.08084151, -1.6283286,  0.4207004, -1.39177378,
  0.93056541, -0.17926549, -3.0048731, -0.5025910, -1.85303476,
  1.53754853, -0.05472897,  1.5021086,  0.4403185,  2.89744417,
 -0.22788299,  1.39869991, -4.7596375,  0.6701745, -2.64747356,
 -0.03181782,  0.32571153,  2.2134612,  0.9078985,  1.31337035,
 -1.33661279, -0.27854634, -0.5440132,  1.2419935, -1.26225362,
 -0.59533918,  0.00907760,  0.6719274, -2.1275244, -0.18358096,
  0.25999617,  0.82992312,  3.0973596, -0.8380164,  1.51150958
    ), nrow=20, ncol=5, byrow=TRUE)
   
#Factor Scores as computed by R Factor Analysis (factor.pa.ginv) of 54 Amino Acid Indices   
AAMetric = matrix( 
	c(
 -0.60677211, -1.08940440, -0.92798080,  1.4726202,  0.251417050,
 -1.25436934,  0.46126134, -0.35464046, -1.1040471, -0.115742254,
  1.20123120,  0.24369180, -1.47196724, -0.3243624,  2.127034065,
  1.29519207, -1.46537642, -0.72868426,  0.1482127,  2.078732004,
 -0.95950963, -0.43941703,  1.15581211, -0.5634368, -0.175433772,
 -0.52784144,  1.75873553, -2.14269672,  0.7848519, -0.002451762,
  0.32715696, -0.52087575, -0.42326104, -1.4308910, -0.494428309,
 -1.27633548, -0.59487265,  1.03001883,  0.4256597,  0.014711047,
  1.82248357, -0.53137304, -0.02567446, -0.2432211, -1.571322338,
 -0.97854190, -1.12494822,  0.87095396,  1.4290018, -0.183751797,
 -0.72854024, -1.41768155,  0.30142692, -1.0227830, -0.094615066,
  0.89146155,  1.11359534, -0.85232822, -0.4088754,  0.246465610,
  0.21811094,  1.96849543, -0.15443545,  0.3903353,  0.472209055,
  0.90061235, -0.53619651,  0.02364059, -0.4459318,  0.234557443,
  1.52748331, -0.01654753,  0.65573629,  0.2563726, -2.460536571,
  0.03475833,  0.97046643, -0.84138070,  0.9825769,  0.108360575,
 -0.01042641,  0.65604029,  0.01085433,  0.8716081,  0.093194098,
 -1.32574165, -0.37200927,  1.01019249,  1.4746391,  0.121129607,
 -0.61057537,  0.01937737,  1.50913663, -1.8958892, -0.441521905,
  0.06016331,  0.91703885,  1.35527720, -0.7964405, -0.208006781
	), nrow=20, ncol=5, byrow=TRUE)

colnames(AAMetric) = c("pah", "pss", "ms", "cc", "ec")
rownames(AAMetric) = AminoAcids

colnames(AAMetric.Atchley) = c("pah", "pss", "ms", "cc", "ec")
rownames(AAMetric.Atchley) = AminoAcids

#####################

Loadings.variation = 
  function(sdev, digits = 5){

	var = sdev^2

	PTV = var / sum(var)		#proportion var
	CTV = cumsum(PTV)				#cumulative var
	TV = rbind(Lambda=round(var, digits), PTV=round(PTV, digits), CTV=round(CTV, digits))
	TV
}



FactorTransform = 
 function(Source, Search= AminoAcids, Replace = AAMetric.Atchley, Factor = 1, bycol = TRUE, SeqName = NULL, alignment=FALSE, fillblank=NA){

	if(is.matrix(Replace)){
		if(bycol){
			if(ncol(Replace) < Factor) 
				stop("Column for Factor is greater than size of Replacement Matrix")
			Replace = as.vector(t(Replace[,Factor]))
		}
		else{
			if(nrow(Replace) < Factor) 
				stop("Row for Factor is greater than size of Replacement Matrix")
			Replace = as.vector(Replace[Factor,])
		}
	}


	if(is.list(Source)){

		if(!missing(SeqName)){
			if(length(SeqName) != length(Source)){
				print("Length of SeqName not the same as length of Source. Using default.")
				SeqName = NULL
			}
			else{ SeqNames = SeqName}
		}
		else if (is.null(SeqName)){
			if(is.null(names(Source)))
				SeqNames = apply(array(1:length(Source)), 1, function(i) paste("Seq",i, sep="")) 			
			else
				SeqNames = names(Source)		
		}
		
		result = FactorTransform.vector(as.character(Source), Search, Replace)			
		names(result) = SeqNames
	}
	else if(is.vector(Source)){
		if(!missing(SeqName)){
			if(length(SeqName) != length(Source)){
				print("Length of SeqName not the same as length of Source. Using default.")
				SeqName = NULL
			}
			else{ SeqNames = SeqName}
		}
		else if (is.null(SeqName)){
			SeqNames = apply(array(1:length(Source)), 1, function(i) paste("Seq",i, sep="")) 
		}
		
		result = FactorTransform.vector(Source, Search, Replace)
		if(length(SeqNames)==1)
			result = list(Seq1=result)
		else
			names(result) = SeqNames
	}
	else if(is.array(Source)){
		result = FactorTransform.vector(Source, Search, Replace)
		rownames(result) = rownames(Source)
		colnames(result) = colnames(Source)	
	}
	else if(is.matrix(Source)){
		result= matrix(nrow= nrow(Source), ncol = ncol(Source))
		for(row in 1:nrow(Source)){
			Seq = Source[row,]
			result[row,] = FactorTransform.default(Seq, Search, Replace)
		}
		rownames(result) = rownames(Source)
		colnames(result) = colnames(Source)
	}
	else{
		stop("Source must be a list, vector, matrix, or array")
	}

	if(alignment){	#create matrix
		maxlength = max(sapply(result, length))
		result = t(sapply(result, function(x){append(x, rep(fillblank, times = maxlength-length(x)))}))
##		result = matrix(unlist(result), nrow = length(result), byrow = TRUE, dimnames = list(names(result)))
	}
	
	result
}

FactorTransform.vector = 
function(Source, Search, Replace){

	if(length(Source) ==1){
		Source = unlist(strsplit(Source, NULL, fixed=TRUE))
		SeqMetric = FactorTransform.default(Source, Search, Replace)
	}
	else{
		SeqMetric = lapply(Source, FactorTransform.vector, Search, Replace)			#transform each seq in vector
	}
	SeqMetric
}

FactorTransform.default=
function(Source, Search=AminoAcids, Replace=AAMetric.Atchley)
{
	if(!is.vector(Source))
		stop("Source is not of type vector\n")

    if (length(Search) != length(Replace))
      	stop("Search and Replace Must Have Equal Number of Items\n")
	
	Changed <- as.character(Source)

	for (i in 1:length(Search))
   		Changed <- replace(Changed, Changed == Search[i], Replace[i])
    
  as.numeric(Changed)
}
## adapted from Marc Schwartz code at https://stat.ethz.ch/pipermail/r-help/2006-July/108829.html


#################################
factor.pa.ginv = function (r, nfactors = 1, residuals = FALSE, prerotate= FALSE, rotate = "varimax",  m=4,
				   n.obs = NA, scores = c("none", "regression", "Bartlett"), force=FALSE, SMC = TRUE,  missing = FALSE, 
			       impute = "median",  min.err = 0.001,  digits = 2, max.iter = 50, symmetric = TRUE, warnings = TRUE
				  ) 
{
	if(missing(scores)) scores = "none"
	x.matrix = NULL
	
    n <- dim(r)[2]
    if (n != dim(r)[1]) {
        n.obs <- dim(r)[1]
        if (scores !="none") {
            x.matrix <- r
            if (missing) {
                miss <- which(is.na(x.matrix), arr.ind = TRUE)
                if (impute == "mean") {
                  item.means <- colMeans(x.matrix, na.rm = TRUE)
                  x.matrix[miss] <- item.means[miss[, 2]]
                }
                else {
                  item.med <- apply(x.matrix, 2, median, na.rm = TRUE)
                  x.matrix[miss] <- item.med[miss[, 2]]
                }
            }
        }
        r <- cor(r, use = "pairwise")
    }
 
   else{
   
	   if(scores != "none"){
 	  	if(force){						#data matrix square, can still compute scores
   			n.obs <- dim(r)[1]
        	x.matrix <- r
	        if (missing) {
				miss <- which(is.na(x.matrix), arr.ind = TRUE)
        	    if (impute == "mean") {
            		item.means <- colMeans(x.matrix, na.rm = TRUE)
                	x.matrix[miss] <- item.means[miss[, 2]]
	            }
    	        else {
        	    	item.med <- apply(x.matrix, 2, median, na.rm = TRUE)
            	    x.matrix[miss] <- item.med[miss[, 2]]
            	}
	        }
    	    r <- cor(r, use = "pairwise")
   		}
	   	else{
    		if(warnings) message("Data structure is square and assumed to be covariance matrix. Cannot compute scores.")
    		scores = "none"
    		}
  	 }
   
   	if (!is.matrix(r)) {
       r <- as.matrix(r)
    }  
    
    sds <- sqrt(diag(r))
    r <- r/(sds %o% sds)											# convert to correlation matrix
  }
  
     
    r.mat <- r
    Phi <- NULL
    colnames(r.mat) <- rownames(r.mat) <- colnames(r)

	orig = diag(r.mat)							# remove invariable variables
    diag(r.mat) = NA
    allNA = apply(is.na(r.mat), 1, all)
    allNA.names = names(allNA[allNA==TRUE])
    if(length(allNA.names)>0){
    	if(warnings) message("Invariables removed from analysis: ", paste(allNA.names, collapse=", "))
		r.mat <- r.mat[!apply(is.na(r.mat), 1, all),]
		r.mat <- r.mat[,!apply(is.na(r.mat), 2, all)]
		r <- r[!allNA,]
		r <- r[,!allNA]
		orig <- orig[!allNA]
		n = n - length(allNA.names)
		if(!is.null(x.matrix)){x.matrix= x.matrix[,!allNA]}
	}
	diag(r.mat) = orig			

  if (!residuals) {
        result <- list(values = c(rep(0, n)), rotation = rotate, 
            n.obs = n.obs, communality = c(rep(0, n)), loadings = matrix(rep(0, 
                n * n), ncol = n), fit = 0)
    }
    else {
        result <- list(values = c(rep(0, n)), rotation = rotate, 
            n.obs = n.obs, communality = c(rep(0, n)), loadings = matrix(rep(0, 
                n * n), ncol = n), residual = matrix(rep(0, n * 
                n), ncol = n), fit = 0)
    }

    if (SMC) {
        if (nfactors < n/2) {
            diag(r.mat) <- smc(r)
        }
        else {
            if (warnings) 
                message("too many factors requested for this number of variables to use SMC, 1s used instead")
        }
    }
    orig <- diag(r)
    comm <- sum(diag(r.mat))
    err <- comm
    i <- 1
    comm.list <- list()
    while (err > min.err) {
        eigens <- eigen(r.mat, symmetric = symmetric)
        if (nfactors > 1) {
            loadings <- eigens$vectors[, 1:nfactors] %*% diag(sqrt(eigens$values[1:nfactors]))
        }
        else {
            loadings <- eigens$vectors[, 1] * sqrt(eigens$values[1])
        }
        model <- loadings %*% t(loadings)
        new <- diag(model)
        comm1 <- sum(new)
        diag(r.mat) <- new
        err <- abs(comm - comm1)
        if (is.na(err)) {
            warning("imaginary eigen value condition encountered in factor.pa,\n Try again with SMC=FALSE \n exiting factor.pa")
            break
        }
        comm <- comm1
        comm.list[[i]] <- comm1
        i <- i + 1
        if (i > max.iter) {
            if (warnings) {
                message("maximum iteration exceeded")
            }
            err <- 0
        }
    }
    if (!is.numeric(loadings)) {
        warning("the matrix has produced imaginary results -- proceed with caution")
        loadings <- matrix(as.double(loadings), ncol = nfactors)
    }
    
    if (nfactors > 1) {
        sign.tot <- vector(mode = "numeric", length = nfactors)
        sign.tot <- colSums(loadings)
        sign.tot = sign(sign.tot)
        
        loadings <- loadings %*% diag(sign.tot)
        
    }
    else {
        if (sum(loadings) < 0) {
            loadings <- -as.matrix(loadings)
        }
        else {
            loadings <- as.matrix(loadings)
        }
        colnames(loadings) <- "F1"
    }
    colnames(loadings) <- paste("F", 1:nfactors, sep = "")
    rownames(loadings) <- rownames(r)
    loadings[loadings == 0] <- 10^-15
    rotmat = diag(1, nfactors)
    model <- loadings %*% t(loadings)
    result$communality <- round(diag(model), digits)
    result$uniquenesses <- round(diag(r - model), digits)
    
    if (rotate != "none") {
        if (nfactors > 1) {
        	if(prerotate){
        		varmax = varimax(loadings) 
        		
        		 ev.rotated <- diag(t(varmax$loadings) %*% varmax$loadings)
			     ev.order <- order(ev.rotated, decreasing = TRUE)
			     loadings <- varmax$loadings[, ev.order]
			     colnames(loadings) <- colnames(loadings[,ev.order])
    	
    			 rotmat <- varmax$rotmat[, ev.order]
              	 Phi <- NULL
        	}
            if ((rotate == "varimax" | rotate == "Varimax")  && !prerotate) {
            	varmax = varimax(loadings)
            	
            	 ev.rotated <- diag(t(varmax$loadings) %*% varmax$loadings)
			     ev.order <- order(ev.rotated, decreasing = TRUE)
			     loadings <- varmax$loadings[, ev.order]
			     colnames(loadings) <- colnames(loadings[,ev.order])
    	
    			 rotmat <- varmax$rotmat[, ev.order]
           	           	
               	Phi <- NULL
            }
            else {
                if ((rotate == "promax") | (rotate == "Promax")) {
                              	
                  pro <- Promax.only(loadings, m, rotmat)
                  loadings <- pro$loadings
                  
                  rotmat = pro$rotmat
                  Phi <- pro$Phi
                }
            }
        }
    }
    if (nfactors > 1) {
        ev.rotated <- diag(t(loadings) %*% loadings)
        ev.order <- order(ev.rotated, decreasing = TRUE)
        loadings <- loadings[, ev.order]
    }
    colnames(loadings) <- paste("F", 1:nfactors, sep = "")
    rownames(loadings) <- colnames(r)
    if (!is.null(Phi)) {
        Phi <- Phi[ev.order, ev.order]
    }
    class(loadings) <- "loadings"
    if (nfactors < 1) 
        nfactors <- n
    residual <- r - model
    r2 <- sum(r * r)
    rstar2 <- sum(residual * residual)
    if (residuals) {
        result$residual <- round(residual, digits)
    }
    r2.off <- r
    diag(r2.off) <- 0
    r2.off <- sum(r2.off^2)
    diag(residual) <- 0
    rstar.off <- sum(residual^2)
    result$fit <- round(1 - rstar2/r2, digits)
    result$fit.off <- round(1 - rstar.off/r2.off, digits)
    result$values <- round(eigens$values, digits)
    diag(model) <- 1
    model.inv <- solve(model)
    m.inv.r <- model.inv %*% r
    if (is.na(n.obs)) {
        result$n.obs = NA
        result$PVAL = NA
    }
    else {
        result$n.obs = n.obs
    }
    
    result$rotmat = rotmat
  
  	result$dof <- n * (n - 1)/2 - n * nfactors + (nfactors * 
        (nfactors - 1)/2)
    result$objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) - 
        n
    result$criteria <- c(objective = result$objective, NA, NA)
    if (!is.na(n.obs)) {
        result$STATISTIC <- result$objective * ((n.obs - 1) - 
            (2 * n + 5)/6 - (2 * nfactors)/3)
        if (!is.nan(result$STATISTIC)) 
            if (result$STATISTIC < 0) {
                result$STATISTIC <- 0
            }
        if (result$dof > 0) {
            result$PVAL <- pchisq(result$STATISTIC, result$dof, 
                lower.tail = FALSE)
        }
        else result$PVAL <- NA
    }
    result$loadings <- round(loadings, digits)
    if (!is.null(Phi)) {
        result$Phi <- Phi
    }
    result$communality.iterations <- round(unlist(comm.list), digits)
    if (scores != "none") {
       
       Lambda <- result$loadings
        zz <- scale(x.matrix, TRUE, TRUE)
        
		switch(scores, 
		 regression = {
		
				InvCov =  try(solve(r, Lambda), TRUE) 
				if(length(InvCov)<=1){
				  	cat("Could not solve for inverse correlation.  Using general inverse ginv(r) \n")
				 	InvCov = ginv(r) %*% Lambda
         	    }
           	    sc <- zz %*% InvCov
            
    	        if (!is.null(Phi)) {
                 	sc <- sc %*% Phi
				}            
     }, Bartlett = {
            d <- 1/(result$uniquenesses)
            tmp <- t(Lambda * d)
            sc <- t(solve(tmp %*% Lambda, tmp %*% t(zz)))
        })

        rownames(sc) <- rownames(x.matrix)
        colnames(sc) <- colnames(Lambda)
        
        result$scores <- sc
    }
    result$factors <- nfactors
    result$fn <- "factor.pa"
    class(result) <- c("psych", "fa")
    return(result)
}

#################

Promax.only = 
function (x, m = 4, rotate.structure=NULL) 
{
    if (!is.matrix(x) & !is.data.frame(x)) {
        if (!is.null(x$loadings)) 
            x <- as.matrix(x$loadings)
    }
    else {
        x <- x
    }
 
    if (ncol(x) < 2) 
        return(x)
    dn <- dimnames(x)
         
    if(missing(rotate.structure)){
    	rotate.structure= diag(1, ncol(x))
    }
         
    Q <- x * abs(x)^(m - 1)
    U <- lm.fit(x, Q)$coefficients
    d <- diag(solve(t(U) %*% U))
    U <- U %*% diag(sqrt(d))
    dimnames(U) <- NULL
    z <- x %*% U
    U <- rotate.structure %*% U
    ui <- solve(U)
    Phi <- ui %*% t(ui)
    dimnames(z) <- dn
    class(z) <- "loadings"
    result <- list(loadings = z, rotmat = U, Phi = Phi)
    class(result) <- c("HDMD")
    return(result)
}


###################################

pairwise.mahalanobis = 
function(x, grouping=NULL, cov =NULL, inverted=FALSE, digits = 5, ...) 
{
	x <- if (is.vector(x)) 									#standardize input data as matrix
    	matrix(x, ncol = length(x))
    else as.matrix(x)

	if(!is.matrix(x))
		stop("x could not be forced into a matrix")

 	if(length(grouping) ==0){								#no group assigned, uses first col
		grouping = t(x[1])
		x = x[2:dim(x)[2]] 	
		cat("assigning grouping\n")
		print(grouping)
 	}

	n <- nrow(x)
	p <- ncol(x)

	if (n != length(grouping)){ 								#grouping and matrix do not correspond
		cat(paste("n: ", n, "and groups: ", length(grouping), "\n"))
        stop("nrow(x) and length(grouping) are different")
	}
    g <- as.factor(grouping)
    g
    lev <- lev1 <- levels(g)								# Groups
    counts <- as.vector(table(g))							# No. of elements in each group

	 if (any(counts == 0)) {								# Remove grouping if not represented in data
        empty <- lev[counts == 0]
        warning(sprintf(ngettext(length(empty), "group %s is empty", 
            "groups %s are empty"), paste(empty, collapse = " ")), 
            domain = NA)
        lev1 <- lev[counts > 0]
        g <- factor(g, levels = lev1)
        counts <- as.vector(table(g))
    }

	ng = length(lev1)

	group.means <- tapply(x, list(rep(g, p), col(x)), mean)		# g x p matrix of group means from x
	
	if(missing(cov)){											#create covariance matrix
		inverted = FALSE
		cov = cor(x)											#standardize into correlation mtx
	}
	else{														#check cov of correct dimension
		if(dim(cov) != c(p,p))
			stop("cov matrix not of dim = (p,p)\n")
	}

	Distance = matrix(nrow=ng, ncol=ng)									#initialize distance matrix	
	dimnames(Distance) = list(names(group.means), names(group.means))

	Means = round(group.means, digits)
	Cov = round(cov, digits)
	Distance = round(Distance, digits)
	
	for(i in 1:ng){
		Distance[i,]= mahalanobis(group.means, group.means[i,], cov, inverted)
	}

	result <- list(means = group.means, cov = cov, distance = Distance)
	result
}


######################

MolecularEntropy = 
function(x, type){

	H=NULL
	
	if(!is.matrix(x)){
		if(is.vector(x)){ 
			x = strsplit(x, "")
			x = matrix(unlist(x), nrow = length(x), byrow = TRUE, dimnames = list(names(x)))
		}
		if(is.list(x)){
			x = sapply(x, strsplit, "")
			x = matrix(unlist(x), nrow = length(x), byrow = TRUE, dimnames = list(names(x)))
		}
	}

    if (type == "DNA") {
    	DNA = c("A", "C", "G", "T")
    	if(!all(x %in% DNA)){print("Warning: Data set contains non-nucleotide elements")}
    	
	    counts = apply(x, 2, function(z){table(factor(z, levels = DNA))})
	
	    freq = apply(counts, 2, function(freqs){if(sum(freqs) > 0){freqs=freqs/sum(freqs) } else{freqs=freqs/1}  })
		
		H = apply(freq, 2, function(freqs){-sum(ifelse(freqs > 0, freqs * log(freqs) , 0)) })
		H = H/log(4)

    }
 
    if (type == "AA") {
    	if(!all(x %in% AminoAcids)){print("Warning: Data set contains non-Amino Acid elements")}

		counts = apply(x, 2, function(z){table(factor(z, levels = AminoAcids))})

	    freq = apply(counts, 2, function(freqs){if(sum(freqs) > 0){freqs=freqs/sum(freqs) } else{freqs=freqs/1}  })
	
		H = apply(freq, 2, function(freqs){-sum(ifelse(freqs > 0, freqs * log(freqs) , 0))} )
		H = H/log(20)

    }
    if (type == "GroupAA") {
    	if(!all(x %in% AminoAcids)){print("Warning: Data set contains non-Amino Acid elements")}
		counts = apply(x, 2, function(z){table(factor(z, levels = AminoAcids))})
    	
    	grpcounts = apply(counts, 2, function(site){c(site[3]+site[4], site[1]+site[6]+site[8]+site[10]+ site[11]+site[18], site[12]+site[14], site[5]+site[19]+site[20], site[7]+site[9]+site[15], site[2], site[16]+site[17],site[13])})
		rownames(grpcounts)=AAGroups

		freq = apply(grpcounts, 2, function(freqs){if(sum(freqs) > 0){freqs=freqs/sum(freqs) } else{freqs=freqs/1}  })
		counts = grpcounts

		H = apply(freq, 2, function(freqs){-sum(ifelse(freqs > 0, freqs * log(freqs) , 0)) })
		H = H/log(8)

    }
     
    result <- list(counts = counts, freq=freq, H=H)
    return(result)
}


MolecularMI = 
function(x, type, normalized){
 
 	if(!is.matrix(x)){
		if(is.vector(x)){ 
			x = strsplit(x, "")
			x = matrix(unlist(x), nrow = length(x), byrow = TRUE, dimnames = list(names(x)))
		}
		if(is.list(x)){
			x = sapply(x, strsplit, "")
			x = matrix(unlist(x), nrow = length(x), byrow = TRUE, dimnames = list(names(x)))
		}
	}
 
 
   n = ncol(x)	
   MI = matrix(,n,n, dimnames=list(dimnames(x)[[2]], dimnames(x)[[2]])) 
   
        if (type == "DNA")     {base = 4
        						counts = apply(x, 2, function(z){table(factor(z, levels = c("A", "C", "G", "T")))})
        						}
   else if (type == "AA")      {base = 20
           						counts = apply(x, 2, function(z){table(factor(z, levels = AminoAcids))})
           						}   
   else if (type == "GroupAA") {base = 8
           					
								for(k in 1:n){
									for(aa in 1:20){
										x[,k] = replace(x[,k],x[,k]==AminoAcids[aa], AAbyGroup[aa])
									}
								}
								counts = apply(x, 2, function(z){table(factor(z, levels = AAGroups))})
							  }


   Totalcounts = colSums(counts)
   
   for(i in 1:n){
   		for(j in 1:n){
   
 	  		Jointcounts= table(interaction(x[,i], x[,j]))
 	  		if(sum(Jointcounts) > 0){freqs=Jointcounts/sum(Jointcounts) } 
 	  		else					{freqs=Jointcounts}  

 	  		if(Totalcounts[i] > 0){freqs.i=as.matrix(counts[,i]/Totalcounts[i])} 
 	  		else				  {freqs.i=as.matrix(counts[,i])}  
 	  		if(Totalcounts[j] > 0){freqs.j=as.matrix(counts[,j]/Totalcounts[j])} 
 	  		else				  {freqs.j=as.matrix(counts[,j])}  
 
 			H.ij =  -sum(ifelse(freqs > 0, freqs * log(freqs) , 0))
	        H.ij = H.ij/log(base)

 			H.i =  -sum(ifelse(freqs.i > 0, freqs.i * log(freqs.i) , 0))
	        H.i = H.i/log(base)
 			H.j =  -sum(ifelse(freqs.j > 0, freqs.j * log(freqs.j) , 0))
	        H.j = H.j/log(base)

			if (missing(normalized)){	MI[i,j] = H.i + H.j - H.ij					}
 			else{					    MI[i,j] = NMI(H.i, H.j, H.ij, normalized)	}
		}
	}
	return(MI)
}

NMI = 
function(Hx, Hy, Hxy, type=c("NULL","marginal", "joint", "min.marginal", "max.marginal", "min.conditional", "max.conditional")){

	Hy.x = Hxy-Hx
	Hx.y = Hxy-Hy

	NMI.val = switch(type, 
		marginal     	= ifelse( Hx + Hy        > 0,  2*( Hx + Hy - Hxy ) / ( Hx + Hy )     , 0), 
		joint	     	= ifelse( Hxy            > 0,  2*( Hx + Hy - Hxy ) / ( Hxy )         , 0), 
		min.marginal	= ifelse( min(Hx,Hy)     > 0,    ( Hx + Hy - Hxy ) / min(Hx,Hy)      , 0), 
		max.marginal 	= ifelse( max(Hx,Hy)     > 0,    ( Hx + Hy - Hxy ) / max(Hx,Hy)      , 0), 
		min.conditional = ifelse( min(Hx.y,Hy.x) > 0,    ( Hx + Hy - Hxy ) / min(Hx.y,Hy.x)  , 0), 
		max.conditional = ifelse( max(Hx.y,Hy.x) > 0,    ( Hx + Hy - Hxy ) / max(Hx.y,Hy.x)  , 0),
		NULL=,
		 				  Hx + Hy - Hxy)
		 
	return(NMI.val)
}
