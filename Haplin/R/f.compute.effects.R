"f.compute.effects"<-
function(x, maternal, reference.method, ref.cat, n.all, info)
{
# COMPUTES ALLELE FREQUENCIES, RELATIVE RISKS ETC. USING THE RESULTS
# FROM f.tri.glm. x IS THE VECTOR OF ESTIMATED COEFFICIENTS, OR POSSIBLY
# A MATRIX WHERE EACH ROW IS A SET OF ESTIMATED COEFFICIENTS
#
#### PREPARE: ##########################
design <- info$model$design
.poo <- info$model$poo
.fn <- function(x) paste(x, 1:n.all, sep = "")
.safe <- T # USES A MORE STABLE (BUT MORE TIME-CONSUMING) COMPUTATION.
		# SHOULD NOT BE NECESSARY IF NONE OF THE ALLELES ARE RARE
#
## MAKE SURE x IS A MATRIX
if(!is.matrix(x)) {
	.resmat <- matrix(x, nrow = 1)
	dimnames(.resmat) <- list(NULL, names(x))
}
else .resmat <- x #
#
## REMOVE CASE-CONTROL AND SEX COLUMN BEFORE COMPUTATIONS
.resmat <- .resmat[, !is.element(colnames(.resmat), c("sex", "cc")), drop = F]
#
## FILL IN ZEROS FOR NON-ESTIMATED EFFECTS
.resmat <- f.fill.effects(resmat = .resmat, info = info)
#
#
.J <- matrix(1, nrow = n.all, ncol = 1) #
#
#### COMPUTE ALLELE FREQUENCIES: #############
.p <- exp(.resmat[, .fn("mf"), drop = F])
.p <- .p/as.numeric(.p %*% .J)	#
#
#### COMPUTE RELATIVE RISKS FOR CHILD: #############
if(.poo){
	## ONLY REFCAT AND POPULATION, THIS FAR
	.Rcm <- exp(.resmat[, .fn("cm"), drop = F])	#
	.Rcf <- exp(.resmat[, .fn("cf"), drop = F])	#
	.Rstar <- exp(.resmat[, .fn("cdd"), drop = F])	#
	## IN REFCAT PARAMETRIZATION:
	.Rstartilde <- .Rcm * .Rcf * .Rstar
	## IN POPULATION REFERENCE:
	if(reference.method == "population") {
		.Bcm <- as.numeric((.Rcm * .p) %*% .J)	# POPULATION BASELINE
		.Bcf <- as.numeric((.Rcf * .p) %*% .J)	# POPULATION BASELINE
		.Rcmtilde <- .Rcm/.Bcm
		.Rcftilde <- .Rcf/.Bcf
		.Rstartilde <- .Rstartilde/(.Bcm * .Bcf)
	}
}else{
	.R <- exp(.resmat[, .fn("c"), drop = F])	#
	.Rstar <- exp(.resmat[, .fn("cdd"), drop = F])	#
	## IN REFCAT PARAMETRIZATION:
	.Rstartilde <- .R^2 * .Rstar
	## IN POPULATION REFERENCE:
	if(reference.method == "population") {
		.B <- as.numeric((.R * .p) %*% .J)	# POPULATION BASELINE
		.Rtilde <- .R/.B
		.Rstartilde <- .Rstartilde/.B^2 # KORRIGERT FRA .B
	}
	## IN RECIPROCAL REFERENCE:
	if(reference.method == "reciprocal"){ #
	#
		.f.minus <- function(x) as.numeric(x %*% rep(1, dim(x)[2])) - x #
	#
		.f.kryssprod <- function(x, f.minus){
			(x * f.minus(x)) %*% rep(1, dim(x)[2])		
		} #
	#
		.f.kryssprod.minus <- function(x){
			.d2 <- dim(x)[2]
			(as.numeric(x %*% rep(1, .d2)) - x)^2 - (as.numeric(x^2 %*% rep(1, .d2)) - x^2)
		} #
	#
		.f.minus.safe <- function(x){
			.d <- dim(x)
			.ut <- x
			.ut[,] <- NA
			for(i in seq(length = .d[2])){
				.ut[,i] <- x[, -i, drop = F] %*% rep(1, .d[2] -1)
			}
			.ut
		} #
	#		
		.f.kryssprod.minus.safe <- function(x, f.minus, f.kryssprod){
			.d <- dim(x)
			.ut <- x
			.ut[,] <- NA
			for(i in seq(length = .d[2])){
				.ut[,i] <- f.kryssprod(x[,-i, drop = F], f.minus)
			}
			.ut
		} #
	#
	#
		.pR <- .p * .R #
	#
		if(.safe){
			.now <- proc.time()[1]
			.Pi.safe <- .R * .f.minus.safe(.pR) / .f.minus.safe(.p)
			.Pi.minus.safe <- .f.kryssprod.minus.safe(.pR, .f.minus.safe, .f.kryssprod) / .f.kryssprod.minus.safe(.p, .f.minus.safe, .f.kryssprod)
			.Fi.safe <- .Pi.safe/.Pi.minus.safe
			.Fi.startilde.safe <- .Rstartilde/.Pi.minus.safe
			f.vis(paste("used", proc.time()[1] - .now, "\n"), vis = F)
			.Fi <- .Fi.safe
			.Fi.startilde <- .Fi.startilde.safe
		}# END if(.safe)
		if(!.safe){
			.now <- proc.time()[1]
			.Pi <- .R * .f.minus(.pR) / .f.minus(.p)
	##		.Pi.minus <- (.f.minus(.pR)^2 - .f.minus(.pR^2))/(.f.minus(.p)^2 - .f.minus(.p^2))
			.Pi.minus <- .f.kryssprod.minus(.pR)/.f.kryssprod.minus(.p)
			if(any(.Pi.minus == 0)){
				.Pi.minus[.Pi.minus == 0] <- 1e-10  
				warning("NAs generated during simulation, perhaps due to rare haplotype!")
			}
			.Fi <- .Pi/.Pi.minus
			.Fi.startilde <- .Rstartilde/.Pi.minus
			f.vis(paste("used", proc.time()[1] - .now, "\n"), vis = F)
		}# END if(!.safe)
	}# END reciprocal
}# END if(!.poo)
#
#
if(maternal) {#
#### COMPUTE RELATIVE RISKS FOR MATERNAL EFFECTS, IF APPLICABLE: #############
	.Rm <- exp(.resmat[, .fn("m"), drop = F])	#
	.Rmstar <- exp(.resmat[, .fn("mdd"), drop = F]) # 
	## IN REFCAT PARAMETRIZATION:
		.Rmstartilde <- .Rm^2 * .Rmstar
	## IN POPULATION REFERENCE:
	if(reference.method == "population") {
		.Bm <- as.numeric((.Rm * .p) %*% .J)	# POPULATION BASELINE
		.Rmtilde <- .Rm/.Bm
		.Rmstartilde <- .Rmstartilde/.Bm^2 # KORRIGERT FRA .Bm
	}
	## IN RECIPROCAL REFERENCE:
	if(reference.method == "reciprocal"){
		.pRm <- .p * .Rm #
		#
		if(.safe){
			.now <- proc.time()[1]
			.Pim.safe <- .Rm * .f.minus.safe(.pRm) / .f.minus.safe(.p)
			.Pim.minus.safe <- .f.kryssprod.minus.safe(.pRm, .f.minus.safe, .f.kryssprod) / .f.kryssprod.minus.safe(.p, .f.minus.safe, .f.kryssprod)
			.Fim.safe <- .Pim.safe/.Pim.minus.safe
			.Fim.startilde.safe <- .Rmstartilde/.Pim.minus.safe
			f.vis(paste("used", proc.time()[1] - .now, "\n"), vis = F)
			.Fim <- .Fim.safe
			.Fim.startilde <- .Fim.startilde.safe
		}# END SAFE
		if(!.safe){
			.now <- proc.time()[1]
			.Pim <- .Rm * .f.minus(.pRm) / .f.minus(.p)
			.Pim.minus <- .f.kryssprod.minus(.pRm)/.f.kryssprod.minus(.p)
			if(any(.Pim.minus == 0)){
				.Pim.minus[.Pim.minus == 0] <- 1e-10  
				warning("NAs generated during simulation, perhaps due to rare haplotype!")
			}
			.Fim <- .Pim/.Pim.minus
			.Fim.startilde <- .Rmstartilde/.Pim.minus
			f.vis(paste("used", proc.time()[1] - .now, "\n"), vis = F)
		}# END if(!.safe)
	}# END RECIPROCAL
}# END if(maternal)
#
#
#
#### OUTPUT: ##########################
if(.poo){
	if(maternal){
		if(reference.method == "ref.cat") {
			.ut <- cbind(.p, .Rcm, .Rcf, .Rcm/.Rcf, .Rstartilde, .Rm, .Rmstartilde)
		}
		if(reference.method == "population") {
			.ut <- cbind(.p, .Rcmtilde, .Rcftilde, .Rcmtilde/.Rcftilde, .Rstartilde, .Rmtilde, .Rmstartilde)	#
		}	
		colnames(.ut) <- c(.fn("p"), .fn("RRcm"), .fn("RRcf"), .fn("RRcm_RRcf"), .fn("RRcdd"), .fn("RRm"), .fn("RRmdd"))
	}else{
		if(reference.method == "ref.cat") {
			.ut <- cbind(.p, .Rcm, .Rcf, .Rcm/.Rcf, .Rstartilde)
		}
		if(reference.method == "population") {
			.ut <- cbind(.p, .Rcmtilde, .Rcftilde, .Rcmtilde/.Rcftilde, .Rstartilde)
		}	
		colnames(.ut) <- c(.fn("p"), .fn("RRcm"), .fn("RRcf"), .fn("RRcm_RRcf"), .fn("RRcdd"))
	}# END if(!maternal)
}else{
	if(maternal){
		if(reference.method == "ref.cat") {
			.ut <- cbind(.p, .R, .Rstartilde, .Rm, .Rmstartilde)
		}
		if(reference.method == "population") {
			.ut <- cbind(.p, .Rtilde, .Rstartilde, .Rmtilde, .Rmstartilde)	#
		}	
		if(reference.method == "reciprocal") {
			.ut <- cbind(.p, .Fi, .Fi.startilde, .Fim, .Fim.startilde)	#
		}
		colnames(.ut) <- c(.fn("p"), .fn("RRc"), .fn("RRcdd"), .fn("RRm"), .fn("RRmdd"))
	}else{
		if(reference.method == "ref.cat") {
			.ut <- cbind(.p, .R, .Rstartilde)
		}
		if(reference.method == "population") {
			.ut <- cbind(.p, .Rtilde, .Rstartilde)	#
		}
		if(reference.method == "reciprocal") {
			.ut <- cbind(.p, .Fi, .Fi.startilde)	#
		}
		colnames(.ut) <- c(.fn("p"), .fn("RRc"), .fn("RRcdd"))
	}# END if(!maternal)
}# END if(!.poo)
#
return(.ut)
}
