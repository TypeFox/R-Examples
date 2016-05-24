f.prob <- function(haplo.freq, RR, RRcm, RRcf, RRstar, maternal, poo, RR.mat, RRstar.mat, n.cases, n.sim, xchrom, BR.girls, .grid = .grid){
##
## NALL GIVE THE NUMBER OF ALLELES IN EACH LOCI
## RR IS A VECTOR OF RELATIVE RISKS 
## RRSTAR IS A VECTOR OF R* IN MODEL
## MATERNAL, ARGUMENT TO INDICATE IF THERE ARE MATERNAL EFFECTS
## RR.MAT AND RRSTAR.MAT MATERNAL RR AND RRSTAR
## N.CASES : THE NUMBER OF CASES
## N.SIM : THE NUMBER OF SIMULATIONS
#
if(xchrom) .sex <- .grid[, "sex"]
#
## HAPLOTYPE FREQUENCIES
.hf1 <- haplo.freq[.grid[,1]]
.hf2 <- haplo.freq[.grid[,2]]
.hf3 <- haplo.freq[.grid[,3]]
if(!xchrom) .hf4 <- haplo.freq[.grid[,4]]
#
## RELATIVE RISKS, SINGLE DOSE
if(!poo){
	.RR1 <- RR[.grid[,"h2.m"]]
	.RR2 <- RR[.grid[,"h2.f"]]
}else{
	.RR1 <- RRcm[.grid[,"h2.m"]]
	.RR2 <- RRcf[.grid[,"h2.f"]]	
}
if(maternal){
	.RR1.mat <- RR.mat[.grid[,"h1.m"]]
	.RR2.mat <- RR.mat[.grid[,"h2.m"]]
}
#
## RELATIVE RISKS, DOUBLE DOSE
.RRstar <- rep(1, length(.RR1))
.ind.double <- .grid[,"h2.m"] == .grid[,"h2.f"]
.RRstar[.ind.double] <- RRstar[.grid[,"h2.m"]][.ind.double]
if(maternal){
	.RRstar.mat <- rep(1, length(.RR1.mat))
	.ind.double.mat <- .grid[,"h1.m"] == .grid[,"h2.m"]
	.RRstar.mat[.ind.double.mat] <- RRstar.mat[.grid[,"h2.m"]][.ind.double.mat]
}
#
## MULTINOMIAL PROBABILITIES, ASSUMING HWE AND MULTIPLICATIVE RISKS
if(xchrom){
	.prob <- .hf1 * .hf2 * .hf3 * .RR1 * ((.sex == 1) + BR.girls * .RR2 * (.sex == 2))
	cat("bare multiplikativ (for now...)\n")
	if(maternal){
		.prob <- .prob * .RR1.mat * .RR2.mat
	}
}else{
	.prob <- .hf1 * .hf2 * .hf3 * .hf4 * .RR1 * .RR2 * .RRstar
	if(maternal){
		.prob <- .prob * .RR1.mat * .RR2.mat * .RRstar.mat
	}
}
#
.prob <- .prob/sum(.prob)
#
return(.prob)
#
}
