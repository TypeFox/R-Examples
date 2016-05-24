f.HWE <- function(data, quiet.warning = F){
##
## TEST A SINGLE LOCUS WITH POSSIBLY MULTIPLE ALLELES FOR HWE
## DATA ARE TWO COLUMNS WITH ALLELE 1 AND ALLELE 2, RESP.
## ONLY NUMERIC DATA, CODED CONSECUTIVELY 1, 2, 3, ...
##
## NOTE: NECESSARY TO MAKE SURE THAT ALL POSSIBLE GENOTYPES ARE COUNTED
#
if(!is.numeric(data) | dim(data)[2] != 2) stop("Wrong data type for HWE testing!", call. = F)
if(nrow(data) == 0){# CAN HAPPEN IF E.G. SELECT comb.sex = "boys" UNDER xchrom = T
	.ut <- list(table = NA, freq = NA, warnings = "No data available for HWE testing", chisq = NA, df = NA, p.value = NA, n.miss.geno = NA, fail = FALSE)
	class(.ut) <- "HWE.test"
	return(.ut)
}

## COUNT AND REMOVE MISSING
.d1 <- dim(data)[1]
.data <- na.omit(data)
.d2 <- dim(.data)[1]
.n.miss <- .d1 - .d2
#
## TOTAL ALLELE COUNT
.a <- c(.data[,1], .data[,2])
.aux <- sort(unique(.a))
if(any(.aux != seq(along = .aux))) stop("Wrong data type for HWE testing!", call. = F)
.frek <- as.numeric(tabulate(.a))
.nall <- length(.frek)
#
## FREQUENCY COUNT
.data <- f.aggregate(.data)
## SET FREQS IN ARRAY
.tab <- matrix(0, ncol = length(.aux), nrow = length(.aux))
.tab[as.matrix(.data[,1:2])] <- .data$freq
## FLIP LOWER TRI TO ADD GENOTYPES
.lower.tri <- .tab * (lower.tri(.tab))
.tab <- .tab + t(.lower.tri)
.tab[lower.tri(.tab)] <- NA
## SET IN LONG FORMAT
.A <- expand.grid(a1 = .aux, a2 = .aux)
.grid <- dframe(.A, freq = as.numeric(.tab))
.grid <- .grid[!lower.tri(.tab),] # REMOVES REDUNDANT GENOTYPES
#
## COMPUTE EXPECTED FREQS
.exp <- .frek[.grid$a1] * .frek[.grid$a2] * (1 + (.grid$a1 != .grid$a2)) / (2*sum(.frek))
.grid$exp <- .exp
if(abs(sum(.grid$freq) - sum(.grid$exp)) > 1e-5) stop("Problem with HWE test")
#
## COMPUTE CHI-SQUARED
.warn <- "None"
.less5 <- any(.grid$exp < 5)

if(.less5){
	.warn <- "Expected frequency less than 5"
	if(!quiet.warning) warning("Expected frequency less than 5 in HWE test!")
}
.grid$chisq <- (.grid$freq - .grid$exp)^2/.grid$exp
.chisq <- sum(.grid$chisq)
.df <- .nall * (.nall - 1)/2
.chisq.pvalue <- 1 - pchisq(.chisq, df = .df)
#
## PREPARE OUTPUT
.ut <- list(table = .grid, freq = .frek, warnings = .warn, chisq = .chisq, df = .df, p.value = .chisq.pvalue, n.miss.geno = .n.miss, fail = FALSE)
class(.ut) <- "HWE.test"


return(.ut)
}
