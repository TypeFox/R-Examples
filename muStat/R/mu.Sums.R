`mu.Sums` <-
function(GE, dsgn=1, wght=TRUE)
{
    ICW <- function(GE,dsgn) {
        wgt <- colSums(GE|t(GE))        # 0: orig NA, 1: like NA, >1: no NA
        sqrt(wgt*(wgt>1)/if (dsgn!=1) colSums(dsgn) else nrow(GE))
    }
    GE <- sq.matrix(GE)
    nE <- colSums(GE*t(GE)) 
    nA <- colSums(GE) - nE
    nB <- rowSums(GE) - nE
	weight <- ifelse1(wght, ICW(GE,dsgn), rep(1, length(nE)))
	list (
		score  = (nB-nA) * ifelse(weight==0,NA,1),
		weight = weight,
		nBelow = nB,
		nAbove = nA,
		nEqual = nE
	)   
}

