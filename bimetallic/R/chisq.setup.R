chisq.setup <-
function(R, gam_ca, gam_co, ppv, npv, homRR, N_co, maf, prev, model, ...){
setup = list()
#genotypic parameters -- need to restrict rr's s.t. prev is achieved
	if(suppressWarnings(!is.na(as.numeric(model)))){
		rrAa = as.numeric(model)
		model = "independent"
	}
	model = match.arg(model, c("multiplicative", "dominant", "recessive", "independent"))
	rrAA = homRR
	if(model=="multiplicative"){
			 rrAa = sqrt(rrAA) 
			}
		else if(model=="dominant"){
			rrAa = rrAA
			}
		else if(model=="recessive"){
			rrAa = 1
			}

#confusion parameters
cmat = matrix(c(ppv, 1-ppv, 1-npv, npv), nrow=2, byrow=TRUE)

#case/control parameters
numP_CA = ceiling(N_co/R)
numP_CO = N_co
numM_CA = ceiling(numP_CA*gam_ca)
numM_CO = ceiling(numP_CO*gam_co)

AA = maf^2
Aa = 2*maf*(1-maf)
aa = (1-maf)^2
gtFreq = c(AA, Aa, aa)
affaa = prev/(rrAA*AA + rrAa*Aa + aa)
affAa = rrAa*affaa
affAA = rrAA*affaa
affGtFreq = c(affAA, affAa, affaa)
popGtJtFreq = rbind(gtFreq*affGtFreq, gtFreq*(1-affGtFreq))
popGtFreq = rbind(gtFreq*affGtFreq/prev, gtFreq*(1-affGtFreq)/(1-prev))

perGtFreq = cmat%*%popGtFreq
combGtFreq = (c(numM_CA, numM_CO) * perGtFreq + c(numP_CA, numP_CO) * popGtFreq)/c(numM_CA+numP_CA, numM_CO+numP_CO)
combOR = allelicOR(combGtFreq)
trueOR = allelicOR(popGtFreq)
list(popGtFreq=popGtFreq, perGtFreq=perGtFreq, combGtFreq=combGtFreq, numP_CO= numP_CO, numP_CA= numP_CA, numM_CA= numM_CA, numM_CO= numM_CO, combOR = combOR, trueOR=trueOR)
}

