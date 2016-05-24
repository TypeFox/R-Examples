
calcATHBpmv <- function(trm, psych, ta, tr, vel, rh, met, wme){

	### calc with clo according to fix effects

	HFtomet <- .092
	PContoHF <- (-.549)
	tointopCon <- (-.06264)
	metadaptPhys <- .193 * HFtomet
	metadaptPCon <- PContoHF * tointopCon * HFtomet

	# adjustment of clothing value based on behavioural adaptation
	cloeq <- 1.252594 + trm * (-0.03023063)
	cloeq <- ifelse (cloeq > 1, 1, ifelse (cloeq < .46, .46, cloeq))

	# adjustment of met based on physiological adaptation
	dmetadaptPhys <- ifelse (((trm - 18) * metadaptPhys) < 0, 0, ((trm - 18) * metadaptPhys))

	# adjustment of met based on psychological adaptation ( variable part)
	dmetadaptpsychVar <- ifelse (((ta - 20) * metadaptPCon) < 0, 0, ((ta - 20) * metadaptPCon))

	# adjustment of met based on psychological adaptation (fixed part)
	dmetadaptpsychFix <- psych * PContoHF * HFtomet

	metadapt <- met - dmetadaptPhys + dmetadaptpsychVar + dmetadaptpsychFix

	comfortData <- data.frame(calcPMVPPD(ta, tr, vel, rh, cloeq, metadapt, wme))
	giveDat <- with(comfortData, get("pmv"))
	giveDat
}

calcATHBset <- function(trm, psych, ta, tr, vel, rh, met, wme, pb, ltime, ht, wt){

	### calc with clo according to fix effects

	HFtomet <- .092
	PContoHF <- (-.549)
	tointopCon <- (-.06264)
	metadaptPhys <- .193 * HFtomet
	metadaptPCon <- PContoHF * tointopCon * HFtomet

	# adjustment of clothing value based on behavioural adaptation
	cloeq <- 1.252594 + trm * (-0.03023063)
	cloeq <- ifelse (cloeq > 1, 1, ifelse (cloeq < .46, .46, cloeq))

	# adjustment of met based on physiological adaptation
	dmetadaptPhys <- ifelse (((trm - 18) * metadaptPhys) < 0, 0, ((trm - 18) * metadaptPhys))

	# adjustment of met based on psychological adaptation ( variable part)
	dmetadaptpsychVar <- ifelse (((ta - 20) * metadaptPCon) < 0, 0, ((ta - 20) * metadaptPCon))

	# adjustment of met based on psychological adaptation (fixed part)
	dmetadaptpsychFix <- psych * PContoHF * HFtomet

	metadapt <- met - dmetadaptPhys + dmetadaptpsychVar + dmetadaptpsychFix

	comfortData <- data.frame(calc2Node(ta, tr, vel, rh, cloeq, metadapt, wme, pb, ltime, ht, wt, obj = "set"))
	giveDat <- with(comfortData, get("set"))
	giveDat
}

calcATHBpts <- function(trm, psych, ta, tr, vel, rh, met, wme, pb, ltime, ht, wt){
	
	set <- calcATHBset(trm, psych, ta, tr, vel, rh, met, wme, pb, ltime, ht, wt)
	ATHBpts <- .25 * set - 6.03
	#names(ATHBpts) <- "ATHBpts"
	data.frame(ATHBpts)
	
}

