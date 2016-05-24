# Calculations of objects needed for figures in the Markarian section

if (FALSE) {
	zeitreihen <- list(Mrk421, Mrk501)
	beschreibung <- c("Mrk421", "Mrk501")
	perioden <- list(1:600, 1:400)
	PPs <- vector("list", 1)
	Crits <- vector("list", 2)
	model <- c("sine", "step")
	for(AA in 1:2){
		zr <- zeitreihen[[AA]]
		beschr <- beschreibung[AA]
		periods <- perioden[[AA]]
		Mod <- model[AA]
		PP_sub <- vector("list", 3)
		Crit_sub <- vector("list", 3)
		for(BB in 1:3) {
			Reg <- c("L2", "tau", "huber")[BB]
			PP_sub[[BB]] <- RobPer(zr, periods=periods, model=Mod, regression=Reg, weighting=FALSE)
			shapes <- betaCvMfit(PP_sub[[BB]])
			Crit_sub[[BB]] <- qbeta(0.95^(1/length(periods)), shape1=shapes[1], shape2=shapes[2])
		}
		Crits[[AA]] <- Crit_sub
		PPs[[AA]] <- PP_sub
	}
	save(PPs, Crits, perioden, file="MrkAnalysis.Rdata")
}