# Calculations of objects needed for figures in the GROJ0422+32 section

if (FALSE) {
	data(star_groj0422.32)
	zr <- star_groj0422.32
	perioden_groj <- 1:330
	PP_groj  <- vector("list", 3)
	Crit_groj<- c(0,0,0)
	for(BB in 1:3) {
		Reg <- c("L2", "tau", "huber")[BB]
		Mod <-  "sine"
		PP_groj[[BB]] <- RobPer(zr, periods=perioden_groj, model=Mod, regression=Reg, weighting=FALSE)
		shapes <- betaCvMfit(PP_groj[[BB]])
		Crit_groj[BB] <-qbeta(0.95^(1/length(perioden_groj)), shape1=shapes[1], shape2=shapes[2])
	}
	
	# Data and analysis with a sine added
	grojfluc <- star_groj0422.32
	grojfluc[,2] <- star_groj0422.32[,2]+0.005*sin(star_groj0422.32[,1]/30*2*pi)
	PP_grojfluc <- vector("list", 3)
	Crit_grojfluc <- c(0,0,0)
	for(BB in 1:3) {
		Reg <- c("L2", "tau", "huber")[BB]
		Mod <-  "sine"
		PP_grojfluc[[BB]] <- RobPer(grojfluc, periods=perioden_groj, model=Mod, regression=Reg, weighting=FALSE)
		shapes <- betaCvMfit(PP_grojfluc[[BB]])
		Crit_grojfluc[BB] <- qbeta(0.95^(1/length(perioden_groj)), shape1=shapes[1], shape2=shapes[2])
	}
	save(star_groj0422.32, PP_groj, Crit_groj, grojfluc, PP_grojfluc, Crit_grojfluc, file="grojanalysis.Rdata")
}