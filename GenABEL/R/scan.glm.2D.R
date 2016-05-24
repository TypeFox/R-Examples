"scan.glm.2D" <- 
		function (formula, family = gaussian(), data, snpsubset, idsubset, bcast=50) {
	if (!is(data,"gwaa.data")) {
		stop("wrong data class: should be gwaa.data")
	}
	if (!is.character(formula)) stop("formula must be character object (apply \"s)")
	if (is.character(family))
		family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family))
		family <- family()
	if (is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
	}
	if (grep("CRSNP",formula,ignore.case=TRUE)!=1) stop("formula must contain CRSNP variable to be replaced with the analysis SNPs")
	if (missing(snpsubset)) snpsubset <- data@gtdata@snpnames
	if (missing(idsubset)) idsubset <- data@gtdata@idnames
	if (is.logical(snpsubset) || is.numeric(snpsubset)) snpsubset <- data@gtdata@snpnames[snpsubset]
	gtdata <- data[idsubset,snpsubset]@gtdata
	phdata <- data[idsubset,snpsubset]@phdata
	allnams <- gtdata@snpnames
	chsize <- ceiling(length(allnams)/80)
	ch <- list()
	if(chsize != 1) {
		chunks <- ceiling(length(allnams)/chsize)-1
		for (i in 1:chunks) {
			ch[[i]] = allnams[((i-1)*chsize+1):(i*chsize)]
		}
		ch[[chunks+1]] = allnams[(chunks*chsize+1):length(allnams)]
	} else {ch[[1]] <- allnams;chunks=0}
	
	fla2 <- as.formula(sub("CRSNP","as.factor(mygt1)*as.factor(mygt2)",formula,ignore.case=TRUE))
	fla2.i0 <- as.formula(sub("CRSNP","as.factor(mygt1)+as.factor(mygt2)",formula,ignore.case=TRUE))
	fla2.1 <- as.formula(sub("CRSNP","as.factor(mygt1)*mygt2",formula,ignore.case=TRUE))
	fla2.1.i0 <- as.formula(sub("CRSNP","as.factor(mygt1)+mygt2",formula,ignore.case=TRUE))
	fla2.2 <- as.formula(sub("CRSNP","mygt1*as.factor(mygt2)",formula,ignore.case=TRUE))
	fla2.2.i0 <- as.formula(sub("CRSNP","mygt1+as.factor(mygt2)",formula,ignore.case=TRUE))
	fla1 <- as.formula(sub("CRSNP","mygt1*mygt2",formula,ignore.case=TRUE))
	fla1.i0 <- as.formula(sub("CRSNP","mygt1+mygt2",formula,ignore.case=TRUE))
	fla0 <- as.formula(sub("CRSNP","DuMmY1*DuMmY2",formula,ignore.case=TRUE))
	
	nsnps <- gtdata@nsnps
	P1df <- matrix(rep(NA,(nsnps*nsnps)),nrow=nsnps)
	P2df <- P1df
	Pint1df <- P1df
	Pint2df <- P1df
#  print(ch)
	
	donan<-0
	for (i1 in 1:(length(snpsubset)-1)) {
		for (i2 in (i1+1):length(snpsubset)) {
			mygt1 <- as.numeric(gtdata[,snpsubset[i1]])
			mygt2 <- as.numeric(gtdata[,snpsubset[i2]])
			polym1 <- length(levels(as.factor(mygt1)))
			polym2 <- length(levels(as.factor(mygt2)))
			if (polym1<=1 || polym2<=1) {
				cat("One of markers",snpsubset[i1],snpsubset[i2],"is (are) monomorphic; skipping in analysis\n")
				P1df[i1,i2]=1.0
				P2df[i1,i2]=1.0
				Pint1df[i1,i2]=1.0
				Pint2df[i1,i2]=1.0
			}  else {
				DuMmY1 <- rep(0,length(mygt1))
				DuMmY1 <- replace(DuMmY1,is.na(mygt1),NA)
				DuMmY2 <- rep(0,length(mygt2))
				DuMmY2 <- replace(DuMmY2,is.na(mygt2),NA)
				if (family$family != "gaussian") {
					m1  <- glm(fla1,family = family,data=phdata)
					m1.i0  <- glm(fla1.i0,family = family,data=phdata)
					if (polym1>2 && polym2>2) {
						m2  <- glm(fla2,family = family,data=phdata)
						m2.i0  <- glm(fla2.i0,family = family,data=phdata)
					} else if (polym1>2) {
						m2  <- glm(fla2.1,family = family,data=phdata)
						m2.i0  <- glm(fla2.1.i0,family = family,data=phdata)
					} else if (polym2>2) {
						m2  <- glm(fla2.2,family = family,data=phdata)
						m2.i0  <- glm(fla2.2.i0,family = family,data=phdata)
					} else {
						m2 <- m1
						m2.i0 <- m1.i0
					}
					m0  <- glm(fla0,family = family,data=phdata)
					anv1 <- anova(m0,m1,test="Chisq")
					anv2 <- anova(m0,m2,test="Chisq")
					anv1.i <- anova(m1.i0,m1,test="Chisq")
					anv2.i <- anova(m2.i0,m2,test="Chisq")
					P1df[i1,i2] <- anv1[2, grep("^P.*Chi",names(anv1))]
					P2df[i1,i2] <- anv2[2, grep("^P.*Chi",names(anv2))]
					Pint1df[i1,i2] <- anv1.i[2, grep("^P.*Chi",names(anv1.i))]
					Pint2df[i1,i2] <- anv2.i[2, grep("^P.*Chi",names(anv2.i))]
					if (is.na(P1df[i1,i2])) P1df[i1,i2] = 1.0
					if (is.na(P2df[i1,i2])) P2df[i1,i2] = 1.0
					if (is.na(Pint1df[i1,i2])) Pint1df[i1,i2] = 1.0
					if (is.na(Pint2df[i1,i2])) Pint2df[i1,i2] = 1.0
				} else {
					m1  <- lm(fla1,data=phdata)
					m1.i0  <- lm(fla1.i0,data=phdata)
					if (polym1>2 && polym2>2) {
						m2  <- lm(fla2,data=phdata)
						m2.i0  <- lm(fla2.i0,data=phdata)
					} else if (polym1>2) {
						m2  <- lm(fla2.1,data=phdata)
						m2.i0  <- lm(fla2.1.i0,data=phdata)
					} else if (polym2>2) {
						m2  <- lm(fla2.2,data=phdata)
						m2.i0  <- lm(fla2.2.i0,data=phdata)
					} else {
						m2 <- m1
						m2.i0 <- m1.i0
					}
					m0  <- lm(fla0,data=phdata)
					anv1 <- anova(m0,m1,test="Chisq")
					anv2 <- anova(m0,m2,test="Chisq")
					anv1.i <- anova(m1.i0,m1,test="Chisq")
					anv2.i <- anova(m2.i0,m2,test="Chisq")
					P1df[i1,i2] <- anv1[2, grep("^P.*Chi",names(anv1))]
					P2df[i1,i2] <- anv2[2, grep("^P.*Chi",names(anv2))]
					Pint1df[i1,i2] <- anv1.i[2, grep("^P.*Chi",names(anv1.i))]
					Pint2df[i1,i2] <- anv2.i[2, grep("^P.*Chi",names(anv2.i))]
					if (is.na(P1df[i1,i2])) P1df[i1,i2] = 1.0
					if (is.na(P2df[i1,i2])) P2df[i1,i2] = 1.0
					if (is.na(Pint1df[i1,i2])) Pint1df[i1,i2] = 1.0
					if (is.na(Pint2df[i1,i2])) Pint2df[i1,i2] = 1.0
				} 
			}
			donan<-donan+1
			if (bcast && round((donan)/bcast) == (donan)/bcast) {
				cat("\b\b\b\b\b\b\b\b",round(100*donan/((nsnps-1)*nsnps/2),digits=2),"%",sep="");
				flush.console();
			}
		}
	}
	if (bcast && donan>=bcast) cat("\n")
	
	map <- gtdata@map
	chromosome <- gtdata@chromosome
	med1df <- median(qchisq(1.-P1df,df=1))
	med2df <- median(qchisq(1.-P2df,df=2))
	colnames(P1df) <- snpsubset
	rownames(P1df) <- snpsubset #[length(snpsubset):1]
	out <- list(P1df = P1df, Pint1df=Pint1df, P2df=P2df, Pint2df=Pint2df, medChi1df = med1df, medChi2df = med2df, snpnames = snpsubset, idnames = idsubset, formula = match.call(), family = family, map = map, chromosome = chromosome)
	out$Pc1df <- rep(NA,length(P1df))
	class(out) <- "scan.gwaa.2D"
	out
}
