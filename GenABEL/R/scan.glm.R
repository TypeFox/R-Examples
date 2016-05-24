"scan.glm" <- 
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

   fla2 <- as.formula(sub("CRSNP","as.factor(mygt)",formula,ignore.case=TRUE))
   fla1 <- as.formula(sub("CRSNP","mygt",formula,ignore.case=TRUE))
   fla0 <- as.formula(sub("CRSNP","DuMmY",formula,ignore.case=TRUE))

  P1df <- rep(1,length(snpsubset))
  P2df <- rep(1,length(snpsubset))
  effB <- rep(0,length(snpsubset))
  effAB <- rep(0,length(snpsubset))
  effBB <- rep(0,length(snpsubset))
#  print(ch)

  for (w in 1:(chunks+1)) {
  mgtd <- as.double(gtdata[,ch[[w]]])
  for (i in 1:length(ch[[w]])) {
    mygt <- mgtd[,i]
    polym <- length(levels(as.factor(mygt)))
    if (polym<=1) {
       cat("Marker",snpsubset[(w-1)*chsize+i],"is monomorphic; skipping in analysis\n")
       P1df[(w-1)*chsize+i]=1.0
       P2df[(w-1)*chsize+i]=1.0
    }  else {
    DuMmY <- rep(0,length(mygt))
    DuMmY <- replace(DuMmY,is.na(mygt),NA)
    if (family$family != "gaussian") {
      m1  <- glm(fla1,family = family,data=phdata)
      if (polym>2) {
        m2  <- glm(fla2,family = family,data=phdata)
      } else {
	m2 <- m1
      }
      m0  <- glm(fla0,family = family,data=phdata)
      anv1 <- anova(m0,m1,test="Chisq")
      anv2 <- anova(m0,m2,test="Chisq")
      P1df[(w-1)*chsize+i] <- anv1[2, grep("^P.*Chi",names(anv1))]
      P2df[(w-1)*chsize+i] <- anv2[2, grep("^P.*Chi",names(anv2))]
      if (polym>2) {
      	effB[(w-1)*chsize+i] <- m1$coeff["mygt"]
      	effAB[(w-1)*chsize+i] <- m1$coeff["as.factor(mygt)1"]
      	effBB[(w-1)*chsize+i] <- m1$coeff["as.factor(mygt)2"]
      } else {
      	effB[(w-1)*chsize+i] <- m1$coeff["mygt"]
	effAB[(w-1)*chsize+i] <- effB[(w-1)*chsize+i]
	effBB[(w-1)*chsize+i] <- effB[(w-1)*chsize+i]
      }
      if (is.na(P1df[(w-1)*chsize+i])) P1df[(w-1)*chsize+i] = 1.0
      if (is.na(P2df[(w-1)*chsize+i])) P2df[(w-1)*chsize+i] = 1.0
    } else {
      m1  <- lm(fla1,data=phdata)
      if (polym>2) {
        m2  <- lm(fla2,data=phdata)
      } else {
	m2 <- m1
      }
      m0  <- lm(fla0,data=phdata)
      anv1 <- anova(m0,m1,test="Chisq")
      anv2 <- anova(m0,m2,test="Chisq")
	  P1df[(w-1)*chsize+i] <- anv1[2, grep("^P.*Chi",names(anv1))]
	  P2df[(w-1)*chsize+i] <- anv2[2, grep("^P.*Chi",names(anv2))]
      if (polym>2) {
      	effB[(w-1)*chsize+i] <- m1$coeff["mygt"]
      	effAB[(w-1)*chsize+i] <- m1$coeff["as.factor(mygt)1"]
      	effBB[(w-1)*chsize+i] <- m1$coeff["as.factor(mygt)2"]
      } else {
      	effB[(w-1)*chsize+i] <- m1$coeff["mygt"]
	effAB[(w-1)*chsize+i] <- effB[(w-1)*chsize+i]
	effBB[(w-1)*chsize+i] <- effB[(w-1)*chsize+i]
      }
      if (is.na(P1df[(w-1)*chsize+i])) P1df[(w-1)*chsize+i] = 1.0
      if (is.na(P2df[(w-1)*chsize+i])) P2df[(w-1)*chsize+i] = 1.0
    } 
    }
    if (bcast && round(((w-1)*chsize+i)/bcast) == ((w-1)*chsize+i)/bcast) {
		cat("\b\b\b\b\b\b\b\b",round(100*(((w-1)*chsize+i)/gtdata@nsnps),digits=2),"%",sep="");
		flush.console();
	}
  }
  }
  if (bcast && gtdata@nsnps >= bcast) cat("\n")

  map <- gtdata@map
  chromosome <- gtdata@chromosome
  med1df <- median(qchisq(1.-P1df,df=1))
  med2df <- median(qchisq(1.-P2df,df=2))
  out <- list(P1df = P1df, P2df=P2df, medChi1df = med1df, medChi2df = med2df, snpnames = snpsubset, idnames = idsubset, formula = match.call(), family = family, map = map, chromosome = chromosome, effB=effB, effAB=effAB,effBB=effBB)
  out$Pc1df <- rep(NA,length(P1df))
  class(out) <- "scan.gwaa"
  out
}
