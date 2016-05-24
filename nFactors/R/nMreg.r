nMreg <-
function(x, cor=TRUE, model="components", details=TRUE, ...) {
 x       <- eigenComputes(x, cor=cor, model=model, ...)
 nlength <- 3
 detail  <- NULL
 n       <- length(x)
 if (n < 6) stop("The number of variables must be at least 6.")
 i       <- 1
 mreg    <- tmreg <- tmreg2 <-pmreg <- numeric(n-5)
 while (i <= (length(x)-5)) {
  xa        <- c(1:(i+2))
  ya        <- x[1:(i+2)]
  ma        <- lm(ya ~ xa)
  Syx.a     <- sd(ya)*sqrt((1-summary(ma)$r.squared) * ((length(ya)-1)/(length(ya)-2))) # Howell(2008, p. 253)
  compa     <- ma$coef[2]
  seCompa   <- summary(ma)$coef[2,2]

  xb        <- c((i+1+nlength):length(x))
  yb        <- x[(i+1+nlength):length(x)]
  mb        <- lm(yb ~ xb)
  Syx.b     <- sd(yb)*sqrt((1-summary(mb)$r.squared) * ((length(yb)-1)/(length(yb)-2))) # Howell(2008, p. 253)
  compb     <- mb$coef[2]
  seCompb   <- summary(mb)$coef[2,2]

  mreg[i]   <- compb - compa
  semreg    <- sqrt((Syx.a^2)/((length(xa)-1)*sd(xa)^2) + (Syx.b^2)/((length(xb)-1)*sd(xb)^2))     # Se_dif_b -> Howell(2008, p. 259, 266)
  tmreg[i]  <- (compb - compa)/(semreg)
  tmreg2[i] <- (mreg[i])/sqrt(seCompa^2 + seCompb^2) # Il semble, selon moi, qu'il y aurait une erreur dans la formule de Zoski et Just. Et ce serait la bonne formul, comme celle plu shaut, mais plus rapide de calcul.
  pmreg[i]  <- pt(tmreg[i],(length(xa)-1) + (length(xb)-1) - 4, lower=FALSE, log=TRUE)
  i         <- i + 1
  }
 if (details == TRUE) detail  <- data.frame(v=(1:(n-5)),values=x[1:(n-5)], mreg=mreg, tmreg=tmreg, pmreg=pmreg)
 mreg  <- as.numeric(which(mreg ==max( mreg, na.rm=TRUE)) + nlength)
 tmreg <- as.numeric(which(tmreg==max(tmreg, na.rm=TRUE)))
 pmreg <- as.numeric(which(pmreg==min(pmreg, na.rm=TRUE)))
 res   <- list(detail=detail, nFactors=c(b=mreg,t.p=tmreg,p.b=pmreg))
 class(res) <- c("nFactors","list")
 return(res)
 }