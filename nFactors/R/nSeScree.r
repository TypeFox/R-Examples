nSeScree <-
function(x, cor=TRUE, model="components", details=TRUE, r2limen=0.75, ...) {
 x               <- eigenComputes(x, cor=cor, model=model, ...)
 detail          <- NULL
 n               <- length(x)
 criteria        <- 1/n
 seScreeCriteria <- R2Criteria <- 0
 if (n < 3) stop("The number of variables must be at least 3.")
 i               <- 1
 seScree         <- R2 <- numeric(n-3)
 while ((i) <= (n-2)) {
  xa              <- c(i:n)
  ya              <- x[i:n]
  ma              <- lm(ya ~ xa)
  seScree[i]      <- sd(ya)*sqrt((1-summary(ma)$r.squared) * ((length(ya)-1)/(length(ya)-2))) # Howell(2008, p. 253)
  seScreeCriteria <- seScreeCriteria + as.numeric(seScree[i] > criteria)
  R2[i]           <- summary(ma)$r.squared
  R2Criteria      <- R2Criteria + as.numeric(R2[i] < r2limen)
  i               <- i + 1
  }
 if (details == TRUE) detail  <- data.frame(v=(1:(n-2)),values=x[1:(n-2)], seScree, R2)
 seScree <- seScreeCriteria
 R2      <- R2Criteria
 res     <- list(detail=detail, nFactors=c(se=seScree, R2=R2))
 class(res) <- c("nFactors","list")
 return(res)
 }