"lmrdiscord" <-
function(site=NULL, t2=NULL, t3=NULL, t4=NULL,
         Dcrit=NULL, digits=4, lmrdigits=4, sort=TRUE,
         alpha1=0.10, alpha2=0.01, ...) {

  if(alpha1 < alpha2) {
    warning("Inverted significance levels, flipping for you")
    alpha <- alpha1
    alpha1 <- alpha2
    alpha2 <- alpha
  }

  if(is.null(site)) site <- 1:length(t2)
  nsite <- length(site)

  if(nsite < 5) {
     warning("Too few sites (<5), although >10 is preferred")
     return(NULL)
  }

  if(length(t2) != nsite) stop("Tau2 length is inconsistent with site length")
  if(length(t3) != nsite) stop("Tau3 length is inconsistent with site length")
  if(length(t4) != nsite) stop("Tau4 length is inconsistent with site length")

  tabDc <- list(Dc = c(1.3333, 1.6481, 1.9166, 2.1401, 2.3287,
                       2.4906, 2.6321, 2.7573, 2.8694, 2.9709, 3),
                n  = 5:15) # Hosking and Wallis (1997, table 3.1)

  # lmrmeans is Hosking and Wallis (1997, eq. 3.1)
  lmrs <- t(matrix(c(t2, t3, t4), ncol=3))
  lmrmeans <- sapply(1:3, function(i) { mean(lmrs[i, ])        })
       tmp <- sapply(1:3, function(i) { lmrs[i,] - lmrmeans[i] })

  # Build sum of squares and cross products
  # Hosking and Wallis (1997, eq. 3.2)
  SSMAT <- matrix(nrow=3, ncol=3)
  for(j in 1:3) {
    for(k in j:3) {
      SSMAT[j,k] <- 0
      for(i in 1:nsite) {
        SSMAT[j,k] <- SSMAT[j,k] + tmp[i,j] * tmp[i,k]
      }
    }
  }
  SSMAT <- chol2inv(chol(SSMAT)) # Matrix inversion

  # Hosking and Wallis (1997, eq. 3.3)
  D <- rep(0, nsite)
  for(i in 1:nsite) {
    for(j in 1:3) {
      for(k in 1:3) {
        D[i] <- D[i] + tmp[i,j] * tmp[i,k] * SSMAT[j,k]
      }
    }
  }
  D <- nsite*D/3

  Dc <- Dcrit
  if(is.null(Dc)) {
    Dc <- ifelse(nsite >= 15, 3, tabDc$Dc[tabDc$n == nsite])
  }
  isD <- rep(FALSE, nsite)
  isD[D > Dc] <- TRUE

  Z <- qf(alpha1/nsite, 3, nsite - 4, lower.tail=FALSE)
  D1.by.Fdist <- (nsite - 1)*Z/(nsite - 4 + 3*Z)

  Z <- qf(alpha2/nsite, 3, nsite - 4, lower.tail=FALSE)
  D2.by.Fdist <- (nsite - 1)*Z/(nsite - 4 + 3*Z)

  star <- rep("-", nsite)
  star[D > D1.by.Fdist] <- "*"
  star[D > D2.by.Fdist] <- "**"

  z <- data.frame(site=site,
                  t2=round(t2, digits=lmrdigits),
                  t3=round(t3, digits=lmrdigits),
                  t4=round(t4, digits=lmrdigits),
                  Dmax=round((nsite - 1)/3, digits=digits),
                  Dalpha1=round(D1.by.Fdist, digits=digits),
                  Dalpha2=round(D2.by.Fdist, digits=digits),
                  Dcrit=Dc, D=round(D, digits=digits), isD=isD,
                  signif=star)
  if(sort) { ix <- order(z$D, decreasing=TRUE); z <- z[ix,] }

  return(z)
}

