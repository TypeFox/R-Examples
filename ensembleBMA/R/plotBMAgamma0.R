plotBMAgamma0 <-
function (WEIGHTS, MEAN, VAR, PROB0, obs = NULL, exchangeable = NULL, power = 1)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
  k <- length(WEIGHTS)

  W <- WEIGHTS

  Q <- matrix( NA, 2, k + 1)


  n <- if (diff(range(VAR)) < 1.e-4) 100 else 1000
  
  bot <- 1/n
  top <- (n-1)/n

  Q[1,k+1] <- quantBMAgamma0( bot,  W, MEAN, VAR, PROB0)
  Q[2,k+1] <- quantBMAgamma0( top,  W, MEAN, VAR, PROB0)

  for (j in 1:k){
     W[] <- 0
     W[j] <- 1
     Q[1,j] <- quantBMAgamma0( bot, W, MEAN, VAR, PROB0)
     Q[2,j] <- quantBMAgamma0( top, W, MEAN, VAR, PROB0)
  }

  r <- range(Q)


  lo <- quantBMAgamma0( .1,  WEIGHTS, MEAN, VAR, PROB0)
  med <- quantBMAgamma0( .5,  WEIGHTS, MEAN, VAR, PROB0)
  up <- quantBMAgamma0( .9,  WEIGHTS, MEAN, VAR, PROB0)

  n <- 1000

  if (!is.null(obs) && !is.na(obs)) r <- range(c(r,obs))

  if (is.null(exchangeable)) exchangeable <- 1:k

  tex <- table(exchangeable)
  lex <- length(tex)

  FORC <- matrix( NA, n, lex + 1)
  x <- seq(from = r[1], to = r[2], length = n)

  RATE <- MEAN/VAR
  for (l in 1:lex) {
     j <- which(exchangeable == l)[1]
     FORC[,l] <- dgamma( x, shape = RATE[j]*MEAN[j], rate = RATE[j])
     FORC[,l] <- tex[l]*(WEIGHTS[j]*(1-PROB0[j])*FORC[,j])
  }

  for (i in 1:n) FORC[i,lex+1] <- sum(FORC[i,1:lex])

  FORC[is.infinite(FORC)] <- 0

  const <- sum(WEIGHTS*(1-PROB0))/max(FORC[,lex+1])
  
  FORC <- const*FORC
  
#  matplot(FORC,PROB)

ylim <- range(c(FORC,sum(WEIGHTS*PROB0)))
xlim <- range(c(0,x))

#if (!is.null(obs) && (obs <= lo || obs >= up)) return(invisible())

  xlab <- "Precipitation"
  if (power != 1) {
      if (power == 1/3) {
            xlab <- "Cube Root of Precipitation"
          }
        else if (power == 1/2) {
              xlab <- "Square Root of Precipitation"
            }
        else if (power == 1/4) {
              xlab <- "Fourth Root of Precipitation"
            }
        else {
              xlab <- paste( xlab, " to the ", round(power,3), " power")
            }
    }
  
plot( c(0,x), c(0,FORC[,lex+1]), type = "l", col = "black", ylim = ylim,
      xlab = xlab, 
ylab = "Prob No Precip and Scaled PDF for Precip", lwd = 3)

segments( 0, 0, 0, sum(WEIGHTS*PROB0), col = "black", lwd = 3)

abline( v = lo, col = "black", lty = 2)
abline( v = med, col = "black")
abline( v = up, col = "black", lty = 2)

if (!is.null(obs) && !is.na(obs)) abline( v = obs, col = "orange", lwd = 3)

colors <- rainbow(lex)
for (l in 1:lex) {
  lines( x, FORC[,l], col = colors[l], lty = 1)
}

lines( x, FORC[,lex+1], col = "black", lwd = 3)

invisible(FORC)
}

