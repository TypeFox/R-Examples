plotBMAnormal <-
function (WEIGHTS, MEAN, SD, obs = NULL, exchangeable = NULL)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
  k <- length(WEIGHTS)

  Q <- matrix( NA, 2, k + 1)

  W <- WEIGHTS

  n <- 1000

  bot <- 1/n
  top <- (n-1)/n

  Q[1,k+1] <- quantBMAnormal( bot,  WEIGHTS, MEAN, SD)
  Q[2,k+1] <- quantBMAnormal( top,  WEIGHTS, MEAN, SD)

  for (j in 1:k){
     W[] <- 0
     W[j] <- 1
     Q[1,j] <- quantBMAnormal( bot, W, MEAN, SD)
     Q[2,j] <- quantBMAnormal( top, W, MEAN, SD)
  }

  r <- range(Q)

  if (!is.null(obs) && !is.na(obs)) r <- range( c(r, obs))

  if (is.null(exchangeable)) exchangeable <- 1:k

  tex <- table(exchangeable)
  lex <- length(tex)

  FORC <- matrix( NA, n, lex + 1)
  x <- seq(from = r[1], to = r[2], length = n)

  for (l in 1:lex) {
     j <- which(exchangeable == l)[1]
     FORC[,l] <- tex[l]*(WEIGHTS[j]*dnorm( x, MEAN[j], SD[j]))
  }

  for (i in 1:n) FORC[i,lex+1] <- sum(FORC[i,1:lex])

#  matplot(FORC,PROB)

med <- quantBMAnormal( .5,  WEIGHTS, MEAN, SD)
up <- quantBMAnormal( .9,  WEIGHTS, MEAN, SD)
lo <- quantBMAnormal( .1,  WEIGHTS, MEAN, SD)

ylim <- range(FORC)

#if (!is.null(obs) && (obs <= lo || obs >= up)) return(invisible())

plot( x, FORC[,lex+1], type = "l", col = "black", ylim = ylim,
      xlab = "Temperature", ylab = "Probability Density", lwd = 3)

abline( v = med, col = "black")
abline( v = lo, col = "black", lty = 2)
abline( v = up, col = "black", lty = 2)

if (!is.null(obs) && !is.na(obs)) abline( v = obs, col = "orange", lwd = 3)

colors <- rainbow(lex)
for (l in 1:lex) {
  lines( x, FORC[,l], col = colors[l], lty = 1)
}

lines( x, FORC[,lex+1], col = "black", lwd = 3)

#aux <- list(...)
#print(c(lat = aux$lat, lon = aux$lon))

invisible(FORC)
}

