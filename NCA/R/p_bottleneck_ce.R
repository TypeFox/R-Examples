p_bottleneck_fdh <-
function (loop.data, mpy, peers, cutoff, bottleneck.x) {
  bottleneck.x <- p_validate_bottleneck(bottleneck.x, "x")

  if (is.vector(peers) || length(peers) == 2) {
    # if peers is a vector there is only one peer
    if (cutoff == 0) {
      mpx <- matrix(NA, nrow=length(mpy), ncol=1)
    } else if (cutoff == 1) {
      mpx <- matrix(loop.data$x.low, nrow=length(mpy), ncol=1)
    } else if (cutoff == 2) {
      mpx <- matrix(peers[1, 1], nrow=length(mpy), ncol=1)
    }
  } else {
    x.peers <- peers[,1]
    y.peers <- peers[,2]
    mpx <- matrix(nrow=length(mpy), ncol=1)

    for (j in 1:length(mpy)) {
      # search the peer that is closest above to the desired outcome,
      # and select it corresponding x value of that peer
      index <- which(y.peers >= mpy[j,1])[1]
      if (is.na(index)) {
        mpx[j,1]  <- NA
      } else {
        mpx[j,1]  <- x.peers[index]
      }
    }

    if (cutoff == 0) {
      mpx [mpx <= loop.data$x.low]  <- NA
      mpx [mpx >  loop.data$x.high] <- NA
    } else if (cutoff == 1) {
      mpx [mpx <= loop.data$x.low]  <- loop.data$x.low
      mpx [mpx >  loop.data$x.high] <- loop.data$x.high
    }
  }

  if (p_bottleneck_id(bottleneck.x) == 3) {
    max.value <- loop.data$x.high
  } else {
    max.value <- 100
  }

  # Display Xs as percentage (either cutoff or 0-high) or percentile
  if (p_bottleneck_id(bottleneck.x) == 1) {
    mpx <- 100 * (mpx - loop.data$x.low) / (loop.data$x.high - loop.data$x.low)
  } else if (p_bottleneck_id(bottleneck.x) == 2) {
    mpx <- 100 * mpx / loop.data$x.high
  } else if (p_bottleneck_id(bottleneck.x) == 4) {
    percentile <- ecdf(sort(loop.data$x))
    mpx <- matrix(100 * percentile(mpx), ncol=1)
  }

  return (p_pretty_mpx(loop.data, mpx, max.value))
}

p_bottleneck_vrs <-
function (loop.data, mpy, peers, cutoff, bottleneck.x) {
  calculate_x <- function(y, p1, p2) {
    return(p1[1] + (y - p1[2]) * (p1[1] - p2[1]) / (p1[2] - p2[2]))
  }

  bottleneck.x <- p_validate_bottleneck(bottleneck.x, "x")

  if (is.vector(peers) || length(peers) == 2) {
    # if peers is a vector there is only one peer
    if (cutoff == 0) {
      mpx <- matrix(NA, nrow=length(mpy), ncol=1)
    } else if (cutoff == 1) {
      mpx <- matrix(loop.data$x.low, nrow=length(mpy), ncol=1)
    } else if (cutoff == 2) {
      mpx <- matrix(peers[1, 1], nrow=length(mpy), ncol=1)
    }
  } else {
    mpx <- matrix(nrow=length(mpy), ncol=1)
    for (j in 1:length(mpy)) {
      idx1 <- tail(which(peers[,2] <= mpy[j,1]), n=1)
      idx2 <- which(peers[,2] >= mpy[j,1])[1]

      if (length(idx1) == 0 || length(idx2) == 0) {
        mpx[j, 1] <- NA
      } else if (idx1 == idx2) {
        if (j == 1) {
          mpx[j, 1] <- NA
        } else {
          mpx[j, 1] <- peers[idx1, 1]
        }
      } else {
        p1 <- peers[idx1, ]
        p2 <- peers[idx2, ]
        mpx[j, 1] <- calculate_x(mpy[j,1], p1, p2)
      }
    }
  }

  if (p_bottleneck_id(bottleneck.x) == 3) {
    max.value <- loop.data$x.high
  } else {
    max.value <- 100
  }

  # Display Xs as percentage (either cutoff or 0-high) or percentile
  if (p_bottleneck_id(bottleneck.x) == 1) {
    mpx <- 100 * (mpx - loop.data$x.low) / (loop.data$x.high - loop.data$x.low)
  } else if (p_bottleneck_id(bottleneck.x) == 2) {
    mpx <- 100 * mpx / loop.data$x.high
  } else if (p_bottleneck_id(bottleneck.x) == 4) {
    percentile <- ecdf(sort(loop.data$x))
    mpx <- matrix(100 * percentile(mpx), ncol=1)
  }

  return (p_pretty_mpx(loop.data, mpx, max.value))
}
