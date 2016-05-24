rotcumvar <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  plus <- 1:n/(n-1) - cumsum(x^2)/sum(x^2)
  minus <- cumsum(x^2)/sum(x^2) - 0:(n-1)/(n-1)
  pmax(abs(plus), abs(minus))
}

testing.hov <- function(x, wf, J, min.coef=128, debug=FALSE) {
  n <- length(x)
  change.points <- NULL

  x.dwt <- dwt(x, wf, J)
  x.dwt.bw <- brick.wall(x.dwt, wf, method="dwt") 
  x.modwt <- modwt(x, wf, J)
  x.modwt.bw <- brick.wall(x.modwt, wf)

  for(j in 1:J) {
    cat("##### Level ", j, " #####", fill=TRUE)
    Nj <- n/2^j
    dwt.list <- list(dwt = (x.dwt.bw[[j]])[!is.na(x.dwt.bw[[j]])],
                     left = min((1:Nj)[!is.na(x.dwt.bw[[j]])]) + 1,
                     right = sum(!is.na(x.dwt.bw[[j]])))
    modwt.list <- list(modwt = (x.modwt.bw[[j]])[!is.na(x.modwt.bw[[j]])],
                       left = min((1:n)[!is.na(x.modwt.bw[[j]])]) + 1,
                       right = sum(!is.na(x.modwt.bw[[j]])))
    if(debug) cat("Starting recursion; using", dwt.list$left,
                  "to", dwt.list$right - 1, "...  ")
    change.points <-
      rbind(change.points,
            mult.loc(dwt.list, modwt.list, wf, j, min.coef, debug))
  }
  dimnames(change.points) <-
    list(NULL, c("level", "crit.value", "loc.dwt", "loc.modwt"))
  return(change.points)
}

mult.loc <- function(dwt.list, modwt.list, wf, level, min.coef, debug)
{
  Nj <- length(dwt.list$dwt)
  N <- length(modwt.list$modwt)
  crit <- 1.358
  change.points <- NULL
  
  if(Nj > min.coef) {
    ## test statistic using the DWT
    P <- cumsum(dwt.list$dwt^2) / sum(dwt.list$dwt^2)
    test.stat <- pmax((1:Nj) / (Nj-1) - P, P - (1:Nj - 1) / (Nj-1))
    loc.dwt <- (1:Nj)[max(test.stat) == test.stat]
    test.stat <- max(test.stat)

    ## location using the MODWT
    P <- cumsum(modwt.list$modwt^2) / sum(modwt.list$modwt^2)
    loc.stat <- pmax((1:N) / (N-1) - P, P - (1:N - 1) / (N-1))
    loc.modwt <- (1:N)[max(loc.stat) == loc.stat]

    if(test.stat > sqrt(2) * crit / sqrt(Nj)) {
      if(debug) cat("Accepted!", fill=TRUE)
      ## Left
      if(debug) cat("Going left; using", dwt.list$left,
                    "to", loc.dwt + dwt.list$left - 1, "...  ")
      temp.dwt.list <- list(dwt = dwt.list$dwt[1:(loc.dwt-1)],
                            left = dwt.list$left,
                            right = loc.dwt + dwt.list$left - 1)
      temp.modwt.list <- list(modwt = modwt.list$modwt[1:(loc.modwt-1)],
                              left = modwt.list$left,
                              right = loc.modwt + modwt.list$left - 1)
      change.points <-
        rbind(c(level, test.stat, loc.dwt + dwt.list$left,
                loc.modwt + modwt.list$left),
              Recall(temp.dwt.list, temp.modwt.list, wf, level, min.coef, debug))
      ## Right
      if(debug) cat("Going right; using", loc.dwt + dwt.list$left + 1,
                    "to", dwt.list$right, "...  ")
      temp.dwt.list <- list(dwt = dwt.list$dwt[(loc.dwt+1):Nj],
                            left = loc.dwt + dwt.list$left + 1,
                            right = dwt.list$right)
      temp.modwt.list <- list(modwt = modwt.list$modwt[(loc.modwt+1):N],
                              left = loc.modwt + modwt.list$left + 1,
                              right = modwt.list$right)
      change.points <-
        rbind(change.points,
              Recall(temp.dwt.list, temp.modwt.list, wf, level, min.coef, debug))
    }
    else
      if(debug) cat("Rejected!", fill=TRUE)
  }
  else
    if(debug) cat("Sample size does not exceed ", min.coef, "!",
                  sep="", fill=TRUE)

  return(change.points)
}

