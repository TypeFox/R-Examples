# compare blocks qualitatively according to the criterion of Elhaik et al. without eps
# true positive if both borders match
# false positive if either detected border was not true
# false negative if either true border was not detected
# sensitivity rate = true pos. / ( true pos. + false neg. ) = true pos. / true number
# precision rate = true pos. / ( true pos. + false pos. ) = true pos. / detected number
# borders are considered correctly identified if they differ by less than dist = 5e3


# also compare blocks quantitatively
# false negative sensitive localization error: for each true block's mid-point find into which estimated block it falls, and sum distances of the respective borders
# false positive sensitive localization error: the other way round
# total localisation error: sum
compareBlocks <- function(truth, estimate, dist = 5e3) {
  if(is.data.frame(truth)) truth <- list(truth)
  if(is.data.frame(estimate)) estimate <- list(estimate)
  stopifnot(length(truth) == length(estimate))
  
  ret <- data.frame(true.num = rep(NA, length(truth)), est.num = rep(NA, length(truth)), true.pos = rep(NA, length(truth)), false.pos = rep(NA, length(truth)), false.neg = rep(NA, length(truth)), sens.rate = rep(NA, length(truth)), prec.rate = rep(NA, length(truth)), fpsle = rep(NA, length(truth)), fnsle = rep(NA, length(truth)), total.le = rep(NA, length(truth)))
  
  for(i in 1:length(truth)) {
    # qualitatively
    ret$true.num[i] <- nrow(truth[[i]])
    ret$est.num[i] <- nrow(estimate[[i]])
    ret$true.pos[i] <- sum(sapply(1:nrow(estimate[[i]]), function(j) any(abs(estimate[[i]]$rightEnd[j] - truth[[i]]$rightEnd) <= dist & abs(estimate[[i]]$leftEnd[j] - truth[[i]]$leftEnd) <= dist)))
    ret$false.pos[i] <- ret$est.num[i] - ret$true.pos[i]
    ret$prec.rate[i] <- ret$true.pos[i] / ret$est.num[i]
    ret$false.neg[i] <- ret$true.num[i] - ret$true.pos[i]
    ret$sens.rate[i] <- ret$true.pos[i] / ret$true.num[i]
    
    # quantitatively
    truth[[i]]$mid <- ( truth[[i]]$leftEnd + truth[[i]]$rightEnd ) / 2
    estimate[[i]]$mid <- ( estimate[[i]]$leftEnd + estimate[[i]]$rightEnd ) / 2
    truth[[i]]$match <- sapply(truth[[i]]$mid, function(m) min(which(estimate[[i]]$rightEnd >= m)))
    estimate[[i]]$match <- sapply(estimate[[i]]$mid, function(m) min(which(truth[[i]]$rightEnd >= m)))
    ret$fnsle[i] <- sum(abs(truth[[i]]$leftEnd - estimate[[i]]$leftEnd[truth[[i]]$match]) + abs(truth[[i]]$rightEnd - estimate[[i]]$rightEnd[truth[[i]]$match])) / 2
    ret$fpsle[i] <- sum(abs(estimate[[i]]$leftEnd - truth[[i]]$leftEnd[estimate[[i]]$match]) + abs(estimate[[i]]$rightEnd - truth[[i]]$rightEnd[estimate[[i]]$match])) / 2
    ret$total.le[i] <- ret$fnsle[i] + ret$fpsle[i]
  }
  ret
}
