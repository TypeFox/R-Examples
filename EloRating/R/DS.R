# library(EloRating)
# data(adv)
# el <- elo.seq(adv$winner, adv$loser, adv$Date)
# interactionmatrix <- m <- creatematrix(el)

DS <- function(interactionmatrix, prop=c("Pij", "Dij")){
  if(length(intersect(rownames(interactionmatrix), colnames(interactionmatrix))) < ncol(interactionmatrix)) stop("not a square matrix")
  if(length(intersect(rownames(interactionmatrix), colnames(interactionmatrix))) < nrow(interactionmatrix)) stop("not a square matrix")
  
  # create an index for the matrix cells with unobserved dyads and the matrix diagonal
  summatrix <- interactionmatrix + t(interactionmatrix); diag(summatrix) <- 0
  summatrix <- replace(summatrix, summatrix==0, NA); l1 <- which(is.na(summatrix), arr.ind=TRUE)
  
  # create matrix with Dij OR Pij Indices
  if(prop[1]=="Dij") propmatrix <- (interactionmatrix+0.5) / (t(interactionmatrix)+interactionmatrix+1)
  if(prop[1]=="Pij") propmatrix <- interactionmatrix / (t(interactionmatrix)+interactionmatrix)
  
  # replace Dij/Pij-values for the diagonal and unobserved dyads with zero (by definition)
  propmatrix <- replace(propmatrix, l1, 0)
  
  # calculate w, w2, l, l2, and then DS and normDS
  w  <- rowSums(propmatrix)
  w2 <- propmatrix %*% w
  l  <- rowSums(t(propmatrix))
  l2 <- t(propmatrix) %*% l
  DS <- w + w2 - l - l2
  normDS <- ((DS+((length(DS)) * (length(DS)-1))/2)) / length(DS)
  
  res <- data.frame(ID=rownames(interactionmatrix), DS=DS, normDS=normDS)
  res <- res[order(res$DS, decreasing=T), ]; rownames(res) <- NULL
  
  return(res)
}

#DS(m, "Pij")
