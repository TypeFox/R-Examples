`drawnames` <- function(nrow.pt, ncol.pt, include.race.contrib.turnout=FALSE)
{
  dim.SI.0 <- nrow.pt*(ncol.pt-1)
  rownames.def.pt <- c("b", "w", "h", "a", "o")
  colnames.def.pt <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")
  rownames.use.pt <- rownames.def.pt[1:nrow.pt]
  colnames.use.pt <- rep("N", ncol.pt)
  colnames.use.pt[1:(ncol.pt-1)] <- colnames.def.pt[1:(ncol.pt-1)]
  vec1 <- rep(rownames.use.pt, 1, each = (ncol.pt-1))
  vec2 <- rep(colnames.use.pt[1:(ncol.pt-1)], nrow.pt)
  vec.comb <- rep("Y", dim.SI.0)
  for(ii in 1:dim.SI.0){
    vec.comb[ii] <- paste(vec1[ii], vec2[ii], sep = "")
  }
  pos.names <- 1
  length.names.final.temp <- 2*dim.SI.0 + (dim.SI.0+1)*dim.SI.0/2 + nrow.pt*ncol.pt + nrow.pt
  if(include.race.contrib.turnout){
    length.names.final.temp <- length.names.final.temp + nrow.pt
  }
  names.final.temp <- rep("Z", length.names.final.temp)
  for(ii in 1:dim.SI.0){
    names.final.temp[pos.names] <- paste("mu", vec.comb[ii], sep = "")
    pos.names <- pos.names + 1
  }
  for(ii in 1:dim.SI.0){
    names.final.temp[pos.names] <- paste("sg", vec.comb[ii], sep = "")
    pos.names <- pos.names + 1
  }
  for(ii in 1:(dim.SI.0-1)){
    for(jj in (ii+1):dim.SI.0){
      names.final.temp[pos.names] <- paste("rho", vec.comb[ii], vec.comb[jj], sep = "")
      pos.names <- pos.names + 1
    }
  }
  vec3 <- rep(rownames.use.pt, 1, each = ncol.pt)
  vec4 <- rep(colnames.use.pt, nrow.pt)
  for(ii in 1:(nrow.pt*ncol.pt)){
    names.final.temp[pos.names] <- paste("NN", vec3[ii], vec4[ii], sep = "")
    pos.names <- pos.names + 1
  }
  for(ii in 1:dim.SI.0){
    names.final.temp[pos.names] <- paste("Per", vec.comb[ii], sep = "")
    pos.names <- pos.names + 1
  }
  for(ii in 1:nrow.pt){
    names.final.temp[pos.names] <- paste(rownames.use.pt[ii], "TO", sep = "")
    pos.names <- pos.names + 1
  }
  if(include.race.contrib.turnout){
    for(ii in 1:nrow.pt){
      names.final.temp[pos.names] <- paste("PerVtrs=", rownames.use.pt[ii], sep = "")
      pos.names <- pos.names + 1
    }
  }

  names.final.temp
}

