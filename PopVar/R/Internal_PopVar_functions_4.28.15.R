## Group of functions that PopVar uses internally

par.position <- function(crossing.table, par.entries){ # Used when a crossing table is defined
  
  par.pos <- matrix(nrow = nrow(crossing.table), ncol = 2)
  crosses.possible <- matrix(nrow = nrow(crossing.table), ncol = 2)
  for(i in 1:nrow(crossing.table)){
    par1 <- as.character(crossing.table[i,1])
    par2 <- as.character(crossing.table[i,2])
    
    if(par1 %in% par.entries & par2 %in% par.entries){
      par.pos[i,1] <- which(par.entries == par1)
      par.pos[i,2] <- which(par.entries == par2)
      crosses.possible[i,1] <- par1
      crosses.possible[i,2] <- par2  
    }
  }
  
  par.pos <- par.pos[which(!is.na(par.pos[,1])), ]
  crosses.possible <- crosses.possible[which(!is.na(crosses.possible[,1])), ]
  
  for.dup <- paste(par.pos[,1], ".", par.pos[,2], sep=""); duplicated <- which(duplicated(for.dup))
  if(length(duplicated) > 0){
    par.pos <- par.pos[-duplicated, ]
    crosses.possible <- crosses.possible[-duplicated, ]
  }
  
  return(list(parent.positions=par.pos, crosses.possible=crosses.possible))
}

par.name <- function(crossing.mat, par.entries){ ## Used when all combinations of parents are crossed
  crosses.possible <- matrix(nrow = nrow(crossing.mat), ncol = 2)
  for(i in 1:nrow(crossing.mat)){
    crosses.possible[i,1] <- par.entries[crossing.mat[i,1]]
    crosses.possible[i,2] <- par.entries[crossing.mat[i,2]]
  }
  return(crosses.possible)
}

tails <- function(GEBVs, tail.p){ #Calculates means of tails; set tail.p to the proportion of the populaiton you want to take the mean of, default is 10%
  u.top <- mean(GEBVs[which(GEBVs >= quantile(GEBVs, 1-tail.p))], na.rm=T)
  u.bot <- mean(GEBVs[which(GEBVs <= quantile(GEBVs, tail.p))], na.rm=T)
  
  return(rbind(u.top, u.bot))
}

maf.filt <- function(G){
  G_noNA <- G[!is.na(G)]
  len_G <- length(G_noNA)
  maf <- (length(which(G_noNA == 1)) / len_G) + .5*(length(which(G_noNA == 0)) / len_G)
  
  return(maf)  
}