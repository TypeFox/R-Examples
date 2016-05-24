# Implicit functions =============
# Function list:
  # PathBuild
  # SumNotNA
  # MouthSource
  # RelPos
  # DigitWeight
  # RelPosMatrix
  # RowCal
  # TitlePaste


# Path building -----------------------------------------------------------

PathBuild <- function(river, parent, OBN){
  p <- cbind(river, parent)
  for (i in 1:OBN) {
    p <- cbind(p,p[,2][match(p[,1+i],p[,1])])
    if ((sum(!is.na(p[,i+2])) == 0)) break
  }
  p <- p[,1:(ncol(p)-1)]
  p
}


# Calculate the sum without NA -------------------

SumNotNA <- function(x){
  sum(x,na.rm=TRUE)
}

# Calculation the location of river Mouth and river Source ----------------------------

MouthSource <- function(path, riverlayout, OBN){
  p <- path
  for (i in 1:OBN){
    p[p == riverlayout$river[i]] <- riverlayout$distance[i]
  }
  class(p) <- "numeric"
  m <- apply(p,1,SumNotNA)
  s <- m + riverlayout$length
  mouthsource <- cbind(rmouth = m, rsource = s)
  row.names(mouthsource) <- riverlayout$river
  mouthsource
}

# Relative position matrix -------------

RelPos <- function(path, riverlayout, OBN, DIGITMAX){
  p <- path
  for (i in 1:OBN){
    p[p == riverlayout$river[i]] <- riverlayout$position[i]
    p[i,1:sum(!is.na(p[i,]))] <- rev(p[i,][!is.na(p[i,])])
  }
  row.names(p) <- riverlayout$River
  colnames(p) <- paste("digit",c(0:DIGITMAX),sep="")
  p
}

# Digit weight-----------------

DigitWeight <- function(DIGITMAX){
  digitweight <- matrix(NA,DIGITMAX,1)
  for (i in 1:DIGITMAX){
    digitweight[i,1] <- 10^(9-i) # Actually 4 is a much better choice than 10.
  }
  digitweight
}

# Relative position numeric matrix ---------------------  

RelPosMatrix <- function(pos, DIGITMAX){
  
  if (DIGITMAX>1){
    p <- pos[,2:(DIGITMAX+1)]
  }else{
    p <- as.matrix(pos[,2:(DIGITMAX+1)])
  }
  #  for (i in 2:(DIGITMAX-1)){
  p[,1][is.na(p[,1])] <- 0
  p[,1][p[,1] == "L"] <- -1
  p[,1][p[,1] == "R"] <- 1
  if (DIGITMAX>1){
    p[,2:DIGITMAX][p[,2:DIGITMAX] == "L"] <- 1
    p[,2:DIGITMAX][is.na(p[,2:DIGITMAX])] <- 2
    p[,2:DIGITMAX][p[,2:DIGITMAX] == "R"] <- 3
  }
  #  }
  class(p) <- "numeric"
  p
}

# Sort, and optimise row number ----------------

RowCal <- function(posmatrix, digitweight, riverlayout, path, OBN, DIGITMAX){
  if (DIGITMAX > 2){
    s <- posmatrix[,2:DIGITMAX] %*% digitweight * posmatrix[,1]
  }else if(DIGITMAX > 1){
    s <- data.frame(s = posmatrix[,2] * digitweight[2,] * posmatrix[,1])
  }else{
    s <- posmatrix[,1] %*% digitweight
  }
  
  colnames(s) <- "s"
  OBNLEFT <- length(s[s<0])
  OBNRIGHT <- length(s[s>0])
  
  if (OBNLEFT*OBNRIGHT != 0){
    row <- c(c(1:length(s[s<0])),0,c(-1:-length(s[s>0])))
  }else if(OBNLEFT == 0 & OBNRIGHT != 0){
    row <- c(0,c(-1:-length(s[s>0])))
  }else if(OBNRIGHT == 0 & OBNLEFT != 0){
    row <- c(c(1:length(s[s<0])),0)  
  }else if(OBNLEFT == 0 & OBNRIGHT == 0){
    row <- 0
  }
  
  k <- data.frame(MouthSource(path, riverlayout, OBN),s)
  k <- data.frame(k[order(k$s,k$rmouth),],row)
  
  # Optimise row numbers
  
  # Optimise left side
  if (OBNLEFT > 1){
    for (i in 2:OBNLEFT){
      j = i
      while (all(k$rmouth[i] > k$rsource[k$row == (j-1)]) | 
               all(k$rsource[i] < k$rmouth[k$row == (j-1)])){
        k$row[i] <- j - 1
        j = j - 1
      }
    }
    
    for (i in 2:OBNLEFT){
      for (j in 1:k$row[i]){
        if (all(k$rmouth[i] > k$rsource[k$row == j]) | 
              all(k$rsource[i] < k$rmouth[k$row == j])){
          k$row[i] <- j
          break
        }
      }
    }
    
    
    
  }
  
  
  
  # Optimise right side
  if (OBNRIGHT > 1){
    for (i in (OBNLEFT+3):OBN){
      j = i
      while (all(k$rmouth[i] > k$rsource[k$row == (OBNLEFT+2-j)]) | 
               all(k$rsource[i] < k$rmouth[k$row == (OBNLEFT+2-j)])){
        k$row[i] <- OBNLEFT+2-j
        j = j - 1
      }
    } 
    
    for (i in (OBNLEFT+3):OBN){
      for (j in -1:k$row[i]){
        if (all(k$rmouth[i] > k$rsource[k$row == j]) | 
              all(k$rsource[i] < k$rmouth[k$row == j])){
          k$row[i] <- j
          break
        }
      }
      
    } 
  }
  
  #  row <- matrix(k$row, dimnames = list(rownames(k),"Row"))
  row <- data.frame(river = rownames(k), row = -k$row)
  
  row
}

# TitlePaste =========================

TitlePaste <- function(x){
  
  if(length(x)==1){
    x[[1]]
  }else{
    y <- NA
    
    for(i in 1:length(x)){
      y[i] <- paste(x[[i]], collapse = "/")
      
    }
    
    y
  }
  
}

