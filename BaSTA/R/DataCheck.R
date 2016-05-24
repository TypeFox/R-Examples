DataCheck <-
function(object, studyStart, studyEnd, autofix = rep(0, 7), silent=TRUE) {

  if (length(autofix) < 7) {
    stop("Autofix specification should be a numerical vector of length 7.", 
          call. = FALSE)
  }
  Ti               <- studyStart
  Tf               <- studyEnd
  st               <- Ti:Tf
  nt               <- length(st)
  idnames          <- object[, 1]
  n                <- nrow(object)
  bd               <- as.matrix(object[, 2:3])
  Y                <- as.matrix(object[,1:nt+3]); colnames(Y) = st
  Tm               <- matrix(st, n, nt, byrow=TRUE)
	
  if(ncol(object) > nt + 3){
    Z              <- as.matrix(object[, (nt + 4):ncol(object)])
  } else {
    Z              <- matrix(1, n, 1)
  }

  
  # 1. Death before observations start
  type1            <- as.vector(which(bd[, 2] < Ti & bd[, 2] != 0))
  if (length(type1) != 0) {
    cat("The following rows deaths occur before observations start:\n")
    print(type1)
        
    # Actions - remove those rows from bd, Y and Z
    if (autofix[1] == 1) {
      bd           <- bd[-type1, ]
      Y            <- Y[-type1, ]
      idnames      <- idnames[-type1]
      Z            <- Z[-type1, ]
      Tm           <- Tm[-type1, ]
      n            <- nrow(Y)
      cat("These records have been removed from the Dataframe\n")
    }
  }
    
  # 2. No birth/death AND no obervations
  type2            <- as.vector(which(rowSums(bd) + rowSums(Y) == 0))
    if (length(type2) != 0) {
      cat("The following rows have no object (unknown birth, unknown death, and no observations):\n")
      print(type2)
        
      #Actions - remove those rows from bd, Y and Z
      if (autofix[2] == 1) {
      bd           <- bd[-type2, ]
      Y            <- Y[-type2, ]
      idnames      <- idnames[-type2]
      Z            <- Z[-type2, ]
      Tm           <- Tm[-type2, ]
      n            <- nrow(Y)
      cat("These records have been removed from the Dataframe\n")
    }
  }
    
  # 3. Birth after death 
  bd2              <- bd
  bd2              <- cbind(bd2, 1:n)

  type3 = as.vector(which(bd[, 1] > bd[, 2] & bd[, 1] != 0 & bd[, 2] != 0))
  if (length(type3) != 0) {
    cat("The following rows have birth dates that are later than their death dates:\n")
    print(type3)
        
    # Actions - remove the death, birth, both records?
    if (autofix[3] == 1) {
      bd[type3,2] = 0; cat("The death records have been replaced with 0.\n\n")
    } else if (autofix[3] == 2) {
    	bd[type3,1] = 0; cat("The birth records have been replaced with 0\n")
    } else if (autofix[3] == 3) {
    	bd[type3,1:2] = 0; cat("The birth and death records have been replaced with 0\n")
    }
  }
    
  # 4. Observations after death
  # Calculate first and last time observed: 
  st               <- Ti:Tf
  ytemp            <- t(t(Y) * st)
  lastObs          <- c(apply(ytemp, 1, max))
  tempDeath        <- bd[,2]
  tempDeath[tempDeath==0] <- Inf
  type4            <- as.vector(which(lastObs > tempDeath & tempDeath >= Ti))
  rm(tempDeath)
    
  if (length(type4) != 0) {
    cat("The following rows have observations that occur after the year of death:\n")
    print(type4)
        
    # Actions - remove spurious post-death observations
    if (autofix[4] == 1) {
      Ymd          <- ((Tm - bd[, 2]) * Y)[type4, ]
      Ymd[Ymd>0]   <- 0
      Ymd[Ymd<0]   <- 1
      Y[type4,]    <- Ymd
      cat("Observations that post-date year of death have been removed.\n\n")
    }
  }

  # 5. Observations before birth
  ytemp[ytemp == 0]<- Inf
  firstObs         <- c(apply(ytemp, 1, min))
  type5            <- as.vector(which(firstObs < bd[, 1]))
    
  if (length(type5) != 0) {
    cat("The following rows have observations that occur before the year of birth:\n")
    print(type5)
        
    # Actions - remove spurious pre-birth observations
    if (autofix[5] == 1) {
      Ymd          <- ((Tm - bd[, 1]) * Y)[type5, ]
      Ymd[Ymd>0]   <- 1
      Ymd[Ymd<0]   <- 0
      Y[type5,]    <- Ymd
      cat("Observations that pre-date year of birth have been removed.\n\n")
    }
  }
    
  # 6. Year of birth should be a zero in recapture matrix Y
  idb              <- which(bd[, 1] > 0 & bd[, 1] >= Ti & bd[, 1] <= Tf)
  bcol             <- bd[idb, 1] - Ti
  bpos             <- bcol * n + idb
  type6            <- as.vector(idb[which(Y[bpos] == 1)])

  if (length(type6) != 0) {
    cat("The following rows have a one in the recapture matrix in the birth year:\n")
    print(type6)
        
    # Actions - put a zero.
    if (autofix[6] == 1) Y[bpos] = 0
  }
    
  # 7. Year of death should be a zero in recapture matrix Y
#  idd              <- which(bd[, 2] > 0 & bd[, 2] >= Ti)
#  dcol             <- bd[idd, 2] - Ti
#  dpos             <- dcol * n + idd
#  type7            <- as.vector(idd[which(Y[dpos] == 1)])
#  if (length(type7) != 0) {
#    cat("The following rows have a one in the recapture matrix in the death year:\n")
#    print(type7)
        
    # Actions - put a zero.
#    if (autofix[7] == 1) Y[dpos] = 0
#  }   
  type7 <- c()
	n               <- nrow(Y)   

  # All OK:
  if (length(c(type1, type2, type3, type4, type5, type6, type7)) == 0) {
    cat("No problems were detected with the data.\n\n")
  }
#--------1---------2
  ok               <- length(c(type1, type2, type3, type4, 
                             type5, type6, type7)) == 0
    
    
  if (!silent) {
    cat(paste("*DataSummary*\n- Number of individuals         =",  
        format(n, big.mark = ',', width = 7), "\n"))
    cat(paste("- Number with known birth year  =", 
        format(sum(bd[, 1] != 0), big.mark = ',', width = 7), "\n"))
    cat(paste("- Number with known death year  =", 
        format(sum(bd[, 2] != 0), 
        big.mark = ',', width = 7), "\n"))
    cat(paste("- Number with known birth\n",
        "AND death years                =", 
        format(sum(bd[, 2] != 0 & bd[, 1] != 0), big.mark = ",", 
        width = 7), "\n\n"))
    cat(paste("- Total number of detections\n",  
        "in recapture matrix            =", 
        format(sum(Y), big.mark = ",", width = 7), "\n\n"))
    cat(paste("- Earliest detection time       =", 
        format(min(ytemp), width = 7), "\n"))
    cat(paste("- Latest detection time         =", 
        format(max(ytemp[ytemp != Inf]), width = 7), "\n"))
    if(length(which(bd[,1] > 0)) > 0){
      cat(paste("- Earliest recorded birth year  =", 
          format(min(bd[bd[, 1] > 0, 1]), width = 7), "\n"))
      cat(paste("- Latest recorded birth year    =", 
          format(max(bd[, 1]), width = 7), "\n"))
    }
    if(length(which(bd[, 2] > 0)) > 0){
      cat(paste("- Earliest recorded death year  =", 
          format(min(bd[bd[, 2] > 0, 2]), width = 7), "\n"))
		cat(paste("- Latest recorded death year    =", 
		    format(max(bd[, 2]),width = 7), "\n"))
    }
  }

  if (ncol(object)>nt+3) {
    newData        <- data.frame(idnames, bd, Y, Z)
  } else {
    newData        <- data.frame(idnames, bd, Y)
  } 
  return(list(ok       = ok, 
              newData  = newData, 
              type1    = type1, 
              type2    = type2, 
              type3    = type3, 
              type4    = type4, 
              type5    = type5, 
              type6    = type6, 
              type7    = type7))
}

