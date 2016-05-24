diagnose <- function(object, data, id.vars, suspect.var,
                     suspect.range = NULL){ 

  out <- list()

  if(is.matrix(object) || is.data.frame(object)){
    object <- list(as.data.frame(object))
  }

  if(is.matrix(data)){
    data <- as.data.frame(data)
  }

  if(!is.character(id.vars)){
    stop("id.vars is not a character string, or vector of character
strings.  See documentation and respecify id.vars.")
  }

  gp.names <- array(NA)

  data.diag <- data[, (names(data) %in% c(id.vars, suspect.var))]  

  lev <- length(id.vars)

####  for(i in 1:length(object)){  
  for(i in 1:length(object$assg)){  
    ## assignment object
    ####assg.gp <- object[[i]]
    assg.gp <- as.matrix(object$assg[[i]])
    ## data from group
    data.gp <- data.diag[data.diag[,id.vars[lev]] %in%
                         as.matrix(assg.gp),]

    if(nrow(data.gp) != 0){
      tmp <- data.gp[, !(names(data.gp) %in% id.vars)]
      d.mat <- expand.grid(tmp, tmp)

      diffs <- abs(d.mat[,1]-d.mat[,2])
      suspect.vec <- (suspect.range[1] <= diffs) & (diffs <= suspect.range[2])

      storage <- as.data.frame(matrix(NA, 1, 2*lev+1))
      ct <- 0

      for(j in which(suspect.vec)){
        ct <- ct+1
        u1 <- j%%nrow(data.gp)
        u1[u1==0] <- nrow(data.gp)
        u2 <- ceiling(j/nrow(data.gp))
        storage[ct,] <- cbind(data.gp[u1,id.vars],
                              data.gp[u2,id.vars], diffs[j])
      }

      tmp.na <- array(NA)

      for(k in 1:nrow(storage)){
        col1 <- ceiling(which(object$assg[[i]] ==
                              storage[k,1])/nrow(object$assg[[i]]))
        col2 <- ceiling(which(object$assg[[i]] ==
                              storage[k,(lev+1)])/nrow(object$assg[[i]]))
        if(length(col1)==0 || length(col2)==0 || (col1==col2)){
          tmp.na <- append(tmp.na,k)
        }
      }

      if(sum(is.na(tmp.na))!=length(tmp.na)){
        tmp.na <- tmp.na[2:length(tmp.na)]
        storage <- storage[-tmp.na,]
      }

      if(sum(is.na(storage))==prod(dim(storage))){
        storage <- "No units with different treatments are suspect."
      }     

      tmp.dup <- NULL
      if(!is.null(nrow(storage))){
        ## cut duplicates
        tmp.dup <- rep(TRUE, nrow(storage))
        for(ll in 1:(nrow(storage)-1)){
          for(mm in (ll+1):nrow(storage)){
            if(sum(c(storage[ll,1], storage[ll, lev+1]) %in%
                   storage[mm,1:(ncol(storage)-1)]) == 2){
              tmp.dup[mm] <- FALSE
            }
          }
        }
      }

      if(length(tmp.dup)>0){
        storage <- storage[tmp.dup,]
      }

      if(sum(is.na(storage))==prod(dim(storage))){
        storage <- "No units with different treatments are suspect."
      }
    }else{
      storage <- "No units with different treatments are suspect."
    }

    if(!is.null(nrow(storage))){
      ## renumber rows
      rownames(storage) <- 1:nrow(storage)
      ## order columns by value
      tmp1 <- storage[,1:lev]
      storage[,1:lev] <- storage[,(lev+1):(2*lev)]
      storage[,(lev+1):(2*lev)] <- tmp1
      ## sort rows, ascending by distance
      o <- order(storage[,ncol(storage)])
      storage <- storage[o,]
      ## name rows
      rownames(storage) <- 1:nrow(storage)
      ## name columns
      names(storage)[ncol(storage)] <- "Difference"
      reps <- floor(ncol(storage)/2)
      names(storage)[1:(ncol(storage)-1)] <- rep(paste("Unit ",
                                                 1:2, sep=""),
                                                 each=reps)
    }

    gp.names[i] <- names(object$assg)[i]
    out[[i]] <- storage
  }

  names(out) <- gp.names
  output <- list(diagnose = out)
  output$call <- match.call()
  output$suspect.var <- suspect.var
  output$suspect.range <- suspect.range
  class(output) <- "diagnose"
  return(output)
}
