bounds <- function(formula, data, rows, column, excluded=NULL,
                   threshold = 0.9, total=NULL){ 

  "%wo%" <- function(x,y){x[!x %in% y]}
  
  D <- model.frame(formula, data = data)
  G <- D[[2]]
  T <- D[[1]]

  idx.r <- apply(G,1,sum)
  idx.c <- apply(T,1,sum)

  prop.rows <- list()

  countG <- countT <- propG <- propT <- FALSE


  if(!is.null(total)){
    if(!is.numeric(total)){
      if(is.character(total)){ total <- data[[total]]}else{
        total <- data[[deparse(substitute(total))]]}}}
  
  if(all(as.integer(T)==T) && all(T)>=0){    
    countT <- TRUE
  }
  else{propT <- TRUE}
  
  if(all(as.integer(G)==G) & all(G)>=0){
    countG <- TRUE
  }
  else{propG <- TRUE}

  if(propT && is.null(total)){
    stop("columns are proportions but no unit totals are provided -
please respecify data")
  }

  if(propG && is.null(total)){
    stop("rows are proportions but no unit totals are provided -
please respecify data")
  }

  if(propT){
    idx.pc <- apply(T,1,sum)
    if(!all(round(idx.pc, digits=3)==1)){
      stop("column marginals are proportions that do not sum to 1 - please
respecify data")
    }
  }

  if(propG){
    idx.pr <- apply(G,1,sum)
    if(!all(round(idx.pr, digits=3)==1)){
      stop("row marginals are proportions that do not sum to 1 - please
respecify data")
    }
  }
  
  if(countT & propG){
    if(!all(0 <= G && G <= 1)){
      stop("row proportions are not within [0,1] - please respecify
data")
    }
    else{
      G <- G*total
      propG <- FALSE
      countG <- TRUE
    }
  }

  if(propT & countG){
    if(!all(0 <= T && T <= 1)){
      stop("column proportions are not within [0,1] - please respecify
data")
    }
    else{
      T <- T*total
      propT <- FALSE
      countT <- TRUE
    }
  }

  if(propT & propG){
    G <- G*total
    T <- T*total
    propT <- propG <- FALSE
    countT <- countG <- TRUE
  }
  
  if(countT & countG){
    if(all(idx.r == idx.c)){
      for(i in rows){
        prop.rows[[i]] <- G[,i]/idx.r
      }
    }
    else{
      stop("row and column count totals unequal in some precincts -
please respecify data")
    }
  }

  idx <- list()

  for(i in rows){
    idx[[i]] <- which(prop.rows[[i]] >= threshold)
  }

  if(all(lapply(idx, length) == 0)){
    stop("no precincts satisfy homogeneity threshold - try lowering
threshold")
  }    

  bounds.out <- list()
  bound.names <- paste(rows, column, sep=".")
  intersection <- list()
  count <- 1
  
  for(i in rows){
    if(length(idx[[i]])==0){
      bounds.out[[bound.names[count]]] <- NA
    }
    else{
      mat <- matrix(NA, length(idx[[i]]), 4)
      T.tmp <- as.matrix(T[idx[[i]],])
      if(ncol(T.tmp)==1){
        T.tmp <- t(T.tmp)
      }
      rownames(T.tmp) <- idx[[i]]

      if(nrow(T.tmp)>1){
        mat[,1] <- pmax(0,G[idx[[i]],i]- apply(as.matrix(T.tmp[, colnames(T)
                                                     %wo%
                                                     column]),1,sum))
        mat[,2] <- apply(cbind(G[idx[[i]],i],
                               apply(as.matrix(T.tmp[, colnames(T)
                                                 %wo%
                                                 c(column,excluded)]),
                                     1,sum)), 1, min)
        
    mat[,3] <- apply(cbind(G[idx[[i]],i],T.tmp[,column]),1,min)
    mat[,4] <-  pmax(0,G[idx[[i]],i]- apply(as.matrix(T.tmp[,
                                                   colnames(T)[colnames(T) %in%
                                                   c(column,
                                                     excluded)]]),1,sum))
      }
      else{
        mat[,1] <- pmax(0,G[idx[[i]],i]- apply(t(T.tmp[, colnames(T)
                                                     %wo%
                                                     column]),1,sum))
        mat[,2] <- apply(cbind(G[idx[[i]],i],
                               apply(as.matrix(t(T.tmp[, colnames(T)
                                                 %wo%
                                                 c(column,excluded)])),
                                     1,sum)), 1, min)

        mat[,3] <- apply(cbind(G[idx[[i]],i],t(T.tmp[,column])),1,min)
        mat[,4] <-  pmax(0,G[idx[[i]],i]- apply(as.matrix(t(T.tmp[,
                                                   colnames(T)[colnames(T) %in%
                                                   c(column,
                                                     excluded)]])),1,sum))
      }

      bounds.out[[bound.names[count]]] <-
        cbind(mat[,1]/(mat[,1]+mat[,2]), mat[,3]/(mat[,3]+mat[,4]))
      bounds.out[[bound.names[count]]][which(mat[,3]+mat[,4]==0),]<-0
      rownames(bounds.out[[bound.names[count]]]) <- idx[[i]]
      colnames(bounds.out[[bound.names[count]]]) <- c("lower",
                                                      "upper")
      glb <- max(bounds.out[[bound.names[count]]][,"lower"])
      lub <- min(bounds.out[[bound.names[count]]][,"upper"])
    }
    if(glb<lub){
      intersection[[i]] <- c(glb, lub)
    }
    else{intersection[[i]] <- NA}

    count <- count +1
  }
  out <- list(bounds = bounds.out, intersection = intersection,
              threshold = threshold)
  class(out) <- "bounds"
  out
}
