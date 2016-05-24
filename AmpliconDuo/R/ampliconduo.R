ampliconduo <-
function(A , B = NULL, sample.names = NULL, correction = "fdr", ...){
  
  ###internal functions####
  make.duo <-
    function(A, B, sample.name, correction, ...){
      cat(".")
      in.both <- which(A!=0 | B!=0) # amplicons with both 0 are removed
      duo <- data.frame(A, B)
      duo <- duo[in.both,] ## keep row names from original data, get with row.names(duo)
      A.all <- sum(duo[,1])
      B.all <- sum(duo[,2])
      
      tst <- 
        t(
          apply(X=duo[, 1:2], 
                FUN = 
                  function(x) test.AiBi(x, A.all, B.all, ...), 
                MARGIN = 1)
        )
      
      q <- p.adjust(tst[,1], correction)
      
      duo <- cbind(
        duo,
        tst,
        q
      )
      
      duo[[8]] <- rep(FALSE, times = length(A))
      duo[[9]] <- (rep(sample.name, times = length(A)))
      names(duo) <- c("freqA", "freqB", "p", "OR", "CI.low", "CI.up", "q", "rejected",  "sample") 
      return(duo)
    }
  
  test.AiBi <-
    function(x, n.A, n.B, conf.level = 0.95, or = 1, alternative = "two.sided") 
    {
      m <- 
        matrix(
          c(x[1], 
            x[2],
            n.A-x[1],
            n.B-x[2]),
          nrow=2
        )
      t <- fisher.test(m, conf.level = conf.level, or = or, alternative = alternative)
      #gives back p-value, odds-ratio, conf.interval:
      c(t$p.value, 
        t$estimate, 
        t$conf.int) 
    }
  #############################
  
  data <- list()
  if((class(A) == "data.frame" & class(B) == "data.frame") || (class(A)=="list" & class(B)=="list") ){
    
    if(length(A)!= length(B)){ stop("Dimensions of A and B are different.") }
    
    if(class(A)=="list"){
      A <- as.data.frame(A)
      B <- as.data.frame(B)
    }
    
    if(length(sample.names)!= length(A)){
      if(!is.null(sample.names)){
        warning("Count of names does not correspond to sample count. Numbering will be used instead.")
      }  
      sample.names <- c(1:length(A))
    }else{
      if(class(sample.names) == "list"){
        sample.names = unlist(sample.names)
      }
    }
    for (i in 1:length(A)){
      duo <- make.duo(A[i], B[i], sample.names[i], correction, ...) 
      data[[i]] <- duo
      names(data)[i]<- sample.names[i]
    }
  }
  
  else{## data are all in A
    if((class(A)=="list" | class(A)=="data.frame")){
      if(as.integer(length(A))%%2 != 0){stop("Odd number of columns in A.")}
      A <- as.data.frame(A)
      sample.number = length(A) / 2
      if(length(sample.names)!= sample.number){
        if(!is.null(sample.names)){
          warning("Count of names does not correspond to sample count. Numbering will be used instead.")
        } 
        sample.names <- c(1:length(A))
      }
      
      for (i in 1:sample.number){
        duo <- make.duo(A[(i*2-1)], A[(i*2)], sample.names[i], correction, ...)
        data[[i]] <- duo
        names(data)[i]<- sample.names[i]
      }
    }
    else(stop("wrong data format."))
  }
  return(data)
}
