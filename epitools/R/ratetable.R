"ratetable" <-
  function(..., byrow = FALSE,
           rev = c("neither", "rows", "columns", "both")
           ){
    lx <- list(...)
    if(length(lx)==0)
      {
        stop("No arguments provided")
      }
    if(length(lx)==1 && (is.character(lx[[1]]) || is.factor(lx[[1]])))
      {
        stop("Must be numeric vector or matrix.")
      }    
    ## r x 2 table
    if(length(lx)==1 && is.matrix(lx[[1]])
       && nrow(lx[[1]])>=2 && ncol(lx[[1]])==2)
      {
        x <- lx[[1]]
        if(is.null(dimnames(lx[[1]])))
          {
            nr <- nrow(x)
            rn <- paste("Exposed", 1:nr, sep="")
            cn <- c("Count", "Person-time")
            dimnames(x) <- list(Predictor = rn, Outcome = cn)      
          }
      }
    ## 2 vectors
    if(length(lx)==2 && is.vector(lx[[1]]) && is.vector(lx[[2]]))
      {      
        x <- cbind(lx[[1]], lx[[2]])
        if(!is.null(names(lx)))
          {
            colnames(x) <- names(lx)
          }
        if(is.null(rownames(x)) && !is.null(colnames(x)))
          {
            nr <- nrow(x)
            rn <- paste("Exposed", 1:nr, sep="")
            rownames(x) <- rn      
          }
        if(is.null(dimnames(x)))
          {
            nr <- nrow(x)
            rn <- paste("Exposed", 1:nr, sep="")
            cn <- c("Count", "Person-time")
            dimnames(x) <- list(Predictor = rn, Outcome = cn)      
          }
      }
    ## >=4 numbers
    is.even <- function(x){ifelse(x%%2==0, TRUE, FALSE)}
    if(length(lx)>=4 && all(sapply(list(1,2,3,4,5),is.numeric))
       && is.even(length(lx)) && all(sapply(lx,length)==1))
      {
        x <- matrix(sapply(lx,as.vector), ncol = 2, byrow = byrow)
        nr <- nrow(x)
        rn <- paste("Exposed", 1:nr, sep="")
        cn <- c("Cases", "Person-time")
        dimnames(x) <- list(Predictor = rn, Outcome = cn)
      }
    ## 1 vector
    if(length(lx)==1 && is.vector(lx[[1]])
       && is.numeric(lx[[1]]) && is.even(length(lx[[1]])))
      {
        x <- matrix(lx[[1]], ncol = 2, byrow = byrow)
        nr <- nrow(x)
        rn <- paste("Exposed", 1:nr, sep="")
        cn <- c("Cases", "Person-time")
        dimnames(x) <- list(Predictor = rn, Outcome = cn)
      }
    nrx <- nrow(x)
    ncx <- ncol(x)
    reverse <- match.arg(rev)
    if(reverse=="rows") finalx <- x[nrx:1,]
    if(reverse=="columns") finalx <- x[,ncx:1]
    if(reverse=="both") finalx <- x[nrx:1,ncx:1]
    if(reverse=="neither") finalx <- x  
    finalx
  }
