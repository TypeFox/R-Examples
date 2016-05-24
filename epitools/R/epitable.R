"epitable" <-
  function(...,
           ncol = 2, byrow = TRUE,
           rev = c("neither", "rows", "columns", "both")){
    lx <- list(...)
    if(length(lx)==0){
      stop("No arguments provided")
    }
    if(length(lx)==1 && (is.character(lx[[1]]) || is.factor(lx[[1]]))){
      stop("Single factor or character vector not allowed.")
    }
    ## r x c table
    if(length(lx)==1 && is.matrix(lx[[1]]) &&
       nrow(lx[[1]])>=2 && ncol(lx[[1]])>=2)
      {
        x <- lx[[1]] 
        if(is.null(dimnames(lx[[1]])))
          {
            nr <- nrow(x)
            nc <- ncol(x)
            rn <- paste("Exposed", 1:nr, sep="")
            cn <- paste("Disease", 1:nc, sep="")
            dimnames(x) <- list(Predictor = rn, Outcome = cn)      
          }
      }
    ## 2 vectors
    if(length(lx)==2 &&
       (is.vector(lx[[1]]) || is.factor(lx[[1]])) &&
       (is.vector(lx[[2]]) || is.factor(lx[[2]]))
       )
      {
        x <- table(lx[[1]], lx[[2]]) 
        if(nrow(x)<2 || ncol(x)<2)
          {
            stop("must have 2 or more rows and columns")
          }
        if(is.null(names(lx)))
          {
            names(dimnames(x)) <- c("Predictor", "Outcome")
          } else names(dimnames(x)) <- names(lx)
      }
    ## >=4 numbers
    is.even <- function(x){ifelse(x%%2==0, TRUE, FALSE)}
    if(length(lx)>=4 && all(sapply(list(1,2,3,4,5),is.numeric))
       && is.even(length(lx)) && all(sapply(lx,length)==1)) 
      {
        x <- matrix(sapply(lx,as.vector), ncol = ncol, byrow = byrow)
        nr <- nrow(x)
        nc <- ncol(x)
        rn <- paste("Exposed", 1:nr, sep="")
        cn <- paste("Disease", 1:nc, sep="")
        dimnames(x) <- list(Predictor = rn, Outcome = cn)      
      }
    ## 1 vector
    if(length(lx)==1 && is.vector(lx[[1]]) &&
       is.numeric(lx[[1]]) && is.even(length(lx[[1]]))
       )
      {
        x <- matrix(lx[[1]], ncol = ncol, byrow = byrow)
        nr <- nrow(x)
        nc <- ncol(x)
        rn <- paste("Exposed", 1:nr, sep="")
        cn <- paste("Disease", 1:nc, sep="")
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
