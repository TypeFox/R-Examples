.SelectOrderData <- function(data, fct, which.lbs, which.Order){
   ## for data to be plot in performs two selections:
   ## on unordered (original) data (acc. to which.lbs)
   ## on data ordered acc. to fct a selection acc. to which.Order is done
   ## return value: list with elements
   #      data, the selected/thinned out data,
   #      y = fct(data)
   #      ind the indices of the selected data in the original data
   #      ind1 the indices of the data selected by which.lbs in the original data
     dimL <- !is.null(dim(data))
     d1  <- if(dimL) dim(data) else 1
     n   <- if(dimL) nrow(data) else length(data)
     ind <- 1:n
     
     ### selection
     if(is.null(which.lbs)) which.lbs <- 1:n
     which.lbs0 <- (1:n) %in% which.lbs
     n <- sum(which.lbs0)
     which.lbx <- rep(which.lbs0, length.out=length(data))
     data <- data[which.lbx]
     if(dimL) dim(data) <- c(n,d1[-1])
     ind <- ind[which.lbs0]
     ### function evaluation
     y <- if(dimL) apply(data, 1, fct) else sapply(data,fct)
     ## ordering
     oN <- order(y)
     ind1 <- rev(ind[oN])
     
     ## selection of ordered
     if(is.null(which.Order))
          which.Order <- 1:n
     oN <-  oN[(n+1)-which.Order]
     data <- if(dimL) data[oN,] else data[oN]
     y <- y[oN]
     ind <- ind[oN]

     return(list(data=data, y=y, ind=ind, ind1=ind1))
}


if(FALSE){
x <- rnorm(1000)
.SelectOrderData(x, function(x1)x1^2, 1:100, 1:5)
}