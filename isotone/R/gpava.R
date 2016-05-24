`gpava` <-
function(z = NULL, y, weights = NULL, solver = weighted.mean, ties = "primary", p = NA)
{
# y ... response; either a single vector or a list of vectors (blocks)
# z ... predictor (1 predictor only so far, maybe extension to >1, i.e. generalized pava)
# w ... weights; either a single vector or a list if vectors (weights)
# solver ... either weighted.mean, weighted.median, weighted.fractile, or a user-specified function
# ties ... string for tie treatment: either "primary", "secondary", "tertiary".
# a,b... fractiles for weighted.fractile solver, otherwise ignored.
    
    x <- z
    w <- weights
    merger <- c
    if ((is.matrix(y)) || (is.data.frame(y))) y <- c(as.data.frame(t(y))) #generate list
    if ((is.matrix(w)) || (is.data.frame(w))) w <- c(as.data.frame(t(w))) #generate list
    
    if (any(is.na(p))) {               #additional arguments a,b
      moreargs <- NULL
    } else {
      moreargs <- list(p = p)
    }
    
    n <- length(y)
    if(is.null(x)) x <- 1:n
    if(is.null(w)) {
        w <- if(is.list(y))
            lapply(sapply(y, length), function(u) rep.int(1, u))
        else
            rep.int(1, n)
    } else if(is.list(y)) 
        w <- as.list(w)
    w1 <- w
    y1 <- y
    
 #------------- ties -----------------
    if (ties == "primary") {
      if (is.list(y)) {
          o <- order(x, mapply(solver, y, w, MoreArgs = moreargs))    #collapse y-list into val (vector)
      } else { 
        o <- order(x,y)
      }
      r <- order(o)
      y <- y[o]                                
      w <- w[o]
    }
    
    if ((ties == "secondary") || (ties == "tertiary" )) {
      if (is.list(y)) {
        yag <- tapply(y, x, function(cps) {   
                         colMeans(do.call(rbind, cps))})
        wag <- tapply(w, x, function(cps) {    
                         colSums(do.call(rbind, cps))})
      } else {
        wag <- tapply(w,x,sum)              #sum weights (within tie); collapsed
        yag <- tapply(y,x,mean)             #mean response (within tie); collapsed
      }  
      xag <- tapply(x,x,mean)               #mean predictor (within tie); collapsed
      o <- order(xag)
      r <- order(o)
      y <- yag[o]
      w <- wag[o]
    }
    
  #----------- end ties --------------    

    n <- length(y)
    inds <- as.list(seq_len(n))    
    vals <- mapply(solver, y, w, MoreArgs = moreargs)          #applies solver for each list element (e.g. list of weighted means)     
    

    ## Combine blocks i and i + 1.    
    combine <- if(is.list(y)) {           #In the repeated data case, we explicitly merge the data (and weight) lists.
        function(i) {                     #Merge the data and indices, solve, and put things back into position i, dropping position i + 1.
            j <- i + 1L
            y[[i]] <<- merger(y[[i]], y[[j]])        #append observations 
            w[[i]] <<- c(w[[i]], w[[j]])
            
            if (is.null(moreargs)) {
              vals[i] <<- solver(y[[i]], w[[i]])       #z_i; apply function for element[[i]]
            } else {
               vals[i] <<- solver(y[[i]], w[[i]], p=p)       #weighted fractile
            }
            
            inds[[i]] <<- c(inds[[i]], inds[[j]])
            keep <- seq_len(n)[-j]
            y <<- y[keep]
            w <<- w[keep]
            vals <<- vals[keep]
            inds <<- inds[keep]
            n <<- n - 1L
        }
    } else {
        function(i) {                     #In the "simple" case, merge only indices and values.
            j <- i + 1L
            inds[[i]] <<- c(inds[[i]], inds[[j]])    #append index c(i, i+1)
             if (is.null(moreargs)) {
               vals[i] <<- solver(y[inds[[i]]], w[inds[[i]]]) #compute target function (i, i+1)
            } else {
               vals[i] <<- solver(y[inds[[i]]], w[inds[[i]]], p=p) #fractile
            }
            
            keep <- seq_len(n)[-j]
            vals <<- vals[keep]                      #delete (i+1)th y-value 
            inds <<- inds[keep]                      #delete (i+1)th y-value    index                 
            n <<- n - 1L                             #decrease number of observations
        }
    }

  #--------- end combine ----------------

  #--------- PAVA iterations -----------
    i <- 1L
    repeat {
        if(i < n) {
            if((vals[i] > vals[i + 1])) {            #if decreasing
                combine(i)
                while((i > 1L) && (vals[i - 1L] > vals[i])) {  #as long as i>1 and decreasing
                    combine(i - 1L)                  #go back 1 elementstep
                    i <- i - 1L                      
                }
            }
            else
                i <- i + 1L                          #move on in data vector/index
            }
        else break
    }

 #------------- end iterations ----------------

 yfit.notie <- rep.int(vals, sapply(inds, length))
 
 if (ties == "primary") yfit <- yfit.notie[r]
 if (ties == "secondary") yfit <- as.vector(ifelse(outer(x,xag,"=="),1,0)%*%yfit.notie[r])
 if (ties == "tertiary") 
   if (!is.list(y1)) {
     yfit <- as.vector(y1 + ifelse(outer(x,xag,"=="),1,0)%*%(yfit.notie[r]-yag[o]))
   } else {
     yfit <- as.vector((mapply(solver, y1, w1, MoreArgs = moreargs)) + ifelse(outer(x,xag,"=="),1,0)%*%(yfit.notie[r]-(mapply(solver, yag[o], wag[o], MoreArgs = moreargs))))
   }
 
 result <- list(x = yfit, z = x, y = y1, w = w1, solver = solver, call = match.call(), p = p)  
 class(result) <- "gpava"
 result
}

