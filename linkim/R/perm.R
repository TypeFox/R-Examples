perm <-
function (n, r, v = NULL,repetition=NULL,...) {
    if (is.null(v))v = 1:n 
    if (is.null(repetition))repetition=TRUE
    if(length(v)!= n)stop("Numbers of Elements of v Should be n")
    if (r == 1) matrix(v, n, 1) 
    else if (n == 1) matrix(v, 1, r) 
    else { 
      X <- NULL
      if(repetition)for (i in 1:n) X <- rbind(X, cbind(v[i], perm(n, r - 1, v, repetition)))
      else for (i in 1:n) X <- rbind(X, cbind(v[i], perm(n - 1, r - 1, v[-i], repetition))) 
      X 
    } 
 }
