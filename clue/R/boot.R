cl_boot <-
function(x, B, k = NULL, 
         algorithm = if(is.null(k)) "hclust" else "kmeans",
         parameters = list(), resample = FALSE)
{
    clusterings <- if(!resample) {
        x <- rep.int(list(x), B)        
        eval(as.call(c(list(as.name("lapply"), x, algorithm),
                       if(!is.null(k)) list(k),
                       parameters)))
    }
    else {
        replicate(B,
                  expr = {
                      algorithm <- match.fun(algorithm)
                      ## <NOTE>
                      ## This is not quite perfect.  We have
                      ## cl_predict() to encapsulate the process of
                      ## assigning objects to classes, but for sampling
                      ## from the objects we assume that they correspond
                      ## to the *rows* of 'x'.  Argh.
                      ## </NOTE>
                      ind <- sample(NROW(x), replace = TRUE)
                      train <- if(length(dim(x)) == 2) x[ind, ] else x[ind]
                      out <- eval(as.call(c(list(algorithm, train),
                                            if(!is.null(k)) list(k),
                                            parameters)))
                      as.cl_partition(cl_predict(out, x, "memberships"))
                  },
                  simplify = FALSE)
    }

    cl_ensemble(list = clusterings)
}
