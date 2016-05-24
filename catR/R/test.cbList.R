test.cbList<-function (list, cbGroup) 
{
    if (is.null(cbGroup)) 
        res <- list(test = FALSE, message = "Item bank has no subgroup for content balancing!")
    else {
        if (!is.list(list)) 
            res <- list(test = FALSE, message = paste("'", deparse(substitute(list)), 
                "' is not a list", sep = ""))
        else {
            if (length(list) != 2) 
                res <- list(test = FALSE, message = paste("'", 
                  deparse(substitute(list)), "' must have exactly two elements", 
                  sep = ""))
            else {
                if (names(list)[1] != "names") 
                  res <- list(test = FALSE, message = paste("first element of '", 
                    deparse(substitute(list)), "' must be named 'names'", 
                    sep = ""))
                else {
                  if (names(list)[2] != "props") 
                    res <- list(test = FALSE, message = paste("second element of '", 
                      deparse(substitute(list)), "' must be named 'props'", 
                      sep = ""))
                  else {
                    if (length(list$names) != length(list$props)) 
                      res <- list(test = FALSE, message = paste("Elements 'names' and 'props' of '", 
                        deparse(substitute(list)), "' must have the same length", 
                        sep = ""))
                    else {
                      if (length(list$names) != length(unique(cbGroup))) 
                        res <- list(test = FALSE, message = paste("Number of subgroups in '", 
                          deparse(substitute(list)), "' and '", 
                          deparse(substitute(cbGroup)), "' are different", 
                          sep = ""))
                      else {
                        nr <- 0
                        for (i in 1:length(list$names)) nr <- nr + 
                          length(unique(cbGroup)[unique(cbGroup) == 
                            list$names[i]])
                        if (nr != length(unique(cbGroup))) 
                          res <- list(test = FALSE, message = paste("Mismatch in names of subgroups in '", 
                            deparse(substitute(list)), "' and '", 
                            deparse(substitute(itemBank)), "'", 
                            sep = ""))
                        else {
                          if (!is.numeric(list$props)) 
                            res <- list(test = FALSE, message = paste("Element 'props' in '", 
                              deparse(substitute(list)), "' is not numeric", 
                              sep = ""))
                          else {
                            if (min(list$props) < 0) 
                              res <- list(test = FALSE, message = paste("All components of 'props' in '", 
                                deparse(substitute(list)), "' must be positive", 
                                sep = ""))
                            else {
                              res <- list(test = TRUE, message = "ok")
                            }
                          }
                        }
                      }
                    }
                  }
                }
            }
        }
    }
    return(res)
}
