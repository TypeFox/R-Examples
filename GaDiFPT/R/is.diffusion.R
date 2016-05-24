is.diffusion <-
  function (obj) 
{
    if (inherits(obj, "diffusion") & is.list(obj) & (length(obj) == 2)) 
      if (all(lapply(obj, length) == 1) & 
            all(unlist(lapply(obj, is.character)))) 
            return(identical(names(obj), c("mean","var")))                                         
    return(FALSE)
  }
