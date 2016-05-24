make_by <- function(term, object, data)
{
  if(!missing(data) && !is.character(data)) {
    by <- data[[object$by]]
    if(is.factor(by) && nlevels(by) > 1) {
      nocenter <- paste(c("lasso", "nigmix", "ridge", "ra"), ".smooth.spec", sep = "")
      term <- paste(paste(rmf(object$by), rmf(levels(by)), sep = ""), "*", term, sep = "")
      if((k <- length(term)) > 1)
        for(j in 1:k) 
          if(!grepl("center", term[j]) && is.null(object$xt$center))
            if(!class(object) %in% nocenter)
              term[j] <- gsub(")", ",center)", term[j])
    } else term <- paste(rmf(object$by), "*", term, sep = "")
  } else term <- paste(rmf(object$by), "*", term, sep = "")

  return(term)
}

