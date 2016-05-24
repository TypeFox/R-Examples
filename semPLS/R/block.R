#' Blocks of manifest variables belonging to latent
#'
##' @param model the specified model as returned by read.splsm
##' @return a named list containing the block of manifest variables for each latent additionally there is the attribute 'mode' telling the measurement model: 'A' for reflective and 'B' for formative

block <-
function(latent, manifest, measuremod){
  ln <- length(latent)        # number of latent variables
  colnames(measuremod) <- NULL
  blocks <- list()
  
  for (i in 1:ln){
    blocks[[i]] <- measuremod[c(which(measuremod[,1]==latent[i], which(measuremod[,2]==latent[i]))),]
    blocks[[i]] <- append(blocks[[i]],
                          measuremod[c(which(measuremod[,2]==latent[i], which(measuremod[,1]==latent[i]))),])
    blocks[[i]] <- sort(blocks[[i]][blocks[[i]] %in% manifest])
        
    # determine the mode ("A"=reflective, "B"=formative)
    if (all(blocks[[i]] %in% measuremod[,2])){
      attr(blocks[[i]], "mode") <- "A"
    }
    else if (all(blocks[[i]] %in% measuremod[,1])) {
      attr(blocks[[i]], "mode") <- "B"
    }
    else stop("A block must be either formative or reflective, not both!")
  }
  names(blocks) <- latent
  
  return(blocks)
}
