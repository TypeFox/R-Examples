outTeX <- function(block.obj, namesCol = NULL, file.names = NULL, captions = NULL, digits = 2, ...){

  #require("xtable")
  
  ## takes block, assignment, or diagnose object
  if(!is.null(block.obj$blocks)){ 
    block.obj <- block.obj$blocks
  }
  if(!is.null(block.obj$assg)){ 
    block.obj <- block.obj$assg
  }
  
  for(i in 1:length(block.obj)){
    tab <- block.obj[[i]]
    nm <- names(block.obj)[i]

    if(is.null(captions)){
      caption <- paste("Group ", nm, ".", sep="")
    }else{
    	caption <- captions[[i]]
    }

    ncol.tab <- ncol(tab)
    
    ## user-specified column names
    if(!is.null(namesCol)){
      names(tab) <- namesCol
    }

    if(is.null(file.names)){
    	file.name <- paste("Group", nm, ".tex", sep="")
    	lab <- paste("group.", nm, sep = "")
    }else{
    	file.name <- paste(file.names[[i]], ".tex", sep = "")
    	lab <- paste("t:", file.names[[i]], sep = "")
    }
    
    tab.tex <- xtable::xtable(tab, label = lab, caption = caption, align =
                      c(rep("c",ncol(tab)+1)), 
                      digits = rep(digits, ncol(tab)+1), ...)
                                            
    print(tab.tex, file = file.name, caption.placement = "bottom")
  }
}
