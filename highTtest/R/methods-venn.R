if(!isGeneric("vennD")){
  setGeneric(name="vennD", 
             def=function(x, gamma, ...){standardGeneric("vennD")})
}

setMethod("vennD",
  signature = c(x="highTtest", gamma="numeric"),
  definition = function(x, gamma, ...){

    if( !requireNamespace("colorfulVennPlot", quietly=TRUE) ) {
      stop("R package colorfulVennPlot is requireded for this function.", 
           call.=FALSE)
    }

    tst <- x@gammas - gamma
    igamma <- which(tst > -1e-8 & tst < 1e-8)

    if(length(igamma) == 0) {
      cat("Requested gamma value not included in object provided.\n")
      return(NULL)
    }

    if(is.null(x@BH) && is.null(x@ST)){
      cat("A Venn Diagram method is not available for 1 set.\n")
      return(NULL)
    } else if(!is.null(x@BH) && !is.null(x@ST)){
      venn3(x, igamma, ...)
    } else {
      venn2(x, igamma, ...)
    }
  }
)

venn2 <- function(x, igamma,...){

  if(!is.null(x@BH)){
    area <- cbind(x@CK[,igamma], x@BH[,igamma])
    category <- c("CK","BH")
  } else if(!is.null(x@ST)){
    area <- cbind(x@CK[,igamma], x@ST[,igamma])
    category <- c("CK","ST")
  }

  same <- sum(rowSums(area)==2)
  area <- colSums(area)

  cgy <- c(paste(category[1], " (", area[1], ")", sep=""),
           paste(category[2], " (", area[2], ")", sep=""))

  area[1] <- area[1] - same
  area[2] <- area[2] - same

  vec <- c(area, same)
  names(vec) <- c("10","01","11")

  args <- list(...)

  if(is.null(args$Title)){
    args$Title <- paste(paste(category,collapse=" & "),
                        "at level", x@gammas[igamma], sep=" ")
  }

  if(is.null(args$Colors)){
    args$Colors <- c("red","yellow")
  }

  args$labels <- cgy
  args$reverseLabelOrdering <- FALSE
  args$x <- vec

  plot.new()

  do.call(colorfulVennPlot::plotVenn2d,args)

}


venn3 <- function(x, igamma, ...){

  area <- cbind(x@CK[,igamma], x@BH[,igamma], x@ST[,igamma])

  same12 <- sum(rowSums(area[,c(1,2)])==2 & area[,3]==0)
  same13 <- sum(rowSums(area[,c(1,3)])==2 & area[,2]==0)
  same23 <- sum(rowSums(area[,c(2,3)])==2 & area[,1]==0)

  same123 <- sum(rowSums(area) == 3)

  area <- colSums(area)

  category <- c(paste("CK(", area[1], ")", sep=""),
                paste("BH(", area[2], ")", sep=""),
                paste("ST(", area[3], ")", sep=""))

  area[1] <- area[1] - same12 - same13 - same123
  area[2] <- area[2] - same12 - same23 - same123
  area[3] <- area[3] - same13 - same23 - same123

  vec <- c(area, same12, same23, same13, same123)
  names(vec) <- c("100","010","001","110","011","101","111")

  args <- list(...)
  if(is.null(args$Title)){
    args$Title <- paste("CK, BH, & ST at level", x@gammas[igamma], sep=" ")
  }
  if(is.null(args$Colors)){
    args$Colors <- c("red","yellow","orange")
  }

  args$labels <- category
  args$x <- vec

  plot.new()

  do.call(colorfulVennPlot::plotVenn3d,args)


}


