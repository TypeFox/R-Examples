#' Texture Matrix Visualization
#'
#' 
#'@param x Matrix of class "glcm", "glrlm", "glszm" or "mglszm"
#'@param xlab The label for the x-axis
#'@param ylab The label for the y-axis
#'@param col Use viridis scale if available
#'@examples
#' \dontrun{
#' image(psf)
#' image(glszm(psf))
#'}
#'@name image.radiomics
NULL
#> NULL

colscale <- function(n) {
  if (requireNamespace("viridis", quietly = TRUE)) {
    viridis::viridis(n)
  } else {
    heat.colors(n)
  }
}

#'GLCM image
#'@rdname image.radiomics
#'@export
setMethod("image",
          signature(x="glcm"),
          definition = function(x, xlab="Grey Level", ylab="Grey Level",
                                col=colscale(length(unique(c(x@.Data))))
                                ) {
            image(x=c(1:nrow(x@.Data)),
                  y=c(1:ncol(x@.Data)),
                  z=x@.Data,
                  xlab=xlab, ylab=ylab,
                  col=col,
                  axes=FALSE
                  )
            axis(1, at=c(1:nrow(x@.Data)), labels=rownames(x@.Data))
            axis(2, at=c(1:ncol(x@.Data)), labels=colnames(x@.Data))
            
          }
)

#'GLRLM image
#'@rdname image.radiomics
#'
#'@export
setMethod("image",
          signature(x="glrlm"),
          definition = function(x, xlab="Grey Level", ylab="Run Length",
                                col=colscale(length(unique(c(x@.Data))))
          ) {
            image(x=c(1:nrow(x@.Data)),
                  y=c(1:ncol(x@.Data)),
                  z=x@.Data,
                  xlab=xlab, ylab=ylab,
                  col=col,
                  axes=FALSE
            )
            axis(1, at=c(1:nrow(x@.Data)), labels=rownames(x@.Data))
            axis(2, at=c(1:ncol(x@.Data)), labels=colnames(x@.Data))
            
          }
)

#'GLSZM image
#'@rdname image.radiomics
#'
#'@export
setMethod("image",
          signature(x="glszm"),
          definition = function(x, xlab="Grey Level", ylab="Zone Size",
                                col=colscale(length(unique(c(x@.Data))))
          ) {
            image(x=c(1:nrow(x@.Data)),
                  y=c(1:ncol(x@.Data)),
                  z=x@.Data,
                  xlab=xlab, ylab=ylab,
                  col=col,
                  axes=FALSE
            )
            axis(1, at=c(1:nrow(x@.Data)), labels=rownames(x@.Data))
            axis(2, at=c(1:ncol(x@.Data)), labels=colnames(x@.Data))
          }
)

#'MGLSZM image
#'@rdname image.radiomics
#'
#'@export
setMethod("image",
          signature(x="mglszm"),
          definition = function(x, xlab="Grey Level", ylab="Zone Size",
                                col=colscale(length(unique(c(x@.Data))))
          ) {
            image(x=c(1:nrow(x@.Data)),
                  y=c(1:ncol(x@.Data)),
                  z=x@.Data,
                  xlab=xlab, ylab=ylab,
                  col=col,
                  axes=FALSE
            )
            axis(1, at=c(1:nrow(x@.Data)), labels=rownames(x@.Data))
            axis(2, at=c(1:ncol(x@.Data)), labels=colnames(x@.Data))
          }
)


