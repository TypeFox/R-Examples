#' Create a PCA movie
#'
#' Make a basic movie using movie3d()
#'
#' A wrapper around movie3d with some default settings. Used to create a
#' movie on the fly.
#' @param ... Any parameters will be passed to movie3d()
#' @return The value returned by movie3d()
#' @examples
#' \dontrun{
#' data( metabo )
#' pca <- prcomp( metabo[,-1], scale.= TRUE )
#' pca3d(pca, group=metabo[,1])
#' makeMoviePCA()
#'}
#' @export
makeMoviePCA <- function(...) {
  rot <- spin3d( axis= c( 0, 1, 0 ) )
  movie3d( rot, duration= 12, ... )
}


#' Save a 3D PCA snapshot
#'
#' Take a snapshot of the 3D PCA to a file.
#' 
#' This is just a wrapper around rgl.snapshot.
#' @param file Name of the file to save the snapshot to
#' data( metabo )
#' pca <- prcomp( metabo[,-1], scale.= TRUE )
#' pca3d(pca, group=metabo[,1])
#' snapshotPCA3d("testfile.png")
#' @export
snapshotPCA3d <- function(file) {
  rgl.snapshot(filename=file, fmt="png", top=TRUE)
}
