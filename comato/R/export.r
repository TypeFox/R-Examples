#' Saving a concept map to a TGF file
#' 
#' \code{write.tgf} stores the graph underlying a conceptmap object into a file using the "Trivial Grpah Format" (TGF).
#' @param map A conceptmap object.
#' @param file The location including filename where the file should be stored.
#' @param translation If not NULL, a vector of strings of equal length as the number of concepts used in the concept map.
#' Then, the names given in this vector will be used in the file instead of the original concepts. Can be used, for example,
#' to translate the concepts into a different language for export.
#' @return -
#' @examples
#' \dontrun{
#' #Create concept map from a random graph
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' 
#' write.tgf(conceptmap(g1), "~/cmap.tgf", 
#'           translation = c("Node_1", "Node_2", "Node_3", "Node_4", "Node_5"))
#' }
#'@export
write.tgf <- function(map, file, translation=NULL) 
{
  f = file(file, open="wt")
  for (i in 1:length(map$concepts))
  {
    if (!is.null(translation))
      write(paste(i, translation[which(tolower(translation[,1]) == tolower(map$concepts[i])), 2]), f)
    else
      write(paste(i, map$concepts[i]), f)
  }
  write("#", f)
  e = get.edgelist(map$map)
  for (i in 1:dim(e)[1])
  {
    write(paste(which(map$concepts == e[i, 1]), which(map$concepts == e[i, 2])), f)
  }
  close(f)
}


#' Saving a set of concept maps to TGF files
#' 
#' \code{write.tgf.folder} stores the graphs underlying the maps of a conceptmaps object into a folder using the "Trivial Grpah Format" (TGF).
#' The function calls \code{\link{write.tgf}} for each of the maps of a conceptmaps object. The files will be named "1.tgf", "2.tgf" and so on.
#' @param maps A conceptmap object.
#' @param folder The location where the files should be stored. The folder is created, if necessary.
#' @param translation See \code{\link{write.tgf}}.
#' @return -
#' @examples
#' \dontrun{
#' #Create concept maps from three random graphs
#' require("igraph")
#' g1 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g2 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' g3 = set.vertex.attribute(erdos.renyi.game(5, 0.7, type="gnp"), "name", value=1:5)
#' 
#' #Create conceptmaps object from three conceptmap objects
#' simple_cms = conceptmaps(list(conceptmap(g1), conceptmap(g2), conceptmap(g3)))
#' 
#' write.tgf.folder(simple_cms, "~/cmaps")
#' }
#'@export 
write.tgf.folder <- function(maps, folder, translation=NULL) 
{
  dir.create(folder)
  for (m in 1:maps$count)
    write.tgf(maps$maps[[m]], file.path(folder, paste(m, ".tgf", sep="")), translation)
}
