#' Importing a concept map from a GraphML file.
#' 
#' \code{read.yEd} reads a graphML file that was created by yEd and imports the graph as a conceptmap object.
#' @param file The filename and path that should be read.
#' @param strip Passed to the call of \code{\link{conceptmap.igraph}} that is used to create the conceptmap object.
#' @return A conceptmap object.
#' @examples
#' \dontrun{
#' #Assuming that the data is in "~/cmap.graphml"
#' cm = read.yEd("~/cmap.graphml")
#' }
#'@export 
read.yEd <- function(file, strip=TRUE) 
{
	res <- xmlTreeParse(file, useInternalNodes=TRUE)
	g <- xmlChildren(xmlRoot(res))$graph
	
	nodes <- g["node"]
	node.list <- sapply(nodes, xmlAttrs)
	node.labels <- unlist(sapply(nodes, xmlValue))
	
	edges <- g["edge"]
	t <- sapply(edges, xmlAttrs)
	data <- cbind(t[2,], t[3,], unlist(sapply(edges, xmlValue)))
	dimnames(data) <- NULL
	for (i in 1:dim(data)[1]) 
	{
		data[i, 1] = node.labels[which(node.list == data[i, 1])]
		data[i, 2] = node.labels[which(node.list == data[i, 2])]
	}
	return(conceptmap(data, strip=strip))
}


#' Importing a concept map from a TGF file.
#' 
#' \code{read.tgf} reads a TGF file and imports the graph as a conceptmap object.
#' @param file The filename and path that should be read.
#' @param strip Passed to the call of \code{\link{conceptmap.igraph}} that is used to create the conceptmap object.
#' @return A conceptmap object.
#' @examples
#' \dontrun{
#' #Assuming that the data is in "~/cmap.tgf"
#' cm = read.tgf("~/cmap.tgf")
#' }
#'@export
read.tgf <- function(file, strip=TRUE) 
{
	input <- readLines(file)
	empty.concepts = grep(pattern="^[0-9]+ $", input)
  if (length(empty.concepts) > 0)
  input <- input[-empty.concepts]
	delim <- which(input == "#")
  if (delim == 1)
    return(conceptmap(NULL))
	l <- unlist(sapply(input[1:(delim-1)], FUN <- function(str) {strsplit(sub(pattern=" ", replacement="\1", x=str), split="\1")}))
	nodes <- l[(1:(delim-1))*2]
  ids <- l[(1:(delim-1))*2 - 1]
	data <- c()
	if ( delim+1 < length(input))
	{
	  for (i in (delim+1):length(input)) 
	  {
  		line <- strsplit(input[i], " ")[[1]]
	  	v1 <- nodes[which(ids == as.numeric(line[1]))]
  		v2 <- nodes[which(ids == as.numeric(line[2]))]
		  e <- ""
		  for (i in 3:length(line))
  			e <- paste(e, line[i])
		  data <- rbind(data, c(v1, v2, trim(e)))
	  }
  }
	cm <- conceptmap(data)
  if (!strip)
    cm <- modify.concepts(cm, nodes)
	return(cm)
}


#' Importing a set of concept maps from GraphML files.
#' 
#' \code{read.folder.yEd} reads several graphML files that were created by yEd and imports them as a conceptmaps object.
#' @param folder The path of a folder in which every graphML file is read.
#' @param strip Passed to the call of \code{\link{read.yEd}} that is used to import the single concept maps.
#' @return A list consisting of a conceptmaps object and the list of filenames (in the same order as the maps in the conceptmaps object).
#' @examples
#' \dontrun{
#' #Assuming that the data is in the folder "~/cmaps"
#' cm = read.folder.yEd("~/cmaps")
#' }
#'@export
read.folder.yEd <- function(folder, strip=TRUE) 
{
	d <- getwd()
	setwd(folder)
	toDo <- list.files(pattern="*.graphml")
	order = c()
	res <- list()
	for (f in toDo) {
	  cat(basename(f))
	  cat("\n")
    m <- read.yEd(f, strip)
		res <- c(res, list(m))
		order = c(order, strsplit(basename(f), "\\.")[[1]][1])
	}
	setwd(d)
	return(list(conceptmaps(res, filter=F), order))
}


#' Importing a set of concept maps from TGF files.
#' 
#' \code{read.folder.tgf} reads several TGF files and imports them as a conceptmaps object.
#' @param folder The path of a folder in which every TGF file is read.
#' @param strip Passed to the call of \code{\link{read.tgf}} that is used to import the single concept maps.
#' @return A list consisting of a conceptmaps object and the list of filenames (in the same order as the maps in the conceptmaps object).
#' @examples
#' \dontrun{
#' #Assuming that the data is in the folder "~/cmaps"
#' cm = read.folder.tgf("~/cmaps")
#' }
#'@export
read.folder.tgf <- function(folder, strip=TRUE) 
{
  toDo <- list.files(path=folder, pattern="*.tgf", full.names=T)
  order = c()
  res <- list()
  for (f in toDo) {
    cat(basename(f))
    cat("\n")
    m <- read.tgf(f, strip)
    res <- c(res, list(m))
    order = c(order, strsplit(basename(f), "\\.")[[1]][1])
  }
  return(list(conceptmaps(res, filter=F), order))
}