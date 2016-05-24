

#' Read a phylo4d object from files
#'
#' This function creates an object of class \code{phylo4d} by combining
#' a phylogenetic tree and its associated tips data stored in two distinct files.
#' 
#' @param phylo.file the name of the file which the phylogenetic tree is to be read from.
#' @param data.file the name of the file which the tips data are to be read from.
#' @param phylo.format the format of the phylogenetic tree provided.
#' Possible formats are "\code{newick}" and "\code{nexus}".
#' @param data.format the format of the table with tips data.
#' Possible formats are "\code{table}", "\code{csv}", "\code{csv2}", "\code{delim}" and "\code{delim2}".
#' 
#' @details \code{phylo.file} and \code{data.file} can be provided as objects of mode character
#' or double-quoted strings.
#' 
#' The phylogenetic tree can be imported in two formats.
#' \itemize{
#'    \item The \code{newick} format refers
#' to the simple parenthetic format known as the Newick or New Hampshire format. The tree is
#' read by calling the function \code{\link[ape]{read.tree}} of the \pkg{ape} package.
#'    \item The \code{nexus} format refers to NEXUS format. The tree is
#' read by calling the function \code{\link[ape]{read.nexus}} of the \pkg{ape} package.
#' }
#' 
#' Tips data are imported from a table formatted file. The different formats allow to use
#' different separator and decimal characters.
#' They correspond to the variants of \code{\link[utils]{read.table}}:
#' \itemize{
#'   \item \code{table} use \code{\link[utils]{read.table}} with default settings.
#'   \item \code{csv} use \code{\link[utils]{read.csv}} with default settings.
#'   \item \code{csv2} use \code{\link[utils]{read.csv2}} with default settings.
#'   \item \code{delim} use \code{\link[utils]{read.delim}} with default settings.
#'   \item \code{delim2} use \code{\link[utils]{read.delim2}} with default settings.
#' }
#' 
#' @return An object of class \code{phylo4d}.
#' 
#' @seealso \code{\link[phylobase]{phylo4d}} to create a \code{phylo4d} object.
#'
#' @export
read.p4d <- function(phylo.file, data.file, phylo.format = "newick", data.format = "table"){
  
  phylo.format <- match.arg(phylo.format, c("newick", "nexus"))
  data.format <- match.arg(data.format, c("table", "csv", "csv2", "delim", "delim2"))
  
  if(phylo.format == "newick"){
    phy <- read.tree(file = phylo.file)
  }
  if(phylo.format == "nexus"){
    phy <- read.nexus(file = phylo.file)
  }

  read.fun <- paste("read.", data.format, "(file = \"", data.file, "\", header = TRUE, row.names = 1)", sep = "")
  dat <- eval(parse(text = read.fun))
  
  p4d <- phylo4d(x = phy, tip.data = dat)
  return(p4d)
}