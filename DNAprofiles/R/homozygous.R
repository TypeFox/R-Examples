#' For which markers is the profile homozygous/heterozygous?
#' 
#' @aliases heterozygous
#' @param x An integer matrix specifying either a single profile or a number of profiles.
#' @param markers Character vector stating the markers to check. Defaults to all markers of \code{x}.
#' @return logical matrix The value at column \code{m} in row \code{i} denotes whether profile \code{i} is homozygous/heterozygous for marker \code{m}.
#' @examples
#' x <- t(c(1L,1L,1L,2L))
#' colnames(x) <- c("locus1.1","locus1.2","locus2.1","locus2.2")
#'
#' homozygous(x,markers = "locus1") # TRUE
#' homozygous(x,markers = "locus2") # FALSE
#' homozygous(x) # t(c(TRUE,FALSE))

#' y <- t(c(1L,NA,1L,2L))
#' colnames(y) <- c("locus1.1","locus1.2","locus2.1","locus2.2")
#' homozygous(y,markers = "locus1") # NA
#' homozygous(y,markers = "locus2") # FALSE
#' homozygous(y)   # t(c(NA,FALSE))
#' heterozygous(y) # t(c(NA,TRUE))
#' @export
homozygous <- function(x,markers=get.markers(x)){
  x.markers <- get.markers(x) # does a check on the column names of x as well
  if (!all(markers %in% x.markers)){      
    stop("x does not contain marker(s) ",paste(markers[!markers %in% x.markers],collapse=", "))}  
  tmp <- sapply(markers,function(m) unname(x[,paste(m,1,sep=".")]==x[,paste(m,2,sep=".")]),simplify = FALSE)  
  do.call(cbind,tmp)
}
#' @export
heterozygous <- function(x,markers=get.markers(x)){
  !homozygous(x = x,markers = markers)
}