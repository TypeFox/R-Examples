#' Obtain An SSE Sequence Vector From A PDB Object
#'
#' Results are similar to that returned by stride(pdb)$sse and dssp(pdb)$sse. 
#'
#' @details call for its effects.
#' 
#' @param pdb an object of class \code{pdb} as obtained from
#'    function \code{\link{read.pdb}}. 
#' @param verbose logical, if TRUE warnings and other messages will be printed.
#'
#' @return a character vector indicating SSE elements for each amino acide residue.
#'    The 'names' attribute of the vector contains 'resno', 'chain', 'insert', and
#'    'SSE segment number', seperated by the character '_'.
#'
#' @seealso \code{\link{dssp}}, \code{\link{stride}}, \code{\link{bounds.sse}}
#'
#' @author Barry Grant & Xin-Qiu Yao
#' 
#' @examples
#' \donttest{
#'    pdb <- read.pdb("1a7l")
#'    sse <- pdb2sse(pdb)
#'    sse
#' } 
pdb2sse <- function(pdb, verbose = TRUE) {
  ##- Function to obtain an SSE sequence vector from a PDB object
  ##   Result similar to that returned by stride(pdb)$sse and dssp(pdb)$sse
  ##   This could be incorporated into read.pdb() if found to be more generally useful
  ##

  if(is.null(pdb$helix) & is.null(pdb$sheet)) {
    if(verbose) 
       warning("No helix and sheet defined in input 'sse' PDB object: try using dssp()")
    ##ss <- try(dssp(pdb)$sse)
    ## Probably best to get user to do this separately due to possible 'exefile' problems etc..
    return(NULL)
  }

  ## An empty full length SSE vector
  ref <- pdb$atom[pdb$calpha, c("resno", "chain", "insert")]
  ref <- paste(ref$resno, ref$chain, ref$insert, sep="_")
  ss <- rep(" ", length(ref))
  names(ss) <- ref
 
  ## loop over 'Helix' and 'Sheet'
  symbol <- c(helix="H", sheet="E")
  for(i in names(symbol)) {

     sse <- pdb[[i]]

     if(length(sse$start) > 0) {
        for(j in 1:length(sse$start)) {
           chain <- ifelse(sse$chain[j]=="", NA, sse$chain[j])
           insert0 <- ifelse(names(sse$start[j])=="", NA, names(sse$start[j]))
           sse.ref0 <- paste(sse$start[j], chain, insert0, sep = "_")

           insert1 <- ifelse(names(sse$end[j])=="", NA, names(sse$end[j]))
           sse.ref1 <- paste(sse$end[j], chain, insert1, sep = "_")
        
           ii <- match(sse.ref0, ref); jj <- match(sse.ref1, ref)
           if(any(is.na(c(ii, jj)))) {
              if(verbose) 
                  warning(paste("The", i, "No.", j, 
                     "start/end with non-protein residue."))
           } else {
              inds <- seq(ii, jj)
              ss[inds] <- symbol[i]
              names(ss)[inds] <- paste(names(ss)[inds], "_", j, sep="")
           }
        }
     } 
  }
  return(ss)
}
