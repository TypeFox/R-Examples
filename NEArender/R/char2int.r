#' char2int function
#' 
#' This function reads input list objects and maps unique numbers to each member in the network/ genesets. 
#' @seealso \code{\link{nea.render}}
#'
#' 
#' @param net.list  the global network object, pre-created using \code{\link{import.net}}
#' @param gs.list.1 AGS or FGS object that lists members of each individual AGS/FGS, pre-created using \code{\link{samples2ags}}
#' @param gs.list.2 FGS or AGS object that lists members of each individual FGS/AGS, pre-created using \code{\link{samples2ags}}
#' @keywords internal



char2int  <- function (net.list, gs.list.1, gs.list.2 = NULL) {
  
  all.names <- unique(c(names(net.list$links), unlist(net.list$links), unlist(gs.list.1)));
  map.names <- 1:length(all.names);
  names(map.names) <- all.names;
  mapped <- NULL; mapped$net <- NULL; mapped$gs <- NULL;
  for (n in names(net.list$links)) {
    mapped$net[[map.names[n]]] <- as.list(map.names[net.list$links[[n]]]);
  }
  for (i in c("a", "b")) {
    if (i == "a") {gsl = gs.list.1;}
    if (i == "b") { 
      if (!is.null(gs.list.2)) {
        gsl = gs.list.2;
      } else {
        gsl = NULL;
      }}
    if (!is.null(gsl)) {
      for (n in names(gsl)) {
        mapped$gs[[i]][[n]] <- as.list(map.names[gsl[[n]]]);
      }}
  }
  return(mapped);
}