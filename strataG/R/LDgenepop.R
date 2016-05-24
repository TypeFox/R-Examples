#' @title Linkage Disequlibrium 
#' @description Calculate linkage disequilibrium p-values using GENEPOP.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param show.output logical. Show GENEPOP output on console?
#' @param delete.files logical. Delete GENEPOP input and output files when done?
#' @param label character string to use to label GENEPOP input and output files.
#' @param ... other arguments to be passed to \code{\link{genepop}}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
LDgenepop <- function(g, show.output = FALSE, delete.files = TRUE, 
                      label = "linkage.genepop", ...) {
  
  # Run Genepop
  g <- stratify(g, rep("1", nInd(g)))
  output <- genepop(g, output.ext = ".DIS", show.output = show.output, 
                    label = label, other.settings = "MenuOptions=2.1")
  if(!is.list(output)) return(NULL)
  
  result <- scan(output$files["output.fname"], what = "character", quiet = TRUE)
  loc.names <- output$locus.names
  
  # Create empty matrix
  numrows <- ((length(loc.names) ^ 2) - length(loc.names)) / 2
  result.mat <- matrix(as.character(NA), numrows, 5)
  
  # Find starting points
  loc <- grep("Switches", result, value = F) + 7
  first.col <- grep(names(loc.names)[1], result)[1]
  num.skip <- first.col - loc
  row.mask <- c(rep(F, num.skip), rep(T, 5))
  
  # Read matrix
  for(r in 1:numrows) {
    result.mat[r, ] <- result[loc:(loc + num.skip + 5)][row.mask]
    loc <- loc + num.skip + 5
  }
  
  # Convert to data.frame and format columns
  result.df <- data.frame(result.mat, stringsAsFactors = FALSE)
  colnames(result.df) <- c("Locus.1", "Locus.2", "p.value", "std.err", "switches")
  result.df$p.value <- as.numeric(result.df$p.value)
  result.df$std.err <- as.numeric(result.df$std.err)
  result.df$switches <- as.integer(result.df$switches)
  result.df$Locus.1 <- loc.names[result.df$Locus.1]
  result.df$Locus.2 <- loc.names[result.df$Locus.2]    
  
  if(delete.files) for(f in output$files) if(file.exists(f)) file.remove(f) 
  result.df
}