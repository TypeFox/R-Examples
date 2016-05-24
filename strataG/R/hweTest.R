#' @title Hardy-Weinberg Equilibrium
#' @description Calculate Hardy-Weinberg equilibrium p-values.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param use.genepop logical. Use GENEPOP to calculate HWE p-values? 
#'   If \code{FALSE} then \code{\link[pegas]{hw.test}} is used.
#' @param B the number of replicates for the Monte Carlo procedure for 
#'   \code{\link[pegas]{hw.test}}.
#' @param show.output logical. Show output from GENEPOP?
#' @param delete.files logical. Delete GENEPOP files when done?
#' @param label character string to use to label GENEPOP files.
#' @param ... arguments to be passed to \code{\link{genepop}}.
#'  
#' @note If \code{use.genepop = TRUE}, the command line version of GENEPOP v.4 
#'   must be properly installed and available on the command line, so it is 
#'   executable from any directory. On PC's, this requires having it in a 
#'   folder in the PATH environmental variable. On Macs, the executable should 
#'   be installed in a folder like \code{/usr/local/bin}.
#'   
#' @return a vector of p-values for each locus.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{genepop}}, \code{\link[pegas]{hw.test}}
#' 
#' @aliases hwe HWE
#' @importFrom pegas hw.test
#' @export
#' 
hweTest <- function(g, use.genepop = FALSE, B = 1000, show.output = FALSE, 
                    delete.files = TRUE, label = "HWE.genepop", ...) {
  
  if(use.genepop) {
    # Run Genepop  
    g <- stratify(g, rep("1", nInd(g)))
    output <- genepop(g, output.ext = ".D", show.output = show.output, label = label, 
                      other.settings = c("HWtests=MCMC", "MenuOptions=1.1"), ...
    )
    
    result <- scan(output$files["output.fname"], what = "character", quiet = TRUE)
    locus.names <- output$locus.names
    
    # Find HWE p-value
    result[result == "-"] <- NA
    result[result == "No"] <- NA
    i <- which(result %in% names(locus.names))
    hwe.p <- as.numeric(result[i + 1])
    names(hwe.p) <- locus.names[result[i]]
    
    if(delete.files) for(f in output$files) if(file.exists(f)) file.remove(f)
    hwe.p
  } else {
    hwe <- hw.test(gtypes2genind(g), B = B)
    hwe[, ncol(hwe)]
  }
}