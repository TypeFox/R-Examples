#' Nonrandom character sequence with sample
#' 
#' Find the seed necessary to produce a character sequence by using sample
#' 
#' @return \code{\link{cat}}s command into the console that can be copypasted to anyone's R script.
#' @note nameSample may take a lot of time, due to nchar^26 possibilities. That's why it warns about strings longer than 5 characters
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, April 2014
#' @seealso \code{\link{yearSample}} to wish a happy new year, \code{\link{set.seed}}, \code{\link{sample}}, \code{\link{letters}}
#' @keywords character
#' @export
#' @examples
#' 
#' ## Not run in RCMD check as they're very time consuming
#' \dontrun{
#' nameSample("berry")  # After that, you can send the result to colleagues:
#' # Kind regards from
#' set.seed(1248272); paste(sample(letters,5,TRUE), collapse='')
#' 
#' # calculation time
#'                                        # on my slow laptop:   # on PC
#' system.time(nameSample("berr"))        # 25 s # berry: 57 s     10  23
#' system.time(nameSample("berr", FALSE)) # 23 s          53 s      9  20
#' 
#' # let <- sapply(1:4, function(n) apply(replicate(n, letters[sample(15)]), 1, paste, collapse=""))
#' # calctime <- sapply(let, function(x) system.time(nameSample(x, progress=F))[3])
#' # write.table(calctime, "calctime_nameSample.txt")
#' ctfile  <- system.file("extdata/calctime_nameSample.txt",  package="berryFunctions")
#' ctfile2 <- system.file("extdata/calctime_nameSample2.txt", package="berryFunctions")
#' calctime <- read.table(ctfile)
#' # regression result in hours:
#' expReg(nchar(rownames(calctime))-8, calctime[,1], xlim=c(1,7), ylim=c(-3,4),
#'        predict=7)/3600
#' 
#' # For my 3 times faster computer:
#' calctime <- read.table(ctfile2)
#' expReg(nchar(rownames(calctime))-8, calctime[,1], xlim=c(1,7), ylim=c(-3,4),
#'        predict=c(4,7))/c(1,3600)
#' # 4 sec for 4 letters are expected to be 10 hours for 7 letters...
#' 
#' }
#' 
#' @param name Character string. long strings (>>5) will compute a VERY long time!
#' @param progress Logical. Monitor progress by printing a dot every 10000 tries? DEFAULT: TRUE for long names (nchar(name)>3).
#' @param estimatetime Estimate computation time? DEFAULT: nc>4
#' @param continue Continue without asking? DEFAULT:  FALSE
#'
nameSample <- function(
name,
progress=FALSE,
estimatetime=nc>4,
continue=FALSE
)
{
nc <- nchar(name)
name <- tolower(name)
if(any(!strsplit(name,"")[[1]] %in% letters)) stop("'", name, "' contains characters not available in letters.")
if(nc>3) {message("function may take a bit of time."); flush.console() }
# input length checking:
if(estimatetime)
  {
  estcalctime <- 10^(1.36*nc-4.66)/3600
  unit <- " hours"
  if(estcalctime <1) { estcalctime <- estcalctime*60 ; unit=" minutes"}
  answer <- readline(paste0("Estimated time: ", signif(estcalctime, 2), unit,
  ". Do you want a more exact estimate (takes 5-20 s)? y/n : "))
  if(answer == "y")
    {
    randomwords <- function(n) paste(sample(letters,n), collapse="")
    let <- sapply(rep(4:2, c(3,15,50)), randomwords)
    if(requireNamespace("pbapply", quietly=TRUE)) sapply <- pbapply::pbsapply
    suppressMessages(calctime <- sapply(let, function(x) system.time(nameSample(x))[3]))
    nchar_name <- nchar(names(calctime))-8
    estcalctime <- expReg(nchar_name, calctime, xlim=c(1, nc),
                    ylim=c(-3, 1.36*nc-4.66), ylab="calculation time seconds", predictnew=nc)/3600
    unit <- " hours"
    if(estcalctime[1] <1) { estcalctime <- estcalctime*60 ; unit=" minutes"}
    if(estcalctime[1] <1) { estcalctime <- estcalctime*60 ; unit=" seconds"}
    }
  if(continue) answer <- "y" else
     answer <- readline(paste0("Estimated time: ", signif(estcalctime[2], 2), " to ",
          signif(estcalctime[3], 2), unit, ". Do you want to continue? y/n : "))
  if(answer == "n") stop("Function cancelled by user.")
  }
#
# The actual work:
anf <- Sys.time()
seed_is_false <- function(i)
  {
  set.seed(i)
  if(progress) if(i%%10000==0) {cat("."); flush.console()}
  if(progress) if(i%%10000==0)  if((i/10000)%%options()[["width"]]==0) cat("\n")
  paste( sample(letters, nc, replace=TRUE), collapse="") != name
  }
i <- 1
while( seed_is_false(i) ) i <- i+1
output <- paste0("set.seed(", i, "); paste(sample(letters,", nc,
                 ",T), collapse='')\n")
message(if(progress) "\n", "Computation time was: ", round(difftime(Sys.time(),anf, units="secs")/60,2), " minutes.")
message(output)
return(invisible(output))
}
