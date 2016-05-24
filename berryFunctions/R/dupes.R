#' Duplicate lines in file
#' 
#' Number of duplicates per line of (text) file. Per default saved to file
#' which can be loaded into excel / libreoffice. With conditional formatting of
#' the first column, colors show for each line how often it occurs in the file.
#' A LibreOffice file is included. Note: OpenOffice does not provide color
#' scales based on cell values.
#' 
#' @return Either: a data.frame with line numbers of duplicate rows and the number of duplicates\cr 
#'         Or: a file is written with the number of duplicates and the original \code{file} content.
#' @note This has not been tested al that much - feedback is heavily welcome!
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2014
#' @seealso \code{\link{compareFiles}}
#' @keywords IO file character
#' @export
#' @examples
#' 
#' file <- system.file("extdata/doublelines.txt", package="berryFunctions")
#' dupes(file, tofile=FALSE)
#' dupes(file, tofile=FALSE, ignore.empty=TRUE)  
#' 
#' ## These are skipped by rcmd check (opening external places is not allowed):
#' \dontrun{dupes(file)}
#' 
#' # a template file (dupes.ods) for libreOffice Calc is available here:
#' system.file("extdata", package="berryFunctions")
#' 
#' \dontrun{system2("nautilus", system.file("extdata/dupes.ods", package="berryFunctions"))}
#' 
#' # To open folders with system2:
#' # "nautilus" on linux ubuntu
#' # "open" or "dolphin" on mac
#' # "explorer" or "start" on windows
#' 
#' @param file File name (character string)
#' @param ignore.empty Should empty lines be ignored? DEFAULT: TRUE
#' @param ignore.space Should leading/trailing whitespace be ignored? DEFAULT: TRUE
#' @param tofile Logical: should output be directed to a file? 
#'        Otherwise, a dataframe with line numbers and number of duplicates of that line 
#'        will be printed in the console. DEFAULT: missing(n)
#' @param n Show only the first n values if \code{tofile=FALSE}. DEFAULT: length(d)
#' 
dupes <- function(
file,
ignore.empty=TRUE,
ignore.space=TRUE,
tofile=missing(n),
n=length(d)
)
{
# empty lines or lines with only (up to 9) spaces:
spaces <- if(ignore.empty) sapply(0:9, function(i) 
                  paste(rep(" ", i), collapse="")) else FALSE
R <- readLines(file)
R2 <- if(ignore.space) removeSpace(R) else R
# indices of duplicated:
d <- which(duplicated(R2, incomparables=spaces, fromLast=TRUE) |
           duplicated(R2, incomparables=spaces, fromLast=FALSE)  )
# return n entries to console if tofile==F:
if(!tofile) {
nd <- sapply(d, function(i) sum(R2[i]==R2[-i]))
return(head(data.frame(line=d, number=nd), n))
}
# else:
nd <- sapply(1:length(R2), function(i) sum(R2[i]==R2[-i]))
if(ignore.empty) nd[R2==""] <- ""
write.table(data.frame(nd, R), paste(file,"_dupes.txt"), row.names=F, 
            col.names=F, quote=F, sep="\t")
message("Created the file '", file,"_dupes.txt'\nin getwd: ", getwd())
}
