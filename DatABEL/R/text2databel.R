#' converts text file to filevector format
#'
#' The file provides the data to be converted to filevector format.
#' The file may provide the data only (no row and column names)
#' in which case col/row names may be left empty or provided in
#' separate files (in which case it is assumed that names are provided
#' only for the imported columns/rows -- see skip-options).
#' There is an option to skip a number of
#' first ros and columns. The row and column names may also
#' be provided in the file itself, in which case one needs to
#' tell the row/column number providing column/row names.
#' Unless option "R_matrix" is set to TRUE, it is asumed that the
#' number of columns is always the same acorss the file. If above
#' option is provided, it is assumed that both column and row
#' names are provided in the file, and the first line contains
#' one column less than other lines (such is the case with
#' files produced from R using the function
#' \code{write.table(...,col.names=TRUE,row.names=TRUE)}.
#'
#' @param infile input text file name
#' @param outfile output filevector file name; if missing,
#' it is set to infile+".filevector"
#' @param colnames where are the column names stored? If missing,
#' no column names; if integer, this denotes the row of the
#' input file where the column names are specified; if
#' character string then the string specifies the name of the file
#' with column names
#' @param rownames where are the row names stored? If missing,
#' no row names; if integer, this denotes the column of the
#' input file where the row names are specified; if
#' character string then the string specifies the name of the file
#' with row names
#' @param skipcols how many columns of the input file to skip
#' @param skiprows how many rows of the input file to skip
#' @param transpose whether the file is to be transposed
#' @param  R_matrix if true, the file format is assumed to follow
#' the format of R data matrix produced with
#' \code{write.table(...,col.names=TRUE,row.names=TRUE)}
#' @param type data DatABEL type to use ("DOUBLE", "FLOAT", "INT",
#' "UNSIGNED_INT", "UNSIGNED_SHORT_INT", "SHORT_INT", "CHAR", "UNSIGNED_CHAR")
#' @param cachesizeMb cache size for the resulting 'databel-class' object
#' @param readonly whether the resulting 'databel-class' object should
#' be opened in readonly mode
#' @param naString the string used for missing data (default: NA)
#' @param unlinkTmpTransposeFiles Boolean to indicate whether
#' the intermediate "_fvtmp.fvi/d" files should be deleted. Default:
#' TRUE. These intermediate files are generated while transposing the
#' filevector files.
#'
#' @author  Yurii Aulchenko
#'
#' @return The converted file is stored in the file system, a
#' \link{databel-class} object connection to the file is returned.
#' @export
#' @examples
#'
#' cat("this is an example which you can run if you can write to the
#' file system\n")
#'
#' \dontrun{
#'
#' # create matrix
#' NC <- 5
#' NR <- 10
#' data <- matrix(rnorm(NC*NR),ncol=NC,nrow=NR)
#' rownames(data) <- paste("r",1:NR,sep="")
#' colnames(data) <- paste("c",1:NC,sep="")
#' data
#'
#' # create text files
#' write.table(data, file="test_matrix_dimnames.dat", row.names=TRUE,
#'             col.names=TRUE, quote=FALSE)
#' write.table(data, file="test_matrix_colnames.dat", row.names=FALSE,
#'             col.names=TRUE, quote=FALSE)
#' write.table(data, file="test_matrix_rownames.dat", row.names=TRUE,
#'             col.names=FALSE, quote=FALSE)
#' write.table(data, file="test_matrix_NOnames.dat", row.names=FALSE,
#'             col.names=FALSE, quote=FALSE)
#' write(colnames(data), file="test_matrix.colnames")
#' write(rownames(data), file="test_matrix.rownames")
#'
#' # generate identical data
#' text2databel(infile="test_matrix_dimnames.dat",
#'              outfile="test_matrix_dimnames", R_matrix=TRUE)
#' x <- databel("test_matrix_dimnames")
#' data <- as(x, "matrix")
#' data
#'
#' # convert text two filevector format
#'
#' text2databel(infile="test_matrix_NOnames.dat",
#'              outfile="test_matrix_NOnames.fvf",
#'              colnames="test_matrix.colnames",
#'              rownames="test_matrix.rownames")
#' x <- databel("test_matrix_NOnames.fvf")
#' if (!identical(data, as(x, "matrix"))) stop("not identical data")
#'
#' text2databel(infile="test_matrix_NOnames.dat",
#'              outfile="test_matrix_NOnames_T.fvf",
#'              colnames="test_matrix.colnames",
#'              rownames="test_matrix.rownames", transpose=TRUE)
#' x <- databel("test_matrix_NOnames_T.fvf")
#' if (!identical(data, t(as(x, "matrix")))) stop("not identical data")
#'
#' text2databel(infile="test_matrix_rownames.dat",
#'              outfile="test_matrix_rownames.fvf",
#'              rownames=1, colnames="test_matrix.colnames")
#' x <- databel("test_matrix_rownames.fvf")
#' if (!identical(data, as(x, "matrix"))) stop("not identical data")
#'
#' text2databel(infile="test_matrix_colnames.dat",
#'              outfile="test_matrix_colnames.fvf",
#'              colnames=1, rownames="test_matrix.rownames")
#' x <- databel("test_matrix_colnames.fvf")
#' if (!identical(data, as(x, "matrix"))) stop("not identical data")
#'
#' text2databel(infile="test_matrix_dimnames.dat",
#'              outfile="test_matrix_dimnames.fvf", R_matrix=TRUE)
#' x <- databel("test_matrix_dimnames.fvf")
#' if (!identical(data, as(x, "matrix"))) stop("not identical data")
#'
#' # stupid extended matrix in non-R format
#' newmat <- matrix(-100, ncol=NC+3, nr=NR+2)
#' newmat[3:(NR+2), 4:(NC+3)] <- data
#' newmat[2, 4:(NC+3)] <- paste("c", 1:NC, sep="")
#' newmat[3:(NR+2), 3] <- paste("r", 1:NR, sep="")
#' newmat
#' write.table(newmat, file="test_matrix_strange.dat",
#'             col.names=FALSE, row.names=FALSE, quote=FALSE)
#'
#' text2databel(infile="test_matrix_strange.dat",
#'              outfile="test_matrix_strange.fvf",
#'              colnames=2, rownames=3)
#' x <- databel("test_matrix_strange.fvf")
#' if (!identical(data, as(x, "matrix"))) stop("not identical data")
#'
#' }
#'

text2databel <- function(infile, outfile, colnames, rownames,
                         skipcols, skiprows, transpose=FALSE,
                         R_matrix=FALSE, type="DOUBLE",
                         cachesizeMb=64, readonly=TRUE,
                         naString="NA", unlinkTmpTransposeFiles=TRUE)
{
    if (missing(infile)) stop("input file name (infile) is required")
    if (missing(outfile))  outfile <- paste(infile, ".filevector", sep="")
    if (missing(skipcols)) skipcols <- 0
    if (missing(skiprows)) skiprows <- 0

    if (missing(colnames)) {
        cnrow <- 0
        cnfile <- ""
    } else {
        if (is.character(colnames)) {
            cnrow <- 0
            cnfile <- colnames
        }
        else {
            cnrow <- as.integer(colnames)
            cnfile <- ""
        }
    }

    if (missing(rownames)) {
        rncol <- 0
        rnfile <- ""
    } else {
        if (is.character(rownames)) {
            rncol <- 0
            rnfile <- rownames
        }
        else {
            rncol <- as.integer(rownames)
            rnfile <- ""
        }
    }

    if (R_matrix) cnrow <- rncol <- 1
    if (rncol > 0 && skipcols == 0) skipcols <- rncol
    if (cnrow > 0 && skiprows == 0) skiprows <- cnrow

    intype <- filevector_type(type)

    charnames <- as.character(c(infile, outfile, rnfile, cnfile,  naString))
    intnames <- as.integer(c(rncol, cnrow, skiprows, skipcols,
                             transpose, R_matrix, intype))
    #print(charnames)
    #print(intnames)
    tmp <- .Call("text2fvf_R", charnames, intnames, package="DatABEL")

    if (unlinkTmpTransposeFiles==TRUE) {
        unlink(paste0(outfile, "_fvtmp.fvi"))
        unlink(paste0(outfile, "_fvtmp.fvd"))
    }

    return(databel(outfile, cachesizeMb=cachesizeMb, readonly=readonly))
}
