#' write rnk to a file from matrix W.
#' 
#' write a rnk file from matrix W in a returned object of function \code{DNMF}.
#' The rnk format is referred \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{RNK}
#' 
#' @param object a returned object of function \code{DNMF}
#' @param fn the output filename. Default is "./tmp.rnk"
#' @param type type o2m (Default) or o2o. to compare with multi sample labels. 
#' o2m means one Vs others, while o2o means one Vs another one.
#' @export
#' @examples
#' \dontrun{
#' rnk(dnmf_result, fn="tmp.rnk")
#' }
rnk <- function (object, fn="./tmp.rnk", type="o2m"){
    
    type <- match.arg(type, c("o2m", "o2o"))
    r <- object$r
    if (r==2) {
        write.table(object$rnk, fn, sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
    } else if (type=="o2o"){
        for (i in 1:(r-1)) {
            for (j in (i+1):r) {
                write.table(object$W[,j]-object$W[,i], paste0(dirname(fn), "/", j, "_", i, "_", basename(fn)), sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
            }
        }
    } else if (type=="o2m") {
        for (i in 1:r) {
            write.table(object$W[,i]-rowMeans(object$W[,-i]), paste0(dirname(fn), "/", i, "_", basename(fn)), sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
        }
    }
}


