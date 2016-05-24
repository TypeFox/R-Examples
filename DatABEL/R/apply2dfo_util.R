#' 'apply2dfo'-associated functions
#'
#' A number of functions used in conjunction with 'apply2dfo'.
#' Standardly supported apply2dfo's anFUN analysis functions include
#' 'lm', 'glm', 'coxph', 'sum', 'prod', "sum_not_NA" (no. non-missing
#' obs), and "sum_NA" (no. missing obs.).  Pre-defined processing
#' functions include "process_lm_output" (can process functions "lm",
#' "glm", "coxph") and "process_simple_output" (process output from
#' "sum", "prod", "sum_not_NA", "sum_NA")
#'
#' @aliases process_lm_output process_simple_output sum_not_NA sum_NA
#'
#' @usage process_lm_output(lmo, verbosity=2)
#'
#' @param lmo object returned by analysis with "lm", "glm", etc.
#' @param verbosity verbosity
#'
#' @seealso \link{apply2dfo}
#' @export
#'
#' @examples
#' a <- matrix(rnorm(50),10,5)
#' rownames(a) <- paste("id",1:10,sep="")
#' colnames(a) <- paste("snp",1:5,sep="")
#' b <- as(a,"databel")
#' apply(a,FUN="sum",MAR=2)
#' apply2dfo(SNP,dfodata=b,anFUN="sum",procFUN="process_simple_output")
#' apply2dfo(SNP,dfodata=b,anFUN="sum",transpose=FALSE)
#'
#' sex <- 1*(runif(10)>.5)
#' trait <- rnorm(10)+sex+as(b[,2],"vector")+as(b[,2],"vector")*sex*5
#' apply2dfo(trait~SNP*sex,dfodata=b,anFUN="lm",procFUN="process_lm_output")
#'

process_lm_output <- function(lmo, verbosity=2)
{
    if (class(lmo) != "lm" && class(lmo) != "glm") {
        stop(paste("cannot process object of type", class(lmo)))
    }
    if (length(grep("coef", names(lmo))) != 1) {
        stop("weird lmo object")
    }
    nams <- names(lmo$coef)
    #print(lmo)
    lmo <- summary(lmo)
    if (length(grep("coef", names(lmo))) != 1) {
        stop("weird lmo object")
    }
    lmo <- lmo$coef
    #print(lmo)
    #print(nams)
    if (dim(lmo)[1]<length(nams)) {
        lmo <- rbind(lmo, matrix(NA,
                                 ncol=dim(lmo)[2],
                                 nrow=length(nams)-dim(lmo)[1]))
        rownames(lmo) <- nams
    }
    snprows <- grep("SNP", rownames(lmo))
    if (length(snprows) < 1) {
        snprows <- rep(NA, 10)
    }

    if (verbosity <= 0) {
        selcols <- c(1)
    }
    else if (verbosity == 1) {
        selcols <- c(1, 2)
    }
    else if (verbosity == 2) {
        selcols <- c(1, 2, 4)
    }
    else {
        selcols <- c(1:dim(lmo)[2])
    }

    out <- matrix(lmo[snprows, selcols], ncol=length(selcols))
    #print(dimnames(lmo))
    #print(c(snprows))
    #print(c(selcols))
    #print(dim(out))
    #print(out)
    dimnames(out) <- list(dimnames(lmo)[[1]][snprows],
                          dimnames(lmo)[[2]][selcols])
    #print(dimnames(out))

    return(out)
}

#' @export
process_simple_output <- function(o) return(o)

#' @export
sum_not_NA <- function(x) return(sum(!is.na(x)))

#' @export
sum_NA <- function(x) return(sum(is.na(x)))
