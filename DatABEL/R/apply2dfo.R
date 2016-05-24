#' applies a function to a 'databel' object
#'
#' An iterator applying a user-defined function to
#' an object of 'databel-class'.
#'
#'
#' @param dfodata 'databel' object which is iterated over
#' @param anFUN user-defined analysis function
#' @param MAR which margin to iterate over (default = 2, usually
#' these are 'columns' used to store SNP data)
#' @param procFUN function to process the output and present that as a
#' fixed-number-of-columns matrix or fixed-length vector. Can be
#' missing if standard functions listed below are used. Pre-defined
#' processors included are "process_lm_output" (can process functions
#' "lm", "glm", "coxph") and "process_simple_output" (process output
#' from "sum", "prod", "sum_not_NA" [no. non-missing obs], "sum_NA"
#' [no. missing obs.])
#' @param outclass output to ("matrix" or "databel")
#' @param outfile if output class is "databel", the generated object
#' is bound to the \code{outfile}
#' @param type if output class is "databel", what data type to use for
#' storage
#' @param transpose whether to transpose the output
#' @param ... arguments passed to the \code{anFUN}
#'
#' @return A matrix (or 'databel'-matrix) containing results of
#' applying the function
#'
#' @author Yurii Aulchenko
#' @export
#' @examples
#' a <- matrix(rnorm(50), 10, 5)
#' rownames(a) <- paste("id", 1:10, sep="")
#' colnames(a) <- paste("snp", 1:5, sep="")
#' b <- as(a, "databel")
#' apply(a, FUN="sum", MAR=2)
#' apply2dfo(SNP, dfodata=b, anFUN="sum")
#' tA <- apply2dfo(SNP, dfodata=b, anFUN="sum",
#'                 outclass="databel", outfile="tmpA")
#' tA
#' as(tA, "matrix")
#' apply2dfo(SNP, dfodata=b, anFUN="sum", transpose=FALSE)
#' tB <- apply2dfo(SNP, dfodata=b, anFUN="sum", transpose=FALSE,
#'                 outclass="databel", outfile="tmpB")
#' tB
#' as(tB, "matrix")
#'
#' sex <- 1*(runif(10)>.5)
#' trait <- rnorm(10) + sex + as(b[, 2], "vector") +
#'          as(b[, 2], "vector") * sex * 5
#' apply2dfo(trait~SNP*sex, dfodata=b, anFUN="lm")
#' tC <- apply2dfo(trait ~ SNP * sex, dfodata=b, anFUN="lm",
#'                 outclass="databel", outfile="tmpC")
#' tC
#' as(tC, "matrix")
#' apply2dfo(trait ~ SNP * sex, dfodata=b, anFUN="lm", transpose=FALSE)
#' tD <- apply2dfo(trait ~ SNP * sex, dfodata=b, anFUN="lm",
#'                 transpose=FALSE, outclass="databel", outfile="tmpD")
#' tD
#' as(tD, "matrix")
#' rm(tA, tB, tC, tD)
#' gc()
#' unlink("tmp*")
#'


apply2dfo <- function(..., dfodata, anFUN="lm", MAR=2, procFUN,
                      outclass="matrix", outfile, type="DOUBLE", transpose=TRUE)
{
    #print("AAA")
    ### should also implement a case when outfile exists?
    if (missing(procFUN))
        if (any(c("lm", "glm", "coxph") == anFUN)) {
            procFUN <- "process_lm_output"
        }
        else if (any(c("sum", "prod", "sum_not_NA", "sum_NA") ==
                     anFUN)) {
            procFUN <- "process_simple_output"
        }

    procFUN <- match.fun(procFUN)
    anFUN   <- match.fun(anFUN)

    if (missing(dfodata)) stop("dfodata should be supplied")
    if (class(dfodata) != "databel") stop("dfodata should be of calss databel")

    if (MAR == 1) {
        SNP <- as(dfodata[1, ], "matrix")
    }
    else if (MAR ==2) {
        SNP <- as(dfodata[, 1], "matrix")
    }
    else stop("MAR should be 1 or 2")

    #print("BBBB")

    tmpout <- procFUN(eval(substitute(anFUN( ... ))))

    #print(tmpout)
    #print(dimnames(tmpout))
    dimout <- dim(tmpout)
    if (is.null(dimout)) dimout <- c(1, length(tmpout))
    namout <- dimnames(tmpout)
    #print(tmpout)
    #print(dimout)

    if (outclass == "matrix") {
        if (transpose) {
            res <- matrix(ncol=dimout[2],
                          nrow=(dimout[1])*(dim(dfodata)[MAR]))
        }
        else {
            res <- matrix(nrow=dimout[2],
                          ncol=(dimout[1])*(dim(dfodata)[MAR]))
        }
    } else if (outclass == "databel")
        {
            if (missing(outfile)) {
                stop("outfile argument must be provided with outcalss=='databel'")
            }
            # no good -- fast access in wrong direction
            if (transpose)
                res <- make_empty_fvf(outfile,
                                      nvariables=dimout[2],
                                      nobservations=(dimout[1])*(dim(dfodata)[MAR]),
                                      type=type,
                                      readonly=FALSE)
            else
                res <- make_empty_fvf(outfile,
                                      nobservations=dimout[2],
                                      nvariables=(dimout[1])*(dim(dfodata)[MAR]),
                                      type=type,
                                      readonly=FALSE)
        }
    else stop("outclass must be 'matrix' or 'databel'")
    #print(c("aaA", dim(res), class(res)))

    for (i in 1:dim(dfodata)[MAR])
        {
            if (MAR == 1) {
                SNP <- as(dfodata[i, ], "vector")
            }
            else {
                SNP <- as(dfodata[, i], "vector")
            }
            #print(summary(SNP))
            cur_res <- procFUN(eval(substitute(anFUN( ... ))))
            #print(length(cur_res))
            #print((((i-1)*dimout[1]+1):(i*dimout[1])))
            #print(dim(res))
            #print(i)
            #print("jsut before...")
            if (transpose) {
                res[(((i-1)*dimout[1]+1):(i*dimout[1])), ] <- cur_res
            }
            else {
                res[, (((i-1)*dimout[1]+1):(i*dimout[1]))] <-
                    matrix(cur_res, ncol=1)
            }
        }
        #print(c("aaa", dim(res)))
        #print(colnames(tmpout))
        #	if (MAR == 1) oMAR <- 2; else oMAR <- 1;

    if (transpose) {
        if (!is.null(colnames(tmpout))) {
            nms <- list(dimnames(res)[[1]], colnames(tmpout))
            #print("i-go-go")
            #print(nms)
            dimnames(res) <- nms
        }
    } else {
        if (dimout[1]==1) {
            #print("o-ho-ho")
            #print(MAR)
            #print(dim(res))
            #print(dim(dfodata))
            #print(dimnames(dfodata)[[2]])
            if (!is.null(dimnames(dfodata)[[MAR]]))
                dimnames(res) <- list(dimnames(res)[[1]], dimnames(dfodata)[[MAR]])
        } else {
            #print("a-ha-ha")
            dimnames(res) <- list(dimnames(res)[[1]],
                                  paste(as.vector(t(
                                      matrix(rep(dimnames(dfodata)[[MAR]], dimout[1]), ncol=dimout[1])
                                      )),
                                        dimnames(tmpout)[[1]], sep="_")
                                  )
        }
    }
    #print("set dimnames[[2]]")


    #print(c("aaa", dim(res)))
    if (transpose) {
        #print("here")
        #print(dimout[1])
        #print(dimnames(dfodata)[[2]])
        #print(dim(res))
        if (dimout[1]==1) {
            if (!is.null(dimnames(dfodata)[[MAR]])) {
                nms <- list(dimnames(dfodata)[[MAR]], dimnames(res)[[2]])
                #print("i-go-go-0")
                #print(nms)
                dimnames(res) <- nms
            }
        } else {
            #print("there")
            #print(dimnames(dfodata)[[2]])
            #print(dimout[1])
            #print(dimnames(tmpout)[[1]])
            if (!is.null(dimnames(dfodata)[[MAR]])) {
                nms <- list(
                    paste(as.vector(t(
                        matrix(rep(dimnames(dfodata)[[MAR]], dimout[1]), ncol=dimout[1])
                        )),
                          dimnames(tmpout)[[1]], sep="_"),
                    dimnames(res)[[2]])
                #print("i-go-go-1")
                #print(nms)
                dimnames(res) <- nms
            }
        }
    } else {
        if (!is.null(colnames(tmpout))) {
            dimnames(res) <- list(colnames(tmpout), dimnames(res)[[2]])
        }
    }
    #print("set dimnames[[1]]")
    return(res)
}
