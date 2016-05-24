is.int <- function(x,n=NULL,debug=FALSE) {
  opt_in <- (!is.null(data_type) && !is.null(n) &&
             n %in% names(data_type)[data_type=="integer"])
  opt_out <- (!is.null(data_type) && !is.null(n) &&
              n %in% names(data_type)[data_type!="integer"])
  trunc_test <- all(trunc(x)==x)
  mode_test <- storage.mode(x)=="integer"
  if (debug) cat("opt_in",opt_in,"opt_out",opt_out,"trunc_test",trunc_test,"mode_test",mode_test,"\n")
  opt_in || (!opt_out &&
             !intcheck=="none" &&
             ((intcheck=="strict" && mode_test) ||
              (intcheck=="trunc" && trunc_test)))
}

intcheck <- "strict"
data_type <- NULL
i1 <- 1
stopifnot(!is.int(i1))
i2 <- 1L
stopifnot(is.int(i2))

intcheck <- "trunc"
stopifnot(is.int(i1))

data_type <- c(i1="numeric")
stopifnot(!is.int(i1,"i1"))
