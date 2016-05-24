"timealign" <- function(origpos, newpos, how, error.how, matchtol = 0) {
    .Call("time_align", as(origpos, "timeDate"), as(newpos, "timeDate"),
        as(c(how, error.how), "character"), as(matchtol, "numeric")+0)
}
"numalign" <- function(origpos, newpos, how, error.how, matchtol = 0) {
    .Call("num_align",  as(origpos, "numeric")+0, as(newpos, "numeric")+0,
          as(c(how, error.how), "character"), as(matchtol, "numeric")+0)
}
