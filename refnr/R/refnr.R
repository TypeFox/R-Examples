# Copyright 2016 Kim. All right reserved.

#------------------------------------------------------------------------------
# refnr provides a simple functionality that refines a set of indicators
#
# TODO(kim.seonghyun) : add more cases
#
# Author: kim.seonghyun@scipi.net (Kim Seonghyun)

#------------------------------------------------------------------------------
refnr <- function(.data, formulas) {
  stopifnot(names(formulas) == c("Name", "Formula"))

  res <- data.frame(matrix(vector(), nrow(.data), 0),
                     stringsAsFactors=F)

  for (i in 1:nrow(formulas)) {
    tryCatch({
      refined <- eval(parse(text = as.character(formulas[i, "Formula"])),
                      envir = .data)
      res[as.character(formulas[i, "Name"])] <- refined
    }, error = function(e) {
      message(e)
    })
  }

  return(res)
}
