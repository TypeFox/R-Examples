make.fdt.format.classes <- function (x,
                                     right,
                                     pattern)
{
  tmp <- strsplit(x,
                  ',')

  res <- lapply(tmp,
                function (.vals) {
                  .vals[1L] <- sprintf(ifelse(right,
                                              paste('(',
                                                    pattern,
                                                    sep=''),
                                              paste('[', 
                                                    pattern,
                                                    sep='')),
                                       as.numeric(substring(.vals[1L], 2)))
                  .vals[2L] <- sprintf(ifelse (right,
                                               paste(pattern,
                                                     ']',
                                                     sep=''),
                                               paste(pattern,
                                                     ')',
                                                     sep='')),
                                       as.numeric(substring(.vals[2L], 1,
                                                            nchar(.vals[2L]) - 1)))
                  paste(.vals, 
                        collapse=', ')
                })

  invisible(unlist(res))
}

