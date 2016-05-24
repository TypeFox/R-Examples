##' ascii table generator
##'
##' @author David Hajage
##' @export
asciiCbind <- setRefClass("asciiCbind",
                          fields = c("args",
                            "caption",
                            "caption.level",
                            "frame",
                            "grid",
                            "col.width",
                            "width"),

                          methods = list(

                            test.nrow = function(args = .self$args) {

                              nrow.args <- sapply(args, function(x) {
                                if (!is.null(x$tgroup)) {
                                  if (!is.list(x$tgroup))
                                    x$tgroup <- list(x$tgroup)
                                }
                                if (!is.null(x$bgroup)) {
                                  if (!is.list(x$bgroup))
                                    x$bgroup <- list(x$bgroup)
                                }
                                nrow(x$x) + length(x$tgroup) + length(x$bgroup) + x$include.colnames
                              })

                              if (length(unique(nrow.args)) != 1)
                                stop("x and y must have same number of rows", call. = FALSE)
                            },

                            show.asciidoc = function(args = .self$args, caption = .self$caption, caption.level = .self$caption.level, frame = .self$frame, grid = .self$grid, col.width = .self$col.width, width = .self$width) {
                              test.nrow(args)

                              args.output <- lapply(args, function(x) {
                                x$caption <- NULL
                                x$caption.level <- NULL
                                x$frame <- NULL
                                x$grid <- NULL
                                x$col.width <- 1
                                x$width <- 0
                                capture.output(x$show.asciidoc())
                              })

                              args.output[-length(args.output)] <- lapply(args.output[-length(args.output)], function(x) {
                                x[1] <- sub(" $", "", x[1])
                                x[length(x)] <- sub(" $", "", x[length(x)])
                                x
                              })

                              args.output[-1] <- lapply(args.output[-1], function(x) {
                                x[1] <- sub("^\\|", "=", x[1])
                                x[length(x)] <- sub("^\\|", "=", x[length(x)])
                                x[c(-1, -length(x))] <- paste("", x[c(-1, -length(x))])
                                x
                              })

                              cat(header.asciidoc(caption = caption, caption.level = caption.level, frame = frame, grid = grid, col.width = col.width, width = width))
                              cat(do.call("paste", c(args.output, list(sep = ""))), sep = "\n")
                            },

                            show.rest = function(args = .self$args, caption = .self$caption, caption.level = .self$caption.level, frame = .self$frame, grid = .self$grid, col.width = .self$col.width, width = .self$width) {
                              test.nrow(args)

                              args.output <- lapply(args, function(x) {
                                x$caption <- NULL
                                x$caption.level <- NULL
                                capture.output(x$show.rest())[-1]
                              })

                              args.output[-1] <- lapply(args.output[-1], function(x) {
                                x <- sub("^\\||\\+", "", x)
                                x
                              })

                              cat(header.rest(caption = caption, caption.level = caption.level), sep = "\n")
                              cat(do.call("paste", c(args.output, list(sep = ""))), sep = "\n")
                            },

                            show.org = function(args = .self$args, caption = .self$caption, caption.level = .self$caption.level, frame = .self$frame, grid = .self$grid, col.width = .self$col.width, width = .self$width) {
                              test.nrow(args)

                              args.output <- lapply(args, function(x) {
                                x$caption <- NULL
                                x$caption.level <- NULL
                                capture.output(x$show.org())
                              })

                              args.output[-length(args.output)] <- lapply(args.output[-length(args.output)], function(x) {
                                x <- sub("\\|$", "", x)
                                x
                              })

                              args.output[-1] <- lapply(args.output[-1], function(x) {
                                x <- sub("^(\\|)(-+)", "\\+\\2", x)
                                x
                              })

                              cat(header.org(caption = caption, caption.level = caption.level), sep = "\n")
                              cat(do.call("paste", c(args.output, list(sep = ""))), sep = "\n")
                            },

                            show.t2t = function(args = .self$args, caption = .self$caption, caption.level = .self$caption.level, frame = .self$frame, grid = .self$grid, col.width = .self$col.width, width = .self$width) {
                              test.nrow(args)

                              args.output <- args
                              args.output[-length(args.output)] <- lapply(args[-length(args)], function(x) {
                                x$frame <- TRUE
                                x
                              })

                              args.output <- lapply(args, function(x) {
                                x$caption <- NULL
                                x$caption.level <- NULL
                                capture.output(x$show.t2t())
                              })

                              args.output[-1] <- lapply(args.output[-1], function(x) {
                                x <- sub("^\\|+", "", x)
                                x
                              })

                              cat(header.t2t(caption = caption, caption.level = caption.level))
                              cat(do.call("paste", c(args.output, list(sep = ""))), sep = "\n")
                            },

                            show.textile = function(args = .self$args, caption = .self$caption, caption.level = .self$caption.level, frame = .self$frame, grid = .self$grid, col.width = .self$col.width, width = .self$width) {
                              test.nrow(args)

                              args.output <- lapply(args, function(x) {
                                x$caption <- NULL
                                x$caption.level <- NULL
                                x$frame <- NULL
                                x$width <- 0
                                capture.output(x$show.textile())
                              })

                              args.output[-length(args.output)] <- lapply(args.output[-length(args.output)], function(x) {
                                x <- sub("\\|$", "", x)
                                x
                              })

                              cat(header.textile(caption = caption, caption.level = caption.level, frame = frame, width = width))
                              cat(do.call("paste", c(args.output, list(sep = ""))), sep = "\n")
                            },

                            show.pandoc = function(args = .self$args, caption = .self$caption, caption.level = .self$caption.level, frame = .self$frame, grid = .self$grid, col.width = .self$col.width, width = .self$width) {
                              test.nrow(args)

                              args.output <- lapply(args, function(x) {
                                x$caption <- NULL
                                x$caption.level <- NULL
                                capture.output(x$show.pandoc())
                              })

                              args.output[-length(args.output)] <- lapply(args.output[-length(args.output)], function(x) {
                                x <- paste(x, " ", sep = "")
                                x
                              })

                              cat(header.pandoc(caption = caption, caption.level = caption.level))
                              cat(do.call("paste", c(args.output, list(sep = ""))), sep = "\n")
                            }
                            )
                          )

##' Cbind two ascii objects
##'
##' This function binds cols of two ascii table.
##' @title Cbind two ascii objects
##' @param ... ascii objects
##' @param caption see \code{?ascii}
##' @param caption.level see \code{?ascii}
##' @param frame see \code{?ascii}
##' @param grid see \code{?ascii}
##' @param col.width see \code{?ascii}
##' @param width see \code{?ascii}
##' @return An \code{"asciiCbind"} object.
##' @export
##' @author David Hajage
cbind.ascii <- function(..., caption = NULL, caption.level = NULL, frame = NULL, grid = NULL, col.width = 1, width = 0) {
  results <- asciiCbind$new(args = list(...), caption = caption, caption.level = caption.level, frame = frame, grid = grid, col.width = col.width, width = width)
  results
}

##' @rdname print-ascii
##' @export
setMethod("print","asciiCbind",
           function(x, type = getOption("asciiType"), file = NULL, append = FALSE, escape = FALSE, list.escape = c("\\_", "\\^"), ...) {
             if (type == "asciidoc") res <- capture.output(x$show.asciidoc())
             if (type == "rest") res <- capture.output(x$show.rest())
             if (type == "org") res <- capture.output(x$show.org())
             if (type == "t2t") res <- capture.output(x$show.t2t())
             if (type == "textile") res <- capture.output(x$show.textile())
             if (type == "pandoc") res <- capture.output(x$show.pandoc())

             if (escape) {
               for (i in list.escape)
                 res <- gsub(i, paste("\\", i, sep = ""), res)
             }

             if (is.null(file)) {
               cat(res, sep = "\n")
             }
             else {
               if (append) op <- "a" else op <- "w"
               f <- file(file, op)
               writeLines(res, f)
               close(f)
             }
             invisible(x)
           }
           )

##' @rdname print-ascii
##' @export
setMethod("show","asciiCbind",
           function(object) {
             print(object)
           }
          )
