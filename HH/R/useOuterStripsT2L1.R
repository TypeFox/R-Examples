useOuterStripsT2L1 <- function(x, ..., strip.height=.4, strip.names=c(TRUE, TRUE)) {
  if (length(dim(x)) != 3)
    stop("useOuterStripsT2L1 requires a three-dimensional trellis object.",
           call. = FALSE)
  layout <- x$layout
  if (is.null(layout)) {
    layout <- c(prod(dim(x)[1:2]), dim(x)[3])
    x$layout <- layout
  }
  update(x,
         par.settings=list(
           layout.heights=list(strip=c(rep(0, layout[2]-1), 2*strip.height)),
           layout.widths=list(strip.left=c(1*strip.height, rep(0, layout[1]-1)))),
         strip=function(which.given, which.panel, var.name, ...) {
           ## if (which.given==1) cat(which.given, which.panel, var.name, packet.number(), "\n")
           if (which.given == 1)
             strip.default(strip.names=strip.names,
                           which.given = 1,
                           which.panel = which.panel[1:2],
                           var.name = var.name[1:2],
                           ...)
           ## if (which.given==2) cat(which.given, which.panel, var.name, packet.number(), "\n")
           if (which.given == 2){
             strip.default(strip.names=strip.names,
                           which.given = 2,
                           which.panel = which.panel[1:2],
                           var.name = var.name[1:2],
                           ...)
           }
         },
         strip.left=function(which.given, which.panel, var.name, ...) {
           ## if (which.given==3) cat(which.given, which.panel, var.name, packet.number(), "\n\n")
           if (which.given == 3)
             strip.default(strip.names=strip.names,
                           which.given = 1,
                           which.panel = which.panel[3],
                           var.name = var.name[3],
                           ...)
         }
         )
}
