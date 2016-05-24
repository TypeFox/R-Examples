useOuterStrips2 <-
function (x, strip.lines = 1, strip.left.lines = 1)
{
  par.settings <- modifyList(x$par.settings,list(
                     layout.heights = list(strip = c(rep(0, dim(x)[3] - 1), strip.lines)),
                     layout.widths = list(strip.left = c(strip.left.lines-0.5,rep(0, dim(x)[2] * 2 - 1)))
                     )
                    )
  new.strip <- function(which.given, which.panel, var.name, ...) {
                if (which.given == 1) strip.default(which.given = 1,
                                            which.panel = which.panel[c(1,2)],
                                            var.name = var.name[2],
                                            ...)
                if (which.given == 2) strip.default(which.given = 2,
                                            which.panel = which.panel[c(1,2)],
                                            var.name = var.name[2],
                                            ...)
               }
  new.strip.left <- function(which.given, which.panel, var.name, ...) {
                      if (which.given == 3) strip.default(which.given = 1,
                                                  which.panel = which.panel[3],
                                                  var.name = var.name[3],
                                                  ...)
                    }
  update(x, par.settings = par.settings, strip = new.strip,
         strip.left = new.strip.left, par.strip.text = list(cex = 1.1, lines = 0.75))
}
