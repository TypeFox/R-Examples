#############################################################################

.mclust <- structure(list(
  emModelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI",
                   "EEE", "EVE", "VEE", "VVE", 
                   "EEV", "VEV", "EVV", "VVV"),
  hcModelNames = c("VVV", "EEE", "VII", "EII"),
  hcUse = "VARS",
  bicPlotSymbols = structure(c(17, 2, 16, 10, 13, 1,
                               15, 5, 8, 9,
                               12, 7, 14, 0,
                               17, 2),
                             .Names = c("EII", "VII", "EEI", "EVI", "VEI", "VVI",
                                        "EEE", "EVE", "VEE", "VVE",
                                        "EEV", "VEV", "EVV", "VVV",
                                        "E", "V")), 
  bicPlotColors = structure(
    { pal <-  grDevices:::colorRampPalette(c("forestgreen", 
                                             "royalblue1", 
                                             "red3"), space = "Lab")
      c("gray", "black", pal(12), "gray", "black")
    },
                            .Names = c("EII", "VII", "EEI", "EVI", "VEI", "VVI",
                                       "EEE", "EVE", "VEE", "VVE",
                                       "EEV", "VEV", "EVV", "VVV",
                                       "E", "V")),
  classPlotSymbols = c(16, 0, 17, 3, 15, 4, 1, 8, 2, 7, 
                       5, 9, 6, 10, 11, 18, 12, 13, 14),
  classPlotColors = c("dodgerblue2", "red3", "green3", "slateblue", 
                      "darkorange", "skyblue1", "violetred4", "forestgreen",
                      "steelblue4", "slategrey", "brown", "black",
                      "darkseagreen", "darkgoldenrod3", "olivedrab",
                      "royalblue", "tomato4", "cyan2", "springgreen2"),
  warn = FALSE))

mclust.options <- function(...)
{
  current <- .mclust
  if(nargs() == 0) return(current)
  args <- list(...)
  if(length(args) == 1 && is.null(names(args))) 
  { arg <- args[[1]]
    switch(mode(arg),
           list = args <- arg,
           character = return(.mclust[[arg]]),
           stop("invalid argument: ", dQuote(arg)))
  }
  if(length(args) == 0) return(current)
  n <- names(args)
  if (is.null(n)) stop("options must be given by name")
  changed <- current[n]
  current[n] <- args
  if(sys.parent() == 0) env <- asNamespace("mclust") else env <- parent.frame()
  assign(".mclust", current, envir = env)
  invisible(current)
}
