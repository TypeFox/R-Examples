
# Generate a list of viewports that correspond to the current base
# inner region, figure region, and plot region
baseViewports <- function() {
  omi <- par("omi")
  innervp <- viewport(x=unit(omi[2], "inches"),
                      y=unit(omi[1], "inches"),
                      width=unit(1, "npc") -
                        unit(omi[1], "inches") -
                        unit(omi[3], "inches"),
                      height=unit(1, "npc") -
                        unit(omi[2], "inches") -
                        unit(omi[4], "inches"),
                      just=c("left", "bottom"))
  fig <- par("fig")
  figurevp <- viewport(x=unit(fig[1], "npc"),
                       y=unit(fig[3], "npc"),
                       width=unit(fig[2] - fig[1], "npc"),
                       height=unit(fig[4] - fig[3], "npc"),
                       just=c("left", "bottom"))
  plt <- par("plt")
  usr <- par("usr")
  logscale <- FALSE
  if (par("xlog") || par("ylog")) {
    warning("viewport scales NOT set to user coordinates")
    logscales <- TRUE
  }
  plotvp <- viewport(x=unit(plt[1], "npc"),
                     y=unit(plt[3], "npc"),
                     width=unit(plt[2] - plt[1], "npc"),
                     height=unit(plt[4] - plt[3], "npc"),
                     just=c("left", "bottom"),
                     xscale=c(usr[1], usr[2]),
                     yscale=c(usr[3], usr[4]))
  list(inner=innervp, figure=figurevp, plot=plotvp)
}
