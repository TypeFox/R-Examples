"col.hh" <-
function () 
  list(background = list(col = "transparent"),
       bar.fill = list(col = "#c8ffc8"), 
       box.rectangle = list(col = "royalblue"),
       box.umbrella = list(col = "royalblue", lty=1),
       box.dot = list(col="royalblue"),
       dot.line = list(col = "#e8e8e8"),
       dot.symbol = list(col = "royalblue"), 
       plot.line = list(col = "royalblue", pch=16),
       plot.symbol = list(col = "royalblue"), 
       regions = list(col = heat.colors(100)),
       strip.shingle =
       list(col = c("#ff7f00", "#00ff00", "#00ffff", "#0080ff",
              "#ff00ff", "#ff0000", "#ffff00")),
       strip.background =
       list(col = c("#ffe5cc", "#ccffcc", "#ccffff", "#cce6ff",
              "#ffccff", "#ffcccc", "#ffffcc")),
       reference.line = list(col = "#e8e8e8"), 
       superpose.line =
       list(col = c("royalblue", "red", "darkgreen", "brown",
              "orange", "turquoise", "orchid"),
            lty = 1:7), 
       superpose.symbol =
       list(pch = c(1, 3, 6, 0, 5, 16, 17),
            cex = rep(0.7, 7),
            col = c("royalblue", "red", "darkgreen", "brown",
              "orange", "turquoise", "orchid")))

