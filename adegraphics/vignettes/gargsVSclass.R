library(grid)
t <- read.csv("gargsVSclass.csv", sep = ",", header = TRUE, check.names = FALSE)
row.names(t) <- t[, 1]
t <- t[, -1]
t[is.na(t)] <- 0

table.value(t, plegend.drawKey = FALSE, ppoints.cex = 0.2, symbol = "circle", axis.text = list(cex = 0.7), pgrid.draw = TRUE,
            ptable.margin = list(bottom = 5, left = 15, top = 15, right = 5),
            ptable.x = list(tck = 5), ptable.y = list(tck = 5, srt = 45, pos = "left"))
            

