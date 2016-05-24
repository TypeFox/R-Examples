t <- read.csv("tableparamVSfunction.csv", sep = ",", header = TRUE, check.names = FALSE)
row.names(t) <- t[, 1]
t <- t[, -1]
t[is.na(t)] <- 0

table.value(t, plegend.drawKey = FALSE, ppoints.cex = 0.2, symbol = "circle", axis.text = list(cex = 0.8), pgrid.draw = TRUE, 
            ptable.y = list(srt = 45, pos = "left"), 
            ptable.margin = list(bottom = 2, left = 15, top = 15, right = 2))
            
            
