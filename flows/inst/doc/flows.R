## ------------------------------------------------------------------------
library(flows)
# Import data
data(nav)
head(nav)
# Prepare data
myflows <- prepflows(mat = nav, i = "i", j = "j", fij = "fij")
myflows[1:4,1:4]

## ---- fig.width=5, fig.height=5------------------------------------------
# Import data
data(nav)
myflows <- prepflows(mat = nav, i = "i", j = "j", fij = "fij")

# Get statistics about the matrix
statmat(mat = myflows, output = "none", verbose = TRUE)

# Plot Lorenz curve only
statmat(mat = myflows, output = "lorenz", verbose = FALSE)

## ---- fig.width=7, fig.height=7------------------------------------------
# Graphics only
statmat(mat = myflows, output = "all", verbose = FALSE)

# Statistics only
mystats <- statmat(mat = myflows, output = "none", verbose = FALSE)
str(mystats)
# Sum of flows
mystats$sumflows


## ---- fig.height = 4, fig.width=4----------------------------------------
# Import data
data(nav)
myflows <- prepflows(mat = nav, i = "i", j = "j", fij = "fij")

# Remove the matrix diagonal
diag(myflows) <- 0

# Selection of flows > 500
flowSel1 <- firstflowsg(mat = myflows, method = "xfirst", k = 500)
# Selection of flows > 1000
flowSel2 <- firstflowsg(mat = myflows, method = "xfirst", k = 1000)

# Compare initial matrix and selected matrices
compmat(mat1 = myflows, mat2 = myflows * flowSel1, digits = 1)
compmat(mat1 = myflows, mat2 = myflows * flowSel2, digits = 1)

## ---- fig.height = 6, fig.width=6----------------------------------------
# Import data
data(nav)
myflows <- prepflows(mat = nav, i = "i", j = "j", fij = "fij")

# Remove the matrix diagonal
diag(myflows) <- 0

# Percentage of each outgoing flows
myflows2 <- myflows / rowSums(myflows) * 100

# Select flows that represent at least 20% of the sum of outgoing flows for 
# each urban area.
flowSel <- firstflows(mat = myflows2, method = "xfirst", k = 20)

# Compare initial and selected matrices
compmat(mat1 = myflows,mat2 = flowSel * myflows)


## ---- fig.height=7, fig.width=7, eval = TRUE-----------------------------
# Import data
data(nav)
myflows <- prepflows(mat = nav, i = "i", j = "j", fij = "fij")

# Remove the matrix diagonal
diag(myflows) <- 0

# Select flows that represent at least 20% of the sum of outgoing flows for 
# each urban area.
flowSel1 <- firstflows(mat = myflows/rowSums(myflows)*100, method = "xfirst", 
                       k = 20)


# Select the dominant flows (incoming flows criterion)
flowSel2 <- domflows(mat = myflows, w = colSums(myflows), k = 1)

# Combine selections
flowSel <- myflows * flowSel1 * flowSel2

# Node weights
inflows <- data.frame(id = colnames(myflows), w = colSums(myflows))

# Plot dominant flows map
opar <- par(mar = c(0,0,2,0))
sp::plot(GE, col = "#cceae7", border = NA)
plotMapDomFlows(mat = flowSel, spdf = UA, spdfid = "ID", w = inflows, wid = "id",
                wvar = "w", wcex = 0.05, add = TRUE,
                legend.flows.pos = "topright",
                legend.flows.title = "Nb. of commuters")
title("Dominant Flows of Commuters")
mtext(text = "INSEE, 2011", side = 4, line = -1, adj = 0.01, cex = 0.8)
par(opar)


# Statistics on major urban areas
inflows <- data.frame(id = colnames(flowSel), w = colSums(flowSel))
UA.df <- unique(data.frame(id = c(nav$i, nav$j),name = c(nav$namei, nav$namej)))
UAindegreew <- merge(inflows, UA.df, by = "id", all.x = TRUE)
UAindegreew[order(UAindegreew$w, decreasing = TRUE),][1:10,]


