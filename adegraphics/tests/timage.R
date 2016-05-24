library(adegraphics)
pdf("timage.pdf")

## ex1
x <- 1:4
y <- 1:4
df <- data.frame(as.matrix(cbind(x, y)))
g1 <- table.image(df, col = 2:4)
update(g1, plegend.drawColorKey = TRUE)

## ex2
df <- matrix(0, 10, 10)
df[1:3, 1:3] <- 5        
g2 <- table.image(df)
g3 <- table.image(df, breaks = c(5, 2, 0))

## ex3
data(rpjdl, package = "ade4")
X <- data.frame(t(rpjdl$fau))
Y <- data.frame(t(rpjdl$mil))
coa1 <- ade4::dudi.coa(X, scan = FALSE)
x <- rank(coa1$co[, 1])
y <- rank(coa1$li[, 1])
g4 <- table.image(Y, coordsx = x, coordsy = 1:8, axis.text = list(alpha = 0), pgrid.col = "black", pgrid.lwd = 0.8, col = c("white", "black"), plegend.drawKey = FALSE)
g5 <- table.image(X, coordsx = x, coordsy = y, ptable = list(x = list(tck = 0), y = list(tck = 4)), pleg.drawKey = FALSE, labelsy = paste(" ", row.names(X), sep = ""))
g6 <- ADEgS(list(g4, g5), positions = rbind(c(0, 0, 1, 0.3), c(0, 0.4, 1, 1)))
