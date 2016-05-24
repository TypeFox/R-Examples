## =============================================================================
## Diagram of the "measel" model, BVP chapter
## Figure 12.6 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diagram)
pm <- par(mar = c(0, 0, 0, 0))
cex1 <- 1.5
openplotmat()
names <- c(expression(Susceptible~(y[1])), "", expression(Latent~(y[2])),
           expression(Infective~(y[3])))

pos <- matrix(ncol = 2, byrow = TRUE,
              data = c(0.25, 0.75, 0.45, 0.5,
                       0.85, 0.5,  0.25, 0.25))
M <- matrix(nrow = 4, byrow = TRUE, data = c(0, 0, 0, 0,
                                             1, 0, 0, 1,
                                             0, 1, 0, 0,
                                             0, 0, 1, 0))
arr.length <- M
M <- as.data.frame(M)
M[2,1] <- M[2,4] <-""
M[3,2] <-"beta * y[1]* y[3]"
M[4,3] <-"y[2]/lambda"

pp <- straightarrow (pos[1,] + c(0,0.2), pos[1,],
                     arr.pos = 0.4, arr.type = "triangle")
text(pp[1]+0.035, pp[2], expression(mu), cex = cex1)

mv <- pos[4,]
pp <- bentarrow(mv, mv - c(0.12, 0.15),
                arr.pos = 1, arr.type = "triangle")
text(pp[1]+0.05, pp[2], expression(y3/eta), cex = cex1)

arr.length[] <- 0.4
arr.length[[2,1]] <- arr.length[[2, 4]] <- 0

pp<-plotmat(M, pos = pos, curve = 0, name = names, lwd = 1, cex = cex1,
            box.lwd = 2, box.size = c(0.135, 0, 0.135,0.135),
            box.type = c("square", "none", "square","square"),
            arr.type = "triangle", dr = 0.001,
            box.prop = c(0.35, 0.35, 0.35, 0.35),
            arr.lwd = 2, arr.length = arr.length,
            arr.pos = 0.4, add = TRUE, box.cex = 1.45)
par(mar = pm)


