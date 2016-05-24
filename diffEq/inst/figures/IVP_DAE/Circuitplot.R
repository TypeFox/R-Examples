## =============================================================================
## Schematic representation of the electrical circuit problem 
## Figure 6.4 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diagram)    # Need at least version 1.6

layoutmat <- matrix(data = c(rep(1, 12), 2, 3, 4, 5),
                    nrow = 4, ncol = 4, byrow = TRUE)
layoutmat
nf <- layout(layoutmat, respect = FALSE)
#layout.show(nf)

par(lwd = 1.5)
par(mar = c(0, 0, 2, 0))
emptyplot(main = "transistor Amplifier",
          ylim = c(-0.1, 1),
          xlim = c(-0.1, 1.1), asp = FALSE)

x1 <- 0.08; x2 <- 0.2; x3 <- 0.3; x4 <- 0.4; x5 <- 0.5; x6 <- 0.6;
x7 <- 0.7;  x8 <- 0.8; x9 <- 0.92;
y1 <- 0.05; y2 <- 0.4; y3 <- 0.5; y4 <- 0.6; y5 <- 0.95 
x23 <- (x2 + x3)/2
x56 <- (x5 + x6)/2
lines(c(x2, x9, x9, x2, x2, x1, x1, x23, x3, x3),
      c(y1, y1, y5, y5, y1, y1, y3, y3,  y4, y5))
lines(c(x23, x3, x3),
      c(y3,  y2, y1))
lines(c(x3,  x4,  x4),
      c(y2,  y2,  y1))
lines(c(x5,  x5),
      c(y1,  y5))
lines(c(x3,  x5),
      c(y4,  y4))
lines(c(x5,  x56),
      c(y3,  y3))
lines(c(x56, x6, x6),
      c(y3,  y4, y5))
lines(c(x56, x6, x6),
      c(y3,  y2, y1))
lines(c(x6,  x7, x7),
      c(y2,  y2, y1))
lines(c(x6,  x8, x8),
      c(y4,  y4, y1))
en.Amplifier(c(x23, y3), r = 0.035)
en.Amplifier(c(x56, y3), r = 0.035)

en.Signal(c(x1, 0.15), lab = expression("U"["in"]))
en.Signal(c(x9, y2), lab = expression("U"["b"]))
straightarrow(c(x1 - 0.065, 0.18), c(x1 - 0.065, 0.12),
                arr.pos = 1, arr.type = "triangle", lwd = 1)
straightarrow(c(x9 + 0.065, y2 + 0.03), c(x9 + 0.065, y2 - 0.03),
                arr.pos =1, arr.type = "triangle", lwd = 1)

en.Node(c(x1, y3), lab = "u1")
en.Node(c(x2, y3), lab = "u2")
en.Node(c(x3, y2), lab = "u3", pos = 2)
en.Node(c(x3, y4), lab = "u4", pos = 2)
en.Node(c(x5, y4), lab = "u5")
en.Node(c(x6, y2), lab = "u6")
en.Node(c(x6, y4), lab = "u7")
en.Node(c(x8, y4), lab = "u8")

en.Capacitator(c(0.5*(x1+x2), y3), lab = "C1", vert = FALSE)
en.Capacitator(c(x4, y4), lab = "C3", vert = FALSE)
en.Capacitator(c(x4, 0.5*(y1+y2)), lab = "C2", vert = TRUE)
en.Capacitator(c(x7, y4), lab = "C5", vert = FALSE)
en.Capacitator(c(x7, 0.5*(y1+y2)), lab = "C4", vert = TRUE)
en.Capacitator(c(0.5*(x1+x2), y3), vert = FALSE)

en.Resistor(c(x1, y2), lab = "R0")
en.Resistor(c(x2, 0.5*(y1+y2)), lab = "R1")
en.Resistor(c(x2, 0.5*(y4+y5)), lab = "R2")
en.Resistor(c(x3, 0.5*(y4+y5)), lab = "R4")
en.Resistor(c(x5, 0.5*(y4+y5)), lab = "R6")
en.Resistor(c(x6, 0.5*(y4+y5)), lab = "R8")
en.Resistor(c(x3, 0.5*(y1+y2)), lab = "R3")
en.Resistor(c(x5, 0.5*(y1+y2)), lab = "R5")
en.Resistor(c(x6, 0.5*(y1+y2)), lab = "R7")
en.Resistor(c(x8, 0.5*(y1+y2)), lab = "R9")

en.Ground(c(0.92, 0.05))
#lines(c(-1,2),c(0,0),col="grey")

par(mar=c(2, 2, 2, 2))

emptyplot(main = "transistor")
lines(c(0.1, 0.5,0.9), c(0.5, 0.5, 0.9))
lines(c(0.5, 0.9), c(0.5, 0.1))
lines(c(0.5, 0.5), c(0.4, 0.6))
text(0.2, 0.4, "Gate", font = 3)
text(0.8, 0.9, "Drain", font = 3,adj = 1)
text(0.8, 0.1, "Source", font = 3,adj = 1)
en.Amplifier(c(0.5, 0.5), r = 0.15)
box(col = "grey")

emptyplot(main = "en.Capacitator")
straightarrow(c(0.5, 0.9), c(0.5, 0.1),
              arr.pos = 0.3, arr.length = 0.25, arr.type = "triangle")
en.Capacitator(c(0.5, 0.5), width = 0.075, length = 0.5, vert = TRUE)
text(0.4, 0.65, "i", font = 3, cex = 2)
straightarrow(c(0.8, 0.3), c(0.8, 0.77), arr.pos = 1,
              arr.length = 0.25, arr.type = "triangle", lwd = 1)
text(0.925, 0.65, "v", font = 3, cex = 2)
text(0.15, 0.5, "C", font = 3, cex = 2)
box(col = "grey")

emptyplot(main = "resistor")
straightarrow(c(0.5, 0.9), c(0.5, 0.1), arr.pos = 0.2,
              arr.length = 0.25, arr.type = "triangle", lwd = 1)
text(0.4, 0.85, "i", font = 3, cex = 2)
en.Resistor(c(0.5, 0.5), width = 0.25, length = 0.35 )
straightarrow(c(0.8, 0.3), c(0.8, 0.77), arr.pos = 1,
              arr.length = 0.25, arr.type = "triangle", lwd = 1)
text(0.925, 0.65, "v", font = 3, cex = 2)
text(0.5, 0.5, "R", font = 3, cex = 2)
box(col = "grey")

emptyplot(main = "voltage source")
lines(c(0.5, 0.5), c(0.1, 0.9))
en.Signal(c(0.5, 0.5), r = 0.15)
straightarrow(c(0.8, 0.3), c(0.8, 0.77), arr.pos = 1,
              arr.length = 0.25, arr.type = "triangle", lwd = 1)
text(0.925, 0.65, "v", font = 3, cex = 2)
box(col = "grey")
