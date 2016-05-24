p <- 0.275
SE <- 0.03
MoE <- 2 * SE
p - MoE
p + MoE
dotPlot(~ prop, width = .005, groups = (0.215 <= prop & prop <= 0.335), data = Sampledist.deg)

