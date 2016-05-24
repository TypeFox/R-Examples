Sampledist.deg <- do(1000) * rflip(200, 0.275) # 1000 samples, each of size 200 and proportion 0.275
head(Sampledist.deg, 3)
dotPlot(~ prop, width = .005, data = Sampledist.deg)

