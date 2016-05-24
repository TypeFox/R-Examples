cdata(~ mean, 0.90, data = Bootstrap)
dotPlot(~ mean, width = .1, groups = (27.70 <= mean & mean <= 30.71), data = Bootstrap)

cdata( ~ mean, 0.99, data = Bootstrap)
dotPlot(~ mean, width = .1, groups = (26.98 <= mean & mean <= 31.63), data = Bootstrap)

