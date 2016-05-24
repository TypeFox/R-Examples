ink.model <- lm(Price ~ PPM, data = InkjetPrinters)
dotPlot( ~ resid(ink.model), cex = .05, nint = 40)

