Boot.Ink <- do(1000) * lm(Price ~ PPM, data = resample(InkjetPrinters))
favstats( ~ PPM, data = Boot.Ink)
dotPlot( ~ PPM, width = 2, data = Boot.Ink)
Rand.Ink <- do(1000) * lm(Price ~ shuffle(PPM), data = InkjetPrinters)
favstats( ~ PPM, data = Rand.Ink)
dotPlot( ~ PPM, width = 2, data = Rand.Ink)

