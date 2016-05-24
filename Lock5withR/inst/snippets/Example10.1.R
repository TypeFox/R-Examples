lm(Price ~ PPM + CostBW, data = InkjetPrinters)
Ink.Price <- makeFun(lm(Price ~ PPM + CostBW, data = InkjetPrinters))
Ink.Price(PPM = 3.0, CostBW = 3.7)

