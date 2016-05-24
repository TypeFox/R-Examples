msummary(lm(Price ~ PPM, data = InkjetPrinters)) 
confint(lm(Price ~ PPM, data = InkjetPrinters) , "PPM")

