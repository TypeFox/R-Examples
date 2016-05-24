data(Catheter)

## make the data for internal use
outdata = makeData()


dev.new()
par(mfrow=c(2,2))
Bym.map(result1$summary.random$region.struct$mean)
Bym.map(result2$summary.random$region.struct$mean)
Bym.map(result3$summary.random$region.struct$mean)




