library(oaPlots)
library(RColorBrewer)
library(ggplot2)

display.brewer.all()


# modify the data
dsub <- subset(diamonds, x > 5 & x < 6 & y > 5 & y < 6)
dsub <- dsub[-which(dsub$z > 4), ]
dsub <- dsub[-which(dsub$z < 3), ]


colorPalette <- brewer.pal(9, "Blues")[4:9]
colorObj <- splitColorVar(colorVar = dsub$z, colorPalette)
colorVec <- colorObj$colorVec
breaks <- colorObj$breaks

prepLegend(side = "right", proportion = 0.3)
oaTemplate(xlim = range(dsub$x), ylim = range(dsub$y), 
		main = "Diamond Length by Width \n Colored by Depth",
		xlab = "Length (mm)", ylab = "Width (mm)")
points(x = dsub$x, y = dsub$y, col = colorVec, pch = 19, cex = 0.6)

densityLegend(x = dsub$z, colorPalette = colorPalette, side = "right",
		main = "Diamond Depth", colorBreaks = breaks)



head(mtcars)
colorPalette <- brewer.pal(9, "YlOrRd")[4:9]
scatterplotDL(x = mtcars$mpg, y = mtcars$wt, 
		colorVar = mtcars$hp, legendTitle = "Horse Power", colorPalette = colorPalette, pch = 19,
		xlab = "MPG (miles per gallon)", ylab = "Weight (tonnes)",  
		main = "MPG by Weight in Cars \n Colored by Horse Power")

scatterplotDL(x = mtcars$mpg, y = mtcars$wt, colorVar = mtcars$disp, 
		colorPalette = colorPalette, side = "bottom",
		pch = 19, 
		legendTitle = "Displacement (cubic inches)", 
		xlab = "MPG (miles per gallon)", ylab = "Weight (tonnes)", 
		main = "MPG by Weight in Cars \n Colored by Displacement")
