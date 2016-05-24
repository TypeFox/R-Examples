###############################################################################
#
# webvis: An R package to create web visualizations.
# author: Shane Conway <shane.conway@gmail.com>
#
# demo example from http://vis.stanford.edu/protovis/ex/wheat.html
#
# This is released under a BSD license.
#
# Documentation was created using roxygen:
# roxygenize('webvis', roxygen.dir='webvis', copy.package=FALSE, unlink.target=FALSE)
#
###############################################################################


# Create the monarchy data
monarch <- read.csv(con <- textConnection(
"name, start, end, commonwealth
Elizabeth, 1565, 1603, 0
James I, 1603, 1625, 0
Charles I, 1625, 1649, 0
Cromwell, 1649, 1660, 1
Charles II, 1660, 1685, 0
James II, 1685, 1689, 0
W&M, 1689, 1702, 0
Anne, 1702, 1714, 0
George I, 1714, 1727, 0
George II, 1727, 1760, 0
George III, 1760, 1820, 0
George IV, 1820, 1821, 0"), header=TRUE, as.is=TRUE)
close(con)
monarch$name <- sapply(monarch$name, function(x) collapse("'", x, "'"))

# Create the wheat data
wheat <- read.csv(con <- textConnection(
				"year, wheat, wages
						1565, 41, 5
						1570, 45, 5.05
						1575, 42, 5.08
						1580, 49, 5.12
						1585, 41.5, 5.15
						1590, 47, 5.25
						1595, 64, 5.54
						1600, 27, 5.61
						1605, 33, 5.69
						1610, 32, 5.78
						1615, 33, 5.94
						1620, 35, 6.01
						1625, 33, 6.12
						1630, 45, 6.22
						1635, 33, 6.3
						1640, 39, 6.37
						1645, 53, 6.45
						1650, 42, 6.5
						1655, 40.5, 6.6
						1660, 46.5, 6.75
						1665, 32, 6.8
						1670, 37, 6.9
						1675, 43, 7
						1680, 35, 7.3
						1685, 27, 7.6
						1690, 40, 8
						1695, 50, 8.5
						1700, 30, 9
						1705, 32, 10
						1710, 44, 11
						1715, 33, 11.75
						1720, 29, 12.5
						1725, 39, 13
						1730, 26, 13.3
						1735, 32, 13.6
						1740, 27, 14
						1745, 27.5, 14.5
						1750, 31, 15
						1755, 35.5, 15.7
						1760, 31, 16.5
						1765, 43, 17.6
						1770, 47, 18.5
						1775, 44, 19.5
						1780, 46, 21
						1785, 42, 23
						1790, 47.5, 25.5
						1795, 76, 27.5
						1800, 79, 28.5
						1805, 81, 29.5
						1810, 99, 30
						1815, 78, 30
						1820, 54, 30
						1821, 54, 30"), header=TRUE)
close(con)
wheat$wages <- as.numeric(wheat$wages)
#save(monarch, wheat, file="c:/Programming/src/R/webvis/data/pw.demo.Rda")

# Produce the chart for the wheat/wages data
wv <- pv.panel(right=60, top=20, bottom=20, width=800, height=445)

wt <- pv.area(wv=wv, data=wheat, interpolate="step-after", height.name="wheat", left.name="year", left.scale="linear.year.x", ymin=0, fill.style="#aaa", stroke.style="#999", scale.min=0)
wt <- wt + pv.rule()
wv <- wv + wt

wg <- pv.area(wv=wv, data=wheat, height.name="wages", left.name="year", fill.style="hsla(195, 50%, 80%, .75)", scale.min=0, scale.max=100)
wg <- wg + pv.line(line.width=4, stroke.style="lightcoral", anchor="top", bottom.name=NULL, left.name=NULL, fill.style="null")
wg <- wg + pv.line(line.width=1.5, stroke.style="black", anchor="top", bottom.name=NULL, left.name=NULL, fill.style="null")
wv <- wv + wg

wv <- wv + pv.label(left=130, bottom=31, font="italic 10px serif", text="Weekly Wages of a Good Mechanic")

# y-axis
wr <- pv.rule(bottom=-0.5)
wr2 <- pv.rule(wv=wv, data=data.frame(y=seq(0,100,10)), bottom.name="y", bottom.scale="linear.y.y", ymin=0, stroke.style="rgba(255, 255, 255, .2)")
wr2 <- wr2 + pv.label(anchor="right", text.name="y")
wr <- wr + wr2
wv <- wv + wr

# x-axis
wr <- pv.rule(wv=wv, data=data.frame(y=seq(1560, 1830, 10)), bottom=0, left.scale="linear.y.x", left.name="y", height=-4)
wr2 <- pv.rule(wv=wv, data=data.frame(y=seq(1600, 1850, 50)), height="null", top=0, left.name="y", left.scale="linear.y.x", stroke.style="rgba(0, 0, 0, .2)")
wr2 <- wr2 + pv.label(anchor="bottom", text.name="y", text.margin=8)
wr <- wr + wr2
wv <- wv + wr

# Example of using the pv.mark() function to create a custom object.
monarch2 <- monarch
monarch2$top <- ifelse(monarch2$commonwealth == 0 & as.numeric(rownames(monarch2)) %% 2 == 0, 15, 10)
monarch2$fill <- ifelse(monarch2$commonwealth == 0, "'#000'", "'#fff'")
monarch2$reign <- (monarch2$end-monarch2$start) * (wv$width/(max(monarch2$end)-min(monarch2$start)))
vm <- new.webvis(root=pv.mark(wv=wv, data=monarch2, type="Bar",  
		pv.param(name="data", value="d"), 
		pv.param(name="height", value=5),
		pv.param(name="strokeStyle", value="#000"),
		pv.param(name="left", data.name="start", scale="linear.start.x"),
		pv.param(name="top", data.name="top"),
		pv.param(name="width", data.name="reign"),
		pv.param(name="fillStyle", data.name="fill")
))
vm <- vm + pv.mark(wv=wv, data=monarch2, type="Label", 
		pv.param(name="font", value="italic 10px serif"), 
		pv.param(name="text", data.name="name"), 
		pv.param(name="textMargin", value=6), 
		pv.param(name="textBaseline", value="top"), 
		anchor="center")
wv2 <- wv + vm

# Render the output with a browser
render.webvis(wv=wv2, vis.name="playfairs_wheat")
