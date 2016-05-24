pdf(file = "Examples/graphicMargins.pdf", width = 4,  height = 4)

plot.new()
box()

mtext("Bottom Margin: side=1", side=1, line = 1, font = 3)
mtext("Left Margin: side=2", side=2, line = 1, font = 3)
mtext("Top Margin: side=3", side=3, line = 1, font = 3)
mtext("Right Margin: side=4", side=4, line = 1, font = 3)

text("Plot  or Panel \n Region", x=.5, y=.5, adj=.5, cex=2, srt=45.)

dev.off()


