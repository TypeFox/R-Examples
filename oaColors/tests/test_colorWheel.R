library(oaColors)

pdf("paletteTest.pdf", width = 6, height = 6)

for( i in 1:9){
	myColors <- oaPalette(i, alpha = 0.7)
	colorWheel(myColors)
}

dev.off()


