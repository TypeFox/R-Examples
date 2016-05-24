ColorPalette<-function(colors=c("red","green"),ncols=5){
	my_palette=grDevices::colorRampPalette(colors)(ncols)
	
	return(my_palette)
	
}
