
xcolors<-function(max_rank=-1) {
	if (max_rank>0 && max_rank<NROW(.color_data))
		.color_data$color_name[seq_len(max_rank)]
	else
		.color_data$color_name	
	}

name2color<-function(name,exact=TRUE,hex_only=TRUE,n=-1){
	if(exact){
		d<-.color_data[match(name,.color_data$color_name),]
	} else {
	   d<-.color_data[grep(name,.color_data$color_name),]
	   if (n>0 && nrow(d)>n) d<-d[seq_len(n),]
	}
	if (hex_only) {
		d$hex
	} else{
		d
	}
}

nearest_named<-function(color, hex_only=FALSE,max_rank=-1,Lab=TRUE){
	if (NCOL(color)==3){
		rgbcol<-color
	} else if (is.character(color)){
		rgbcol<-t(col2rgb(color))
	} else if (is.factor(color)){#sigh
		rgbcol<-t(col2rgb(as.character(color)))
	}
	if (max_rank>0 & max_rank<nrow(.color_data)) 
		ranks<-seq_len(max_rank) 
	else	
		ranks<-1:nrow(.color_data)
		
	if (Lab){
		labcol<-convertColor(rgbcol/255,from="sRGB",to="Lab")
		nearest<-knnx.index(as.matrix(.color_data[ranks,c("L","a","b")]),query=labcol,k=1)
	} else {
		nearest<-knnx.index(as.matrix(.color_data[ranks,c("red","green","blue")]),query=rgbcol,k=1)
	}
	if (hex_only)
		.color_data[nearest,"hex"]
	else
		.color_data[nearest,]
}