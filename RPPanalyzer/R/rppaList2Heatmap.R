`rppaList2Heatmap` <-
function (x,sampledescription="sample",side.color="tissue"
			, remove=c("blank","protein","Abmix"),distance = "eucsq"
         , dendros="both", cutoff=0.005, fileName=NULL
         , cols=colorpanel(100, low="blue",mid="yellow",high="red"), 
         hclust.method="ward", scale="row"){

       
   data <- select.measurements(x)

	mat <- data[[1]]

	rownames(mat) <- data[[4]][,sampledescription]
	colnames(mat) <- data[[3]]["target",]

	sel.cols <- which(is.na(match(colnames(mat),remove)))

	mat <- mat[,sel.cols]
	
	groups <-unique(data[[4]][,side.color])

   rsc <- match(data[[4]][,side.color],groups)

	#colors <- rainbow(length(groups))
    colors <- colorRampPalette(c("red", "orange", "blue"),space = "Lab")(length(groups))

	for (i in seq(along=groups)){
		rsc[rsc==i]=colors[i]
	}
	
	if(!is.null(fileName)) {	
		pdf(file=fileName)
	}
		
      plotHeatmap(t(mat), distance = distance, dendros=dendros, cutoff=cutoff
                  , toFile=FALSE, fileName=fileName
                  , cols=cols
                  , ColSideColors=rsc
                  , hclust.method=hclust.method
                  , scale=scale)

     legend(x=0,y=0.8,legend=groups,col=colors,pch=15)
     
     if(!is.null(fileName)) {
		dev.off()
	}
                  
   }

