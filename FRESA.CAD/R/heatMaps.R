heatMaps <- function(variableList,varRank=NULL,Outcome,data,title="Heat Map",hCluster=FALSE,prediction=NULL,outcomeGain=0,theFiveColors=c("blue","cyan","black","yellow","red"),...) 
  {
    
    par(pty='m')
    if (!requireNamespace("gplots", quietly = TRUE)) {
      install.packages("gplots", dependencies = TRUE)
    } 
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      install.packages("RColorBrewer", dependencies = TRUE)
    } 
    
    # creates a own color palette from red to blue
    my_palette <- colorRampPalette(theFiveColors)(n = 224)
    #	my_palette <- colorRampPalette(c("blue","light blue","black","pink","red"))(n = 224)
    
    # defines the color breaks manually for a "skewed" color transition
    col_breaks = c(
      seq(-2.0,-1.101,length=50),  		# for blue
      seq(-1.10,-0.201,length=50),          # for cyan
      seq(-0.20,0.20,length=25),          # for black
      seq(0.201,1.10,length=50),            # for yellows
      seq(1.101,2.0,length=50))             # for red
    
    
    if (class(variableList)=="data.frame")
    {
      vnames <- as.vector(variableList[,1]);
    }
    else
    {
      vnames <- names(variableList);
    }
    frm <- paste(Outcome," ~ ",1);
    added = 0;
    for (i in 1:length(vnames))
    {
      if ((vnames[i] != "")&&(vnames[i] != "1"))
      {
        frm <- paste(frm," + ",vnames[i])
        added = added + 1;
      }
    }
    modelFrame <- model.frame(formula(frm),data);
    rownames(modelFrame) <- rownames(data)
    
    if (is.null(varRank)) 
    {
      topvarID <- seq(1, added, 1)
      hits <- seq(1, added, 1)
    }
    else
    {
      topvarID <- as.numeric(rownames(varRank));
      hits <- as.vector(as.matrix(varRank))
    }
    if (outcomeGain==0)
    {
      orderData <- cbind(data[,Outcome]);
    }
    else
    {
      orderData <- cbind(outcomeGain*data[,Outcome]-outcomeGain/2);
    }
    orderData <- as.data.frame(orderData);
    rownames(orderData) <- rownames(data);
    colnames(orderData) <- c(Outcome);
    cn = 1; 
    if (!is.null(prediction))
    {
      cn = cn + 1;
      orderData <- cbind(orderData,as.vector(prediction));
      colnames(orderData)[cn] <- "prediction";
    }
    for ( i in 1:length(topvarID))
    {
      if (hits[i]>0)
      {
        orderData <- cbind(orderData,modelFrame[,topvarID[i]+1])
        cn = cn + 1;
        colnames(orderData)[cn] <- colnames(modelFrame)[topvarID[i]+1];
      }
    }
    
    dataMat <- as.matrix(orderData);
    
    
    if (!is.null(prediction))
    {
      orderData <- eval(parse(text=paste("with(orderData,orderData[order(-",Outcome,",-orderData[,2]),])")))
    }
    else
    {
      orderData <- eval(parse(text=paste("with(orderData,orderData[order(-",Outcome,",-orderData[,2]),])")))
      cm <- as.data.frame(cor(t(orderData),method="spearman"))
      rownames (cm) <- c(1:nrow(cm))
      colnames (cm) <- c(1:nrow(cm))
      indx <- 1
      ix <- 1
      rm <- cm;
      for (i in 1:(nrow(cm)-2))
      {
        rm <- cm[ix,-indx]
        ix <- as.integer(names(which.max(rm)))
        indx <- append(indx,ix)
        #		  cat (ix,",")
      }
      orderData <- orderData[indx,];
      orderData <- eval(parse(text=paste("with(orderData,orderData[order(-",Outcome,"),])")))
    }
    
    
    orderData <- as.matrix(orderData);
    
    if (outcomeGain==0) orderData <- scale(orderData);

    index <- ((orderData[,1]+2)/4)
    index <- index*(index>0);
    index <- index*(index<1)+0.999*(index>=1);
    index <- floor(nrow(orderData)*index)+1;  		
    rowcolors <- colorRampPalette(theFiveColors)(n = nrow(orderData))[index];
    orderData <- orderData[,-1];
    
    if (is.null(prediction))
    {
      if (hCluster==TRUE)
      {
        heatMap <- gplots::heatmap.2(orderData, 
                                     main = title, 			# heat map title
                                     notecol="black",      	# change font color of cell labels to black
                                     density.info="none",  	# turns off density plot inside color legend
                                     trace="none",         	# turns off trace lines inside the heat map
                                     margins =c(12,9),     	# widens margins around plot
                                     col=my_palette,       	# use on color palette defined earlier 
                                     breaks=col_breaks,    	# enable color transition at specified limits
                                     dendrogram="both", 		# raw and column dendrogram
                                     RowSideColors=rowcolors,
                                     lmat=rbind(c(0,5,4,0,0), c(0,3,2,1,0)),
                                     lhei=c(1.5,4.5),
                                     lwid=c(0.25,2,6,0.25,0.5),
                                     ...) 				# turn off column clustering
      }
      else
      {
        #			rownames(orderData) <- format(orderData[,1],digits=4);

        heatMap <- gplots::heatmap.2(orderData, 
                                     main = title, 			# heat map title
                                     notecol="black",      	# change font color of cell labels to black
                                     density.info="none",  	# turns off density plot inside color legend
                                     trace="none",         	# turns off trace lines inside the heat map
                                     margins =c(12,9),     	# widens margins around plot
                                     col=my_palette,       	# use on color palette defined earlier 
                                     breaks=col_breaks,    	# enable color transition at specified limits
                                     dendrogram="col", 		# cluster by column
                                     RowSideColors=rowcolors,
                                     lmat=rbind(c(0,5,4,0,0), c(0,3,2,1,0)),
                                     lhei=c(1.5,4.5),
                                     lwid=c(0.25,2,6,0.25,0.5),
                                     Rowv="NA",
                                     ...) 				# turn off row clustering				
      }
    }
    else
    {
      if (hCluster==FALSE)
      {
        heatMap <- gplots::heatmap.2(orderData, 
                                     main = title, 			# heat map title
                                     notecol="black",      	# change font color of cell labels to black
                                     density.info="none",  	# turns off density plot inside color legend
                                     trace="none",         	# turns off trace lines inside the heat map
                                     margins =c(12,9),     	# widens margins around plot
                                     col=my_palette,       	# use on color palette defined earlier 
                                     breaks=col_breaks,    	# enable color transition at specified limits
                                     dendrogram="none", 		# sorted by prediction
                                     Colv="NA",				# turn off column clustering
                                     Rowv="NA",				# turn off row clustering
                                     RowSideColors=rowcolors,
                                     colsep=1,
                                     lmat=rbind(c(0,5,4,0,0), c(0,3,2,1,0)),
                                     lhei=c(1.5,4.5),
                                     lwid=c(0.25,2,6,0.25,0.5),
                                     ...) 				
      }
      else
      {
        #			rownames(orderData) <- format(orderData[,1],digits=4);
#        index <- ((orderData[,1]+2)/4)
#        index <- index*(index>0);
#        index <- index*(index<1)+0.999*(index>=1);
#        index <- floor(nrow(orderData)*index)+1;			
#        rowcolors <- colorRampPalette(theFiveColors)(n = nrow(orderData))[index];
#        orderData <- orderData[,-1];
        heatMap <- gplots::heatmap.2(orderData, 
                                     main = title, 			# heat map title
                                     notecol="black",      	# change font color of cell labels to black
                                     density.info="none",  	# turns off density plot inside color legend
                                     trace="none",         	# turns off trace lines inside the heat map
                                     margins =c(12,9),     	# widens margins around plot
                                     col=my_palette,       	# use on color palette defined earlier 
                                     breaks=col_breaks,    	# enable color transition at specified limits
                                     dendrogram="col", 		# 
                                     RowSideColors=rowcolors,
                                     lmat=rbind(c(0,5,4,0,0), c(0,3,2,1,0)),
                                     lhei=c(1.5,4.5),
                                     lwid=c(0.25,2,6,0.25,0.5),
                                     Rowv="NA",				# turn off row clustering
                                     ...)
      }
    }
    par(pty='m',par(mfrow=c(1,1)))

    result <- list(dataMatrix=dataMat,
                   orderMatrix=orderData,
                   heatMap=heatMap)
    return (result);
  }
