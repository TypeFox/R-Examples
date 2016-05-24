##' @export
`IdfMltSV.default` <- function(MultiData,LgLim,SVSize){

  MultiDt <-   (subset(MultiData$MultiSvLg2,abs(MultiData$MultiSvLg2$PopLg2)>=LgLim))
  MultiDt <-   (subset(MultiDt,abs(MultiDt$PopLg2)!="Inf"))
  MultiDt$MidPos <- round((MultiDt$end+MultiDt$start)/2)  
	MultiDt$MultiSV <- 0
	MultiDt$MultiSVSize <- NA
	MultiDt$MultiSVlog2 <- NA
	MultiDt$MultiSVPval <- NA
	ComBin = 1
	PrevRegion = 0
	MltSv <- 1
	StepSize = mean((MultiDt$end-MultiDt$start))+1
	Brdg = StepSize*3
	Lengths <- abs(MultiDt[2:nrow(MultiDt),'MidPos']-MultiDt[1:(nrow(MultiDt)-1),'MidPos'])
	MultiEnd <- c(which(Lengths>Brdg),nrow(MultiDt))
	
	for(MSV in MultiEnd){
  		
  		if(MSV-MltSv+1>=ComBin){
    		PrevRegion <- PrevRegion+1
    		MultiDt[MltSv:MSV,'MultiSV'] <- PrevRegion
  		}
  		MltSv <- MSV+1
	}
	
	MltSV = 0

	for(MltSV in seq(1,max(MultiDt$MultiSV))){
  		MltSvSbSt <- subset(MultiDt, MultiDt$MultiSV==MltSV)
  		Strt <- min(MltSvSbSt$MidPos)-StepSize
  		End <- max(MltSvSbSt$MidPos)+StepSize
		Size <- (End - Strt)
		Chr <- as.character(unique(MltSvSbSt$chr))
  
  		if(length(MltSvSbSt$PopLg2) > 1){
    		if(mean(MltSvSbSt$PopLg2) > 0){
      			Lg2 <- max(MltSvSbSt$PopLg2)
    		}else if(mean(MltSvSbSt$PopLg2) < 0){
      			Lg2 <- min(MltSvSbSt$PopLg2)
    		}
  		}else if( length(MltSvSbSt$PopLg2) == 1){
    		Lg2 <- unique(MltSvSbSt$PopLg2)
  		}else if(length(MltSvSbSt$PopLg2) == 0){
    		Lg2 <- 0
  		}
	
		RwNms <- which(MultiDt$MultiSV==MltSV)
		MultiDt[RwNms,'MultiSVlog2'] <- Lg2
		MultiDt[RwNms,'MultiSVSize'] <- Size
  	MultiDt[RwNms,'MultiSVPval'] <- CmptMltPvl((subset(MultiData$MultiSvScl, MultiData$MultiSvScl$chr==Chr & 
                                              MultiData$MultiSvScl$start>=Strt & 
                                              MultiData$MultiSvScl$end<=End)),MultiData$Conf)
	}

	MultiSVData <- as.data.frame(matrix(nrow = max(MultiDt$MultiSV), ncol = 6))
	colnames(MultiSVData) <- c("Chr","Start","End","Size","log2","Pval")
	MltSV = 0

	for(MltSV in seq(max(min(MultiDt$MultiSV),1), max(MultiDt$MultiSV))){
  		MltSvSbSt <- subset(MultiDt, MultiDt$MultiSV==MltSV)

  		if (unique(MltSvSbSt$MultiSVlog2) != "NaN"){
			MultiSVData[MltSV,1] <-  as.character((unique(MltSvSbSt$chr))[1])
    		MultiSVData[MltSV,2] <-  min(MltSvSbSt$MidPos)-StepSize+1
			MultiSVData[MltSV,3] <-  max(MltSvSbSt$MidPos)+StepSize
    		MultiSVData[MltSV,4] <- unique(MltSvSbSt$MultiSVSize)
    		MultiSVData[MltSV,5] <- unique(MltSvSbSt$MultiSVlog2)
    		MultiSVData[MltSV,6] <- unique(MltSvSbSt$MultiSVPval)
  		} else {  
    		MultiSVData[MltSV,4] <-  NA
    		MultiSVData[MltSV,5] <- NA
    		MultiSVData[MltSV,6] <- NA
    	}  
  	}
subset(MultiSVData,MultiSVData$Pval<=0.05&MultiSVData$Size>SVSize)
}





