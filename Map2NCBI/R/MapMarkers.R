#############################################
### Code to create function "MapMarkers"  ###
### Hanna and Riley                       ###
#############################################


MapMarkers = function(features,markers,nAut,other=c("X"),destfile,savefiles=TRUE){
		if(missing(features)){
			stop("ERROR: Did not specify list of features to use for mapping.")
		}
		if(missing(markers)){
			stop("ERROR: Did not specify list of markers to be mapped.")
		}
		if(missing(nAut)){
			stop("ERROR: Did not specify the number of autosomes present in the marker file.")
		}
		if(savefiles == TRUE){
			if(missing(destfile)){
				stop("ERROR: No path was specified for the folder to save the output file.")
			}else{
			dest = paste(destfile,"MappedMarkers.txt",sep="")
			}
		}
		if(other == FALSE){
			chr = matrix(1:nAut,ncol=1)
			nchr = nAut
		} else{
			Aut = matrix(1:nAut,ncol=1)
			nchr = nAut
			other = matrix(other,ncol=1)
			Oth = matrix(0,dim(other)[1],1)
			for(i in 1:nrow(other)){
				Oth[i,1] = nchr+1
				nchr = nchr+1
			}
			chr = rbind(Aut,Oth)
		}
		chromosome = NULL
    Colnames = matrix(c(colnames(markers),colnames(features),"Distance","Inside?"),nrow=1)
		nCol = ncol(Colnames)
		MarkMap = matrix(0,1,nCol,byrow=TRUE,dimnames=list(c(1),Colnames))
		GLnCol = ncol(features)
		for(i in 1:nchr){
			if(i > nAut){
				j = i-nAut
				k = other[j,1]
			} else { k = i }
			Chr_Features = subset(features, chromosome==k)
			Chr_Markers = subset(markers, chromosome==i)
			nmarkers = nrow(Chr_Markers)
			nfeatures = nrow(Chr_Features)
			if(nmarkers > 0){
				rownames(Chr_Markers) = 1:nmarkers
				rownames(Chr_Features) = 1:nfeatures
				for(locus in 1:nmarkers){
					MarkerInfo = subset(Chr_Markers[locus,])
					MapPos = MarkerInfo[1,'position']
					FeatureInfo = matrix(0,1,GLnCol,byrow=TRUE,dimnames=list(1,colnames(features)))
					Inside = matrix(0,1,2,byrow=TRUE,dimnames=list(1,c("Distance","Inside?")))
					MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
					MaxDis = 1000000
					for(feature in 1:nfeatures){ 
						Start = as.numeric(Chr_Features[feature,'chr_start'])
						Stop = as.numeric(Chr_Features[feature,'chr_stop'])
						DisStart = MapPos-Start
						DisStop = MapPos-Stop
						if(DisStart >= -2500){
							if(DisStop <= 0){
								if(DisStart >= 0){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(0,"Yes,_Inside_Gene")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
								}
								if(DisStart < 0){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStart),"Marker_is_<=_2500_bp_Before_Feature")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
								}
							}
						}
						if(DisStop <= 2500){
							if(DisStart >= 0){
								if(DisStop <= 0){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(0,"Yes,_Inside_Gene")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
								}
								if(DisStop > 0){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStop),"Marker_is_<=_2500_bp_After_Feature")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
								}
							}
						}
					}
					if(MapIt[,'Inside?'] == 0){
						for(feature in 1:nfeatures){
							Start = as.numeric(Chr_Features[feature,'chr_start'])
							Stop = as.numeric(Chr_Features[feature,'chr_stop'])
							DisStart = MapPos-Start
							DisStop = MapPos-Stop				
							if(DisStop <= 5000){
								if(DisStop > 2500){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStop),"Marker_is_>_2500_bp_<=5000_bp_After_Feature")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
								}
							}
							if(DisStart >= -5000){
								if(DisStart < -2500){		
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStart),"Marker_is_>_2500_bp_<=5000_bp_Before_Feature")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
								}
							}
						}
					}
					if(MapIt[,'Inside?'] == 0){			
						for(feature in 1:nfeatures){
							Start = as.numeric(Chr_Features[feature,'chr_start'])
							Stop = as.numeric(Chr_Features[feature,'chr_stop'])
							DisStart = MapPos-Start
							DisStop = MapPos-Stop
							if(DisStart >= -25000){
								if(DisStart < -5000){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStart),"Marker_is_>_5000_bp_<=25000_bp_Before_Feature")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
								}
							}
							if(DisStop <= 25000){
								if(DisStop > 5000){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStop),"Marker_is_>_5000_bp_<=25000_bp_After_Feature")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
								}
							}
						}
					}
					if(MapIt[,'Inside?'] == 0){
						for(feature in 1:nfeatures){
							Start = as.numeric(Chr_Features[feature,'chr_start'])
							Stop = as.numeric(Chr_Features[feature,'chr_stop'])
							DisStart = MapPos-Start
							DisStop = MapPos-Stop
							if(DisStart > (-1*MaxDis)){
								if(DisStart < -25000){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStart),"Nearest_feature_is_>_25,000_bp_after_marker")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
									MaxDis = abs(DisStart)
								}
							}
							if(DisStop < MaxDis){
								if(DisStop > 25000){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStop),"Nearest_feature_is_>_25,000_bp_before_marker")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
									MaxDis = DisStop
								}
							}
						}
					}
					if(MapIt[,'Inside?'] == 0){
						MaxDis = 30000000
						for(feature in 1:nfeatures){
							Start = as.numeric(Chr_Features[feature,'chr_start'])
							Stop = as.numeric(Chr_Features[feature,'chr_stop'])
							DisStart = MapPos-Start
							DisStop = MapPos-Stop
							if(DisStart > (-1*MaxDis)){
								if(DisStart < -1000000){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStart),"Nearest_feature_is_>_1,000,000_bp_after_marker")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
									MaxDis = abs(DisStart)
								}
							}
							if(DisStop < MaxDis){
								if(DisStop > 1000000){
									FeatureInfo = Chr_Features[feature,]
									Inside[1,] = c(abs(DisStop),"Nearest_feature_is_>_1,000,000_bp_before_marker")
									MapIt = cbind(MarkerInfo,FeatureInfo,Inside)
									MaxDis = DisStop
								}
							}
						}
					}
				MarkMap = rbind(MarkMap,MapIt)
				}	
			}else { next }
		} 
		MarkMapF = MarkMap[2:nrow(MarkMap),]
		if(savefiles == TRUE){
			write.table(MarkMapF,dest,quote=FALSE,sep=" ",row.names=FALSE)
		}
		return(MarkMapF)
}