#make the treatment contrast matrix, this function allow the input of the treatment contrasts.

	makeTreatProjectors <- function(design.df, trtCols, effectsMatrix, trt.contr, contr.matrix){
	  	#Square contrast matrix with the sepcific contrast (trtContr)
      indMatrix1 <-
      function(x,n, trtContr){
      if(x == 1) X <- trtContr
      else if(x == 2)  X <- mI(n)     
      else       X <- mK(n)
      return(X)
      }

	    #Square contrast matrix without  defining the sepcific contrast
      indMatrix <-
      function(x,n){
      if(x == 1) X <- mJ(n)
      else if(x == 2)  X <- mI(n)     
      else       X <- mK(n)
      return(X)
      }

      	#transfrom each treatment contrast metrix to C matrix
      transContrToT <-
      function(fact, contr){
      T = matrix(0, nrow = nlevels(fact), ncol= nlevels(fact))
      if(is.matrix(contr)){
      rownames(contr) = fact
      for(i in 1:ncol(contr)){         
      x = contr[levels(fact),i]
      T = T + (x %*% t(x))/as.numeric(t(x) %*% x)
      }
      
      }  else {
      names(contr) = fact
      x = contr[levels(fact)]
      T = T + (x %*% t(x))/as.numeric(t(x) %*% x)      
      }                                                     
      return(T)
      }
      
		#obatining the numbers of the levels from the design

		if(any(grepl(":", trtCols))){
		  uniqueTrtCols = unique(unlist(strsplit(trtCols, "\\:")))
			nLevels <- sapply(design.df[,uniqueTrtCols], function(x) nlevels(as.factor(x)))

    }else if(length(trtCols) == 1){ 
    	nLevels = nlevels(design.df[,trtCols])
			names(nLevels) =  trtCols
	  
		}else{
			nLevels <- sapply(design.df[,trtCols], function(x) nlevels(as.factor(x)))
		}

		nLevels = nLevels[rownames(effectsMatrix)]
		nEffects = ncol(effectsMatrix)

		if(all(is.na(trt.contr))){
			#Without the contrasts specifically defined

			X <- as.list(rep(1,nEffects))
			names(X)  <- colnames(effectsMatrix)

			for(i in 1:nrow(effectsMatrix)){
				matList <- lapply(effectsMatrix[i,], function(y) indMatrix(y, nLevels[i]))
				for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]
			}
		}else{
			#With the contrasts specifically defined
			names(trt.contr) = trtCols

			if(!any(sapply(trt.contr, class)=="list")){
				X <- as.list(rep(1,nEffects))
				names(X)  <- colnames(effectsMatrix)
				names(trt.contr)<- colnames(effectsMatrix)

				for(i in 1:nrow(effectsMatrix)){
					if(all( is.na(trt.contr[[i]]))){
						matList <- lapply(effectsMatrix[i,], function(y) indMatrix(y, nLevels[i]))

					} else{

						trtContr = transContrToT(design.df[,rownames(effectsMatrix)[i]],
								trt.contr[[rownames(effectsMatrix)[i]]])
						matList <- lapply(effectsMatrix[i,], function(y)
									indMatrix1(y, nLevels[i], trtContr))
					}
					for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]
				}
			} else {

				if(contr.matrix){
					X <- as.list(rep(1,nEffects))
					names(X)  <- colnames(effectsMatrix)

					trt.contr = lapply(trt.contr, function(x) cbind(sapply(x, function(y) y)))

					names(trt.contr)<- colnames(effectsMatrix)

					for(i in 1:nrow(effectsMatrix)){
						if(all( is.na(trt.contr[[i]]))){
							matList <- lapply(effectsMatrix[i,], function(y) indMatrix(y, nLevels[i]))

						} else{

							trtContr = transContrToT(design.df[,rownames(effectsMatrix)[i]],
									trt.contr[[rownames(effectsMatrix)[i]]])
							matList <- lapply(effectsMatrix[i,], function(y)
										indMatrix1(y, nLevels[i], trtContr))
						}
						for(j in 1:nEffects) X[[j]] <- X[[j]] %x% matList[[j]]
					}
				} else {

					#extract the interactions
					if(any(grepl(":", trtCols))){
						interTerms = grep(":", trtCols)

						for(i in 1:length(interTerms)){
							mainTerms = unlist(strsplit(trtCols[interTerms[i]], ":"))

							mainTerms.list = vector(length = length(mainTerms), mode = "list")
							names(mainTerms.list) =  mainTerms

							for(j in 1:length(mainTerms))
								mainTerms.list[[j]] = names(trt.contr[[mainTerms[j]]])

							#for(j in 1:length(mainTerms))
							# for(k in 1:length(trt.contr[[mainTerms[j]]]))
							#   assign(names(trt.contr[[mainTerms[j]]][k] ), trt.contr[[mainTerms[j]]][[k]])

							trt.contr[[interTerms[i]]] =
									vector( length = nlevels(interaction( mainTerms.list)), mode = "list")

							names(trt.contr[[interTerms[i]]] ) = levels(interaction( mainTerms.list))
						}
					}

					X = list()
					count = 1

					totalLength = length(trt.contr[-which(sapply(trt.contr, class)=="list")]) +
							sum(sapply( trt.contr[which(sapply(trt.contr, class)=="list")], length))

					newEffectsMatrix = matrix(0, ncol =totalLength  , nrow = nrow(effectsMatrix))
					rownames(newEffectsMatrix)  =  rownames(effectsMatrix)

					#constract a new effects matrix to use for get the Kronecker products
					for(i in 1:length(trt.contr)){
						if(is.list(trt.contr[[i]])){
							tmpX = rep(1,length(trt.contr[[i]]))
							names(tmpX) = paste(names(trt.contr[i]), names(trt.contr[[i]]),sep=".")

							tmpCounter =  count:(count+length(trt.contr[[i]])-1)

							X = c(X, tmpX)
							rowNames =  unlist(strsplit(names(trt.contr[i]), ":"))
							newEffectsMatrix[rowNames,tmpCounter] = 1
							count = count+length(trt.contr[[i]])
						}else{
							tmpX = 1
							names(tmpX) = names(trt.contr[i])
							X = c(X, tmpX)

							newEffectsMatrix[,count] = effectsMatrix[,i]

							tmpCounter =  count = count + 1
						}
					}

					colnames(newEffectsMatrix) = names(X)

					for(i in 1:nrow(newEffectsMatrix)){
						if(is.list(trt.contr[[rownames(effectsMatrix)[i]]])){
							trtContr = lapply(trt.contr[[rownames(effectsMatrix)[i]]], function(x)
										transContrToT(design.df[,rownames(effectsMatrix)[i]], x))

							matList <- vector(mode = "list", length = totalLength)

							for( j in 1:ncol(newEffectsMatrix)){

								if(newEffectsMatrix[i,j] == 1){
									trtNames = match(unlist(strsplit(colnames(newEffectsMatrix)[j], "[[:punct:]]")),
											names(trtContr))

									matList[[j]] <- trtContr[[trtNames[-which(is.na(trtNames))]]]
								}else if(newEffectsMatrix[i,j] == 2){
									matList[[j]] <-   mI(nLevels[i])
								}else{
									matList[[j]] <- mK(nLevels[i])
								}
							}
						}else{
							if( all(is.na(trt.contr[[i]]))){
								matList <- lapply(newEffectsMatrix[i,], function(y) indMatrix(y, nLevels[i]))

							} else{
								trtContr = transContrToT(design.df[,rownames(effectsMatrix)[i]],
										trt.contr[[rownames(effectsMatrix)[i]]])
								matList <- lapply(newEffectsMatrix[i,], function(y)
											indMatrix1(y, nLevels[i], trtContr))

							}
						}


						for(j in 1:length(X)) X[[j]] <- X[[j]] %x% matList[[j]]

         }
			 }
			}
			}
			return(X)
		}
