

getVcCoefMS.twoPhase =
function(PNTginvATNP, design.df, v.mat, response, table.legend){

   if(all(is.na(response))){
    response = rep(NA, nrow(design.df))
    } 
  
    MS = lapply(PNTginvATNP, function(x) lapply(x, function(q) lapply(q, function(w) t(response) %*%w %*% response)))

    V = v.mat
		VC <-rep("1",length(V)+2)
		names(VC)  = c("DF", c(names(V), "MS"))
		VC = t(data.frame(VC))

	for(i in 1:length(PNTginvATNP)){
			VC = rbind(VC, character(length = length(V) + 2))
			if( names(PNTginvATNP[i]) =="Within"){
				rownames(VC)[nrow(VC)] = paste(names(PNTginvATNP[i]), sep = " ")
			}else{
				rownames(VC)[nrow(VC)] = paste("Between", names(PNTginvATNP[i]), sep = " ")
			}

			for(k in 1:(length(PNTginvATNP[[i]]))){
				tmp <- matrix(0, nrow=length(names(PNTginvATNP[[i]][[k]])) , ncol=(length(V) + 1),
						dimnames=list(names(PNTginvATNP[[i]][[k]]), c(names(V), "MS")))

				for(j in 1:((length(names(PNTginvATNP[[i]][[k]]))))){
					for(z in 1:length(V)){
						tmp[j,z] <- tr(PNTginvATNP[[i]][[k]][[j]] %*% V[[z]])
					}
					tmp[j, (length(V) + 1)] =  as.numeric(MS[[i]][[k]][[j]])
				}

				if(nrow(tmp) == 1 && rownames(tmp) == "Residual"){
					tmp = c(tmp[1], tmp/tmp[1])
					tmp = c(attr(fractions(tmp[-length(tmp)]),"fracs"), as.character(round(tmp[length(tmp)], digits = 5)))
					VC = rbind(VC, tmp)
					if(names(PNTginvATNP[[i]][k]) == "Residual"  && length(PNTginvATNP[[i]]) > 1){
						rownames(VC)[nrow(VC)] =  paste("  ", names(PNTginvATNP[[i]][k]), sep = " ")
					}else if(names(PNTginvATNP[[i]][k]) == "Residual" && length(PNTginvATNP[[i]]) == 1){
						rownames(VC)[nrow(VC)] = paste("Between", names(PNTginvATNP[i]), sep = " ")
						VC = VC[-(nrow(VC)-1),]
					}else{
						rownames(VC)[nrow(VC)] = paste("   Between", names(PNTginvATNP[[i]][k]), sep = " ")
					}
				}else{
					VC = rbind(VC, character(length = length(V)+2))

					if(names(PNTginvATNP[[i]][k]) == "Residual"){
						rownames(VC)[nrow(VC)] =  paste("  ", names(PNTginvATNP[[i]][k]), sep = " ")
					}else{
						rownames(VC)[nrow(VC)] = paste("   Between", names(PNTginvATNP[[i]][k]), sep = " ")
					}
					rownames(tmp) = paste("     ", rownames(tmp), sep = " ")
          tmp = t(apply(tmp, 1, function(x) c(attr(fractions(c(x[1], x[-length(x)]/x[1])),"fracs"),
                    as.character(round(x[length(x)]/x[1], digits = 5)))))

					VC = rbind(VC, tmp)
				}
			}
		}

      if(all(is.na(response))){
      VC = VC[,-ncol(VC)]
    }

		if(length((which(apply(apply(VC[,-1], 2, function(x) x==VC[,"e"]), 2, all))))>1){
			VC = noquote(VC[-1,-2])
		} else{
			VC = noquote(VC[-1,])
		}


		if(table.legend){
			Legend = paste(paste( letters[1:(length(colnames(VC))-1)],colnames(VC)[-1], sep = " = "))
			colnames(VC)[-1] = letters[1:(length(colnames(VC))-1)]
			VC = list(VC = VC, Legend = Legend)
		}
		
		return(VC)
}

