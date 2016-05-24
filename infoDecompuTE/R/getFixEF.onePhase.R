 getFixEF.onePhase  =
 function(effFactors, trt.Coef, T,  Rep, table.legend){

  	trt  = numeric(length(trt.Coef) + ncol(Rep))
		names(trt) = c( names(T), paste("eff",  names(T), sep = ".") )

		for(i in 1:length(effFactors)){
			trt = rbind(trt, character(length = length(T)*2))
			if( names(effFactors[i]) =="Within"){
				rownames(trt)[nrow(trt)] = paste(names(effFactors[i]), sep = " ")
			}else{
				rownames(trt)[nrow(trt)] = paste("Between", names(effFactors[i]), sep = " ")
			}

			for(j in 1:length(effFactors[[i]])){
				if(is.null(effFactors[[i]][[j]])) next
				
				 char.trt = attr(fractions(trt.Coef*unlist(effFactors[[i]][[j]])),"fracs")
				 #char.trt = ifelse(nchar(char.trt) > 8, round(trt.Coef*unlist(effFactors[[i]][[j]]), digits = 7), char.trt)
			   char.trt.eff = attr(fractions(unlist(effFactors[[i]][[j]])),"fracs")
				 #char.trt.eff = ifelse(nchar(char.trt.eff) > 8, round(unlist(effFactors[[i]][[j]]), digits = 7), char.trt.eff)
			   
				trt.temp =c(char.trt, char.trt.eff) 
				
				trt = rbind(trt,trt.temp)
				rownames(trt)[nrow(trt)] =  paste("  ", names(effFactors[[i]][j]), sep = " ")
			}

		}

		trt = trt[-1,]

		trt = noquote(ifelse(trt == "NaN", "", trt))

		if(table.legend){
			Legend = paste(paste( letters[1:(length(colnames(trt)))],colnames(trt), sep = " = "))
			colnames(trt) = letters[1:(length(colnames(trt)))]
			trt = list(trt = trt, Legend = Legend)
		}

    return(trt)

 }
