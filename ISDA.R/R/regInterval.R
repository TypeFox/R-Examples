regInterval <-
function (dados,headersMin,headersMax){
	takeNumColumns <- function(varVectors,Data) {
		name = names(Data)
		numColumns <- 0;
		for (i in 1: length(varVectors)){
			for(j in 1:length(name)){
				if(name[j]== varVectors[i]){
					if(i == 1){
						numColumns = j
					}else{
						numColumns = c(numColumns,j)
					}
				}
			}
		}
		numColumns
	}

	numColumnsMin = takeNumColumns(headersMin,dados )
	numColumnsMax = takeNumColumns(headersMax,dados )
	
	
	matrix_center = matrix(-1,nrow=dim(dados)[1],ncol = length(numColumnsMax))
	for(i in 1:dim(dados)[1]){
		for(j in 1:length(numColumnsMax)){
		minVal = dados[i,numColumnsMin[j]]
		maxVal = dados[i,numColumnsMax[j]]
		
				matrix_center[i,j] = (minVal+maxVal)/2
			
		}
	}
	centerFrame =as.data.frame(matrix_center)
	lm(centerFrame)
}

