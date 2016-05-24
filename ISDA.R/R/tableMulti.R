tableMulti <-
function (variableVector,groupedVar,Data){
	getNumColumns <- function(variableVector,Data) {
		name = names(Data)
		numColumn <- 0;
		for (i in 1: length(variableVector)){
			for(j in 1:length(name)){
				if(name[j]== variableVector[i]){
					if(i == 1){
						numColumn = j
					}else{
						numColumn = c(numColumn,j)
					}
				}
			}
		}
		numColumn
	}
	numColumns = getNumColumns(variableVector,Data )
	numColumnsGrouped = getNumColumns(groupedVar,Data )
	originalMatrix= NULL

	for (i in 1:length(numColumns)) {
		if(i == 1){
			originalMatrix <- Data[numColumns[1]]
			
		}else{
			originalMatrix <- cbind(originalMatrix,Data[numColumns[i]])
			
		}
		
	}
	for(j in 1:length(numColumnsGrouped)){
		originalMatrix <- cbind(originalMatrix,Data[numColumnsGrouped[j]])
	}

	
#	iremos pecorrer o vetor armazenando noutra tabela cada caso diferente e a cada linha
#	da tabela original iremos comparar com a nova para verificar se ela ja existe
#	caso exista adicionaremos um ao contador

	contTEMP = 1
	resultMatrix = NULL
	qntColumns =length(variableVector)
	qntColumnsGrouped = length(groupedVar)
	for (i in 1:dim(originalMatrix)[1]) {
		if(i == 1){
			for(z in 1:qntColumns){
				if(z == 1){
					resultMatrix = originalMatrix[1,1]
				}else{
					resultMatrix = cbind(resultMatrix,originalMatrix[1,z])
				}
			}
			a=qntColumns+1
			b=qntColumns+qntColumnsGrouped
			for(w in a:b){
				resultMatrix = cbind(resultMatrix,originalMatrix[1,w],originalMatrix[1,w])
			}

		}else{

			found = FALSE
			itsEqual = NULL
			newLineCounter =1

			while( (found == FALSE) && (newLineCounter<= dim(resultMatrix)[1])){

				
				itsEqual = TRUE
				columnsCounter = 1

				while(itsEqual && columnsCounter<= (qntColumns)  ){
		
					if(originalMatrix[i,columnsCounter] != resultMatrix[newLineCounter,columnsCounter] ){
						
						itsEqual =FALSE
					
					}
					columnsCounter = columnsCounter +1
				}
				if(itsEqual == TRUE){
					a = qntColumns +1
					b=qntColumns+qntColumnsGrouped
				
					for (z in a:b) {
						if(originalMatrix[i,z]< resultMatrix[newLineCounter,(2*z)-(qntColumns+1)]){
							resultMatrix[newLineCounter,(2*z)-(qntColumns+1)] = originalMatrix[i,z]
						}
						if(originalMatrix[i,z]>resultMatrix[newLineCounter,(2*z)-(qntColumns)]){
							resultMatrix[newLineCounter,(2*z)-(qntColumns)] = originalMatrix[i,z]
						}
					}
				
					found = TRUE
				}else if(newLineCounter == dim(resultMatrix)[1]){
					a=qntColumns+1
					b=qntColumns+qntColumnsGrouped
					lineTemp = NULL
					for(z in 1: qntColumns){
						if(z == 1 ){
							lineTemp = originalMatrix[i,z]
						}else{
							lineTemp = cbind(lineTemp, originalMatrix[i,z])
						}
					}
					for(w in a:b){
						lineTemp = cbind(lineTemp,originalMatrix[i,w],originalMatrix[i,w])
					}
				
						resultMatrix=rbind(resultMatrix,lineTemp) 
				
					
					found = TRUE
				}
				newLineCounter = newLineCounter+1
			}
			
		}
	}
	groupNames = NULL
	for(i in 1:length(groupedVar)){
		nameMin = paste(groupedVar[i],"min")
		nameMax = paste(groupedVar[i],"max")
		groupNames = c(groupNames,nameMin,nameMax)
	}
	colnames(resultMatrix) = c(variableVector,groupNames)
	resultMatrix
}

