# Get Data
# It takes data as dataframe (edges) or as matrix (table) to be exported in proper form to be used by the diversity function.
# data Data to be processed as dataframe or as matrix. 
# category_row TRUE if column analysis is needed.
get_data <- function(data, category_row=FALSE)
{
	if (is.data.frame(data)) {
		if(category_row==TRUE) {
			#diversity <- data.frame(row.names=levels(data[,2]))
			data <- droplevels(data) #delete un used levels
			X <- matrix(0, nrow=nlevels(data[,2]), ncol=nlevels(data[,1]), dimnames=list(levels(data[,2]),levels(data[,1])))
			X[cbind(data[,2], data[,1])] <- data[,3]
		}
		else {
			#diversity <- data.frame(row.names=levels(data[,1]))
			X <- matrix(0, nrow=nlevels(data[,1]), ncol=nlevels(data[,2]), dimnames=list(levels(data[,1]),levels(data[,2])))
			X[cbind(data[,1], data[,2])] <- data[,3]
		}
	}
	else {
		if (category_row==TRUE) {
			X <- t(data)
			#diversity <- data.frame(row.names=rownames(X))
		}
		else {
			X <- data
			#diversity <- data.frame(row.names=rownames(X))
		}  
	}
	return(X)
	
}