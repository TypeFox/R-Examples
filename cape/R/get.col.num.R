get.col.num <-
function(data.mat, col.which = NULL, warn = TRUE){

	if(is.null(col.which)){
		col.num <- 1:dim(data.mat)[2]
		}else{

	if(is.numeric(col.which[1])){
		return(col.which)
		}else{
	
		#get the locations of the phenotype names
		# col.num <- which(colnames(data.mat) %in% col.which) 
		col.mat <- matrix(col.which, ncol = 1)

		find.col <- function(col.name, data.names){
			col.num <- which(data.names == col.name)
			if(length(col.num) == 0){
				return(0)
				}
			return(col.num)
			}
		
		col.num <- unlist(apply(col.mat, 1, function(x) find.col(x, colnames(data.mat)))) #grep allows us to ignore the case

	
		#if we didn't find all the input phenotypes
		#stop and write out which ones we couldn't find
		didnt.find <- which(col.num == 0)
		if(length(didnt.find) > 0){
			if(warn){
				message("The following column headers are not present:")
				cat(col.which[didnt.find], sep = "\n")
				}
			#remove the 0's from the vector, and return only 
			#the column numbers we found
			col.num <- col.num[-didnt.find]
			}
				
		}
	}
		return(col.num)

}
