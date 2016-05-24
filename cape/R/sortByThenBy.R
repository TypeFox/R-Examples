sortByThenBy <-
function(tableX, sort.cols = c(1,2), col.type = c("c", "n"), decreasing = FALSE, return.order = FALSE){
	
	if(length(sort.cols) > 2){
		stop("This script can't handle more than 2 columns")
		}
	
	if(length(sort.cols) != length(col.type)){
		stop("The col.type vector must be the same length as sort.cols")
		}
	
	if(length(sort.cols) == 1){
		
		if(col.type == "n"){
			sorted.col <- sort.int(as.numeric(tableX[,sort.cols]), index.return = TRUE, decreasing = decreasing)
			}else{
				sorted.col <- sort.int(tableX[,sort.cols], index.return = TRUE, decreasing = decreasing)
				}
		
		if(return.order){
			return(sorted.col$ix)
		}else{
			new.table <- tableX[sorted.col$ix,]
			return(new.table)
			}

		}else{
			col.order <- matrix(NA, nrow = nrow(tableX), ncol = ncol(tableX))
			
			#start with the table ordered by the first sort column. We will change
			#chunks of this as we go
			if(col.type[1] == "n"){
				new.order <- sort.int(as.numeric(tableX[,sort.cols[1]]), index.return = TRUE, decreasing = decreasing)$ix	
				}else{
					new.order <- sort.int(tableX[,sort.cols[1]], index.return = TRUE, decreasing = decreasing)$ix
					}
			sorted.table <- tableX[new.order, ]
			col.order[,1] <- new.order
				
			if(col.type[1] == "n"){
				u_col_el <- sort(as.numeric(unique(sorted.table[,sort.cols[1]]))) #find the unique column elements for column 1
			}else{
				u_col_el <- sort(unique(sorted.table[,sort.cols[1]])) #find the unique column elements for column 1
				}
				
			
			new.order <- NULL
				
			for(i in 1:length(u_col_el)){ #go through each of these elements, and sort the second column within the category of the first column
				el.locale <- which(sorted.table[,sort.cols[1]] == u_col_el[i]) #find the entries for element i in column 1
				if(col.type[2] == "n"){
					subset.order <- sort.int(as.numeric(sorted.table[el.locale, sort.cols[2]]), index.return = TRUE, decreasing = decreasing)$ix
					}else{
					subset.order <- sort.int(sorted.table[el.locale, sort.cols[2]], index.return = TRUE, decreasing = decreasing)$ix
					}
				new.order <- c(new.order, el.locale[subset.order])
				}
			
			col.order[,2] <- new.order
					
			if(return.order){
				return(col.order)
				}else{
				sorted.table <- sorted.table[new.order,]			
				return(sorted.table)
				}
		} #end case for if sort.cols has more than one element
}
