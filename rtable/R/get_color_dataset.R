get_color_dataset = function( color, columns, dataset, default_col = "black" ){
	
	color_data = matrix(default_col, ncol = length( columns ), nrow = nrow(dataset) )
	dimnames(color_data) = list(NULL, columns )
	mode(color_data) = "character"
	
	if( length( color ) > 0 ){
		color.columns = character( length( columns ) )
		names(color.columns) = columns
		color.columns[names(color)] = color
		for(j in names(color) ){
			if( color.columns[j] != "" ){
				color_data[,j] = as.character( dataset[, color.columns[j] ] )
			}
		}
	} 
	
	color_data
}


