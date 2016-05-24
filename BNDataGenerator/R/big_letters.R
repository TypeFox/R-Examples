big_letters = function(len)
{
	letters_list = list()
	letters_list[[1]] = letters
	
	while(TRUE)
	{
		num_of_letters = 0;
		len_letters_list = length(letters_list)
		for (i in 1:len_letters_list)
		{
			num_of_letters = num_of_letters + length(letters_list[[i]]);
		}
		
		if (num_of_letters < len)
		{
			merge_mat = merge(letters, letters_list[[len_letters_list]])
			letters_list[[len_letters_list+1]] = sort(paste(merge_mat[,1], merge_mat[,2], sep=""))
		} else {
			break;
		}
	}
	
	result = NULL;
	
	for(i in 1:length(letters_list))
	{
		result = c(result, letters_list[[i]])
	}
	result = result
	
	return(result[1:len])
}