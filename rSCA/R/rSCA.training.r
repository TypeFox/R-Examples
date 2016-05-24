#################################################################
# Filename: 	rSCA.training.r
# Created: 		2010/10/29, Regina, SK, Canada
# Author: 		Xiuquan Wang
# Email:		xiuquan.wang@gmail.com
# Desp: 		Implementation of stepwise cluster analysis
# ===============================================================
# History: 	2010/10/29		created by Xiuquan Wang
#	     	2010/11/10   	modified the way to find X[i,j], by Xiuquan Wang
#	     	2010/11/19   	replaced %in% with two-level order, by Xiuquan Wang
#	     	2010/11/20   	added hash table to avoid duplicated comparisons, by Xiuquan Wang
#	     	2010/11/28   	added optimization for the computation of wilks value, by Xiuquan Wang
#			2014/04/25		added the option of Golden Section Search, by Xiuquan Wang
##################################################################

# ---------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------

#: generate randome filename
get_random_filename = function()
{
	cur_time = Sys.time()
	s_temp = format(cur_time, "%Y%m%d_%H%M%S")
	options(digits=13)
	n_milliseconds = as.numeric(unclass(Sys.time())) * 1000
	n_randId = floor(runif(1, 1000, 9000))
	s_temp = paste(s_temp, "_", n_milliseconds, "_", n_randId, sep = "")
	return(s_temp)
}

#: calculate radius for a leaf node
#: input: column vector from y data
f_cal_radius = function(o_colvec)
{
	n_radius = (max(o_colvec) - min(o_colvec)) / 2
	return(n_radius)
}

#: hash storing function
#: input: a, b, value
f_hash_store = function(a, b, value)
{
	n_rowId = (a + b) %% rSCA.env$n_sample_size + 1
	if (rSCA.env$o_hashtable_matrix[n_rowId, 1] <= 0)
	{
		rSCA.env$o_hashtable_matrix[n_rowId, 1] <- value
		rSCA.env$o_hashtable_matrix[n_rowId, 2] <- a
		rSCA.env$o_hashtable_matrix[n_rowId, 3] <- b
	}
	else
	{
		#: store into hashtable list
		n_steps = 0
		n_endId = n_rowId
		n_pointer = rSCA.env$o_hashtable_matrix[n_endId, 4]
		while (n_pointer > 0)
		{
			n_steps = n_steps + 1
			n_endId = n_pointer
			n_pointer = rSCA.env$o_hashtable_list[[n_endId]]$pointer
		}
	
		#: set the cursor of hashstable list
		rSCA.env$n_hashtable_list_index <- rSCA.env$n_hashtable_list_index + 1
		o_temp_list = list(n_value=value, n_a=a, n_b=b, pointer=0)
		rSCA.env$o_hashtable_list[[rSCA.env$n_hashtable_list_index]] <- o_temp_list

		#: update the previous pointer
		if (n_steps < 1)
		{
			#: only one replicate
			rSCA.env$o_hashtable_matrix[n_endId, 4] <- rSCA.env$n_hashtable_list_index
		}
		else
		{
			rSCA.env$o_hashtable_list[[n_endId]]$pointer <- rSCA.env$n_hashtable_list_index
		}
	}
}

#: hash reading function
#: input: a, b
#: output: ab-value
f_hash_read = function(a, b)
{
	n_value = 0
	n_rowIndex = (a + b) %% rSCA.env$n_sample_size + 1
	n_a = rSCA.env$o_hashtable_matrix[n_rowIndex, 2]
	n_b = rSCA.env$o_hashtable_matrix[n_rowIndex, 3]
	if (n_a == a && n_b == b)
	{
		n_value = as.numeric(rSCA.env$o_hashtable_matrix[n_rowIndex, 1])
	}
	else
	{
		n_pointer = as.numeric(rSCA.env$o_hashtable_matrix[n_rowIndex, 4])
		while(n_pointer > 0)
		{
			n_a = as.numeric(rSCA.env$o_hashtable_list[[n_pointer]]$n_a)
			n_b = as.numeric(rSCA.env$o_hashtable_list[[n_pointer]]$n_b)
			if (n_a == a && n_b == b)
			{
				n_value = as.numeric(rSCA.env$o_hashtable_list[[n_pointer]]$n_value)
				break;
			}
			n_pointer = as.numeric(rSCA.env$o_hashtable_list[[n_pointer]]$pointer)
		}
	}
	return(n_value)
}


#: construct the ordered sub-matrix 
#: input: rowids matrix = [ , 1]
#: output: ordered sub-matrix
f_ordered_submatrix = function(o_rowids_matrix)
{
	#: sub matrix to be returned
	o_ordered_submatrix = matrix(0, nrow(o_rowids_matrix), rSCA.env$n_sample_x_cols + 1)
	o_ordered_submatrix[ , 1] = o_rowids_matrix
	
	for (icol in 2:(rSCA.env$n_sample_x_cols + 1))
	{
		irow = 0
		for (iindex in c(o_rowids_matrix))
		{
			irow = irow + 1
			o_ordered_submatrix[irow, icol] = rSCA.env$o_sorted_matrix[iindex, (icol-1)]
		}
	}
	return(o_ordered_submatrix)
}

#: calculate wilks statistic value
#: input: top matrix & bot matrix
#: output: wilks value
f_wilks_statistic <- function(o_top_matrix, o_bot_matrix)
{
	n_top = nrow(o_top_matrix)
	n_bot = nrow(o_bot_matrix)
	if ((n_top + n_bot) <= (rSCA.env$n_sample_y_cols + 1))
	{
		return(0)
	}

	o_top_mean = matrix(colMeans(o_top_matrix), 1, )
	o_bot_mean = matrix(colMeans(o_bot_matrix), 1, )
	o_between_matrix = (n_top*n_bot)/(n_top+n_bot) * crossprod(o_top_mean - o_bot_mean, o_top_mean - o_bot_mean)
	
	o_top_matrix_mean = matrix(colMeans(o_top_matrix), nrow(o_top_matrix), ncol(o_top_matrix), byrow = TRUE)
	o_bot_matrix_mean = matrix(colMeans(o_bot_matrix), nrow(o_bot_matrix), ncol(o_bot_matrix), byrow = TRUE)
	o_within_top = crossprod(o_top_matrix - o_top_matrix_mean, o_top_matrix - o_top_matrix_mean)
	o_within_bot = crossprod(o_bot_matrix - o_bot_matrix_mean, o_bot_matrix - o_bot_matrix_mean)
	o_within_matrix = o_within_top + o_within_bot	
	
	#: wilks value should be in [0, 1] ==> to assure that F value is in [0, +Inf]
	#: however, wilks value could be negative or even -Inf (due to n_det_total = 0)
	o_wilks_vaule = 0
	n_det_within = det(o_within_matrix)
	n_det_total = det(o_within_matrix + o_between_matrix)
	if (n_det_total == 0)
	{
		if (n_det_within < 0)
			o_wilks_vaule = -1 #: -Inf
		else if (n_det_within > 0)
			o_wilks_vaule = 1 #: +Inf
		else
			o_wilks_vaule = 0
	}
	else
		o_wilks_vaule = n_det_within / n_det_total
	
	return(o_wilks_vaule)
}

#: calculate minimum wilks value for a matrix indentified by a group of row id [original]
#: input: matrix of rowid [nrows, col=1]
#: output: min_wilks_list(min_wilks_value=1, col_id=1, x_value=1, left_rowids=matrix(), right_rowids=matrix())
f_min_wilks <- function(o_matrix_rowid)
{
	if (nrow(o_matrix_rowid) <= (rSCA.env$n_sample_y_cols + 1))
	{
		return(list(min_wilks_value=0, col_id=0, x_value=0, left_rowids=matrix(), right_rowids=matrix()))
	}

	#: matrix to store all wilks value and the corresponding colid and rowid
	#: element format: wilks value || colid || x_value[get by original row id] || sub row id
	#: NOTE: divide the o_matrix_rowid according to sub_row_id!!!
	o_wilks_value_matrix = matrix(0, (nrow(o_matrix_rowid)-1)*rSCA.env$n_sample_x_cols, 4)
	n_wilks_value_matrix_index = 0 

	#: get the ordered sub-matrix
	o_ordered_sub_matrix = f_ordered_submatrix(o_matrix_rowid)

	#: calculate all wilks values
	for (icol in 1:rSCA.env$n_sample_x_cols)
	{
		#: Step-1: contruct the temporary matrix to store the sub sorted matrix
		o_temp_sub_matrix = matrix(o_ordered_sub_matrix[order(o_ordered_sub_matrix[, (icol+1)]), 1], , 1)

		#: Step-2: calculate all wilks value for this column
		for (irow in 1:(nrow(o_matrix_rowid)-1))
		{
			#: o_top_matrix: top y matrix divided by original row id
			#: o_bot_matrix: bot y matrix divided by original row id
			o_top_matrix = matrix(0, irow, rSCA.env$n_sample_y_cols)
			o_top_matrix[1:irow, ] = data.matrix(rSCA.env$o_sample_data_y[c(o_temp_sub_matrix[1:irow,1]), ])

			o_bot_matrix = matrix(0, nrow(o_matrix_rowid)-irow, rSCA.env$n_sample_y_cols)
			o_bot_matrix[1:(nrow(o_matrix_rowid)-irow), ] = data.matrix(rSCA.env$o_sample_data_y[c(o_temp_sub_matrix[(irow+1):nrow(o_matrix_rowid),1]), ])

			n_wilks_value_matrix_index = n_wilks_value_matrix_index + 1
			o_wilks_value_matrix[n_wilks_value_matrix_index, 1] = f_wilks_statistic(o_top_matrix, o_bot_matrix)
			o_wilks_value_matrix[n_wilks_value_matrix_index, 2] = icol
			o_wilks_value_matrix[n_wilks_value_matrix_index, 3] = rSCA.env$o_sample_data_x[o_temp_sub_matrix[irow, 1], icol]
			o_wilks_value_matrix[n_wilks_value_matrix_index, 4] = irow

			#: clear memory
			rm(o_top_matrix)
			rm(o_bot_matrix)
		}
		#: clear memory
		rm(o_temp_sub_matrix)
	}

	#: find the min wilks value by sorting in terms of col 1
	o_wilks_value_matrix = o_wilks_value_matrix[order(o_wilks_value_matrix[,1]),]
	o_min_wilks_list = list(min_wilks_value=0, col_id=0, x_value=0, left_rowids=matrix(), right_rowids=matrix())
	o_min_wilks_list$min_wilks_value = o_wilks_value_matrix[1, 1]
	o_min_wilks_list$col_id = o_wilks_value_matrix[1, 2]
	o_min_wilks_list$x_value = o_wilks_value_matrix[1, 3]	

	#: construct the left and right rowids matrix
	o_temp_rowids_matrix = matrix(o_ordered_sub_matrix[order(o_ordered_sub_matrix[, (o_wilks_value_matrix[1, 2]+1)]), 1], , 1)

	o_min_wilks_list$left_rowids = matrix(o_temp_rowids_matrix[1:o_wilks_value_matrix[1, 4], 1], , 1)
	o_min_wilks_list$right_rowids = matrix(o_temp_rowids_matrix[(o_wilks_value_matrix[1, 4]+1):nrow(o_matrix_rowid), 1], , 1)

	#: clear memory
	rm(o_wilks_value_matrix)
	
	#: return value
	return(o_min_wilks_list)
}

#: Using the Golden Section Search method
#: calculate minimum wilks value for a matrix indentified by a group of row id [original]
#: input: matrix of rowid [nrows, col=1]
#: output: min_wilks_list(min_wilks_value=1, col_id=1, x_value=1, left_rowids=matrix(), right_rowids=matrix())
f_min_wilks_GSS <- function(o_matrix_rowid)
{
	if (nrow(o_matrix_rowid) <= (rSCA.env$n_sample_y_cols + 1))
	{
		return(list(min_wilks_value=0, col_id=0, x_value=0, left_rowids=matrix(), right_rowids=matrix()))
	}

	#: matrix to store all wilks value and the corresponding colid and rowid
	#: element format: wilks value || colid || x_value[get by original row id] || sub row id
	#: NOTE: divide the o_matrix_rowid according to sub_row_id!!!
	o_wilks_value_matrix = matrix(0, rSCA.env$n_sample_x_cols, 4)
	n_wilks_value_matrix_index = 0 

	#: get the ordered sub-matrix
	o_ordered_sub_matrix = f_ordered_submatrix(o_matrix_rowid)

	#: calculate all wilks values
	for (icol in 1:rSCA.env$n_sample_x_cols)
	{
		#: Step-1: contruct the temporary matrix to store the sub sorted matrix
		o_temp_sub_matrix = matrix(o_ordered_sub_matrix[order(o_ordered_sub_matrix[, (icol+1)]), 1], , 1)

		#: Step-2: find the location of minimum wilks value in this column
		n_min_wilks_value = 0
		n_a = 1
		n_b = nrow(o_matrix_rowid)
		n_k1 = floor(0.5 + n_a + 0.618 * (n_b - n_a)) #: right-side point
		n_k2 = floor(0.5 + n_a + 0.382 * (n_b - n_a)) #: left-side point
		while(n_k1 > n_k2)
		{
			#: 1> cut at n_k2 ==> left-side point
			n_top_rows_LS = n_k2 - n_a + 1
			o_top_matrix_LS = matrix(0, n_top_rows_LS, rSCA.env$n_sample_y_cols)
			o_top_matrix_LS[1:n_top_rows_LS, ] = data.matrix(rSCA.env$o_sample_data_y[c(o_temp_sub_matrix[n_a:n_k2,1]), ])
			n_bot_rows_LS = n_b - n_k2
			o_bot_matrix_LS = matrix(0, n_bot_rows_LS, rSCA.env$n_sample_y_cols)
			o_bot_matrix_LS[1:n_bot_rows_LS, ] = data.matrix(rSCA.env$o_sample_data_y[c(o_temp_sub_matrix[(n_k2+1):n_b,1]), ])
			
			n_wilks_LS = f_wilks_statistic(o_top_matrix_LS, o_bot_matrix_LS)
			
			#: update the minimum wilks value
			n_min_wilks_value = n_wilks_LS
			
			#: 2> cut at n_k1 ==> right-side point
			n_top_rows_RS = n_k1 - n_a + 1
			o_top_matrix_RS = matrix(0, n_top_rows_RS, rSCA.env$n_sample_y_cols)
			o_top_matrix_RS[1:n_top_rows_RS, ] = data.matrix(rSCA.env$o_sample_data_y[c(o_temp_sub_matrix[n_a:n_k1,1]), ])
			n_bot_rows_RS = n_b - n_k1
			o_bot_matrix_RS = matrix(0, n_bot_rows_RS, rSCA.env$n_sample_y_cols)
			o_bot_matrix_RS[1:n_bot_rows_RS, ] = data.matrix(rSCA.env$o_sample_data_y[c(o_temp_sub_matrix[(n_k1+1):n_b,1]), ])
			
			n_wilks_RS = f_wilks_statistic(o_top_matrix_RS, o_bot_matrix_RS)
			
			#: compare
			if (n_wilks_LS >= n_wilks_RS)
			{
				#: the right-side point is a good one, so focusing on the right-side
				n_a = n_k2
			}
			else
			{
				#: otherwise, focusing on the left-side
				n_b = n_k1
			}
			
			#: update k1 + k2
			n_k1 = floor(0.5 + n_a + 0.618 * (n_b - n_a)) #: right-side point
			n_k2 = floor(0.5 + n_a + 0.382 * (n_b - n_a)) #: left-side point

			#: clear memory
			rm(o_top_matrix_LS)
			rm(o_bot_matrix_LS)
			rm(o_top_matrix_RS)
			rm(o_bot_matrix_RS)
		}
		
		#: regard n_k1 as the best potential cutting point
		n_wilks_value_matrix_index = n_wilks_value_matrix_index + 1
		o_wilks_value_matrix[n_wilks_value_matrix_index, 1] = n_min_wilks_value
		o_wilks_value_matrix[n_wilks_value_matrix_index, 2] = icol
		o_wilks_value_matrix[n_wilks_value_matrix_index, 3] = rSCA.env$o_sample_data_x[o_temp_sub_matrix[n_k1, 1], icol]
		o_wilks_value_matrix[n_wilks_value_matrix_index, 4] = n_k1
			
		#: clear memory
		rm(o_temp_sub_matrix)
	}

	#: find the min wilks value by sorting in terms of col 1
	o_wilks_value_matrix = o_wilks_value_matrix[order(o_wilks_value_matrix[,1]),]
	o_min_wilks_list = list(min_wilks_value=0, col_id=0, x_value=0, left_rowids=matrix(), right_rowids=matrix())
	o_min_wilks_list$min_wilks_value = o_wilks_value_matrix[1, 1]
	o_min_wilks_list$col_id = o_wilks_value_matrix[1, 2]
	o_min_wilks_list$x_value = o_wilks_value_matrix[1, 3]	

	#: construct the left and right rowids matrix
	o_temp_rowids_matrix = matrix(o_ordered_sub_matrix[order(o_ordered_sub_matrix[, (o_wilks_value_matrix[1, 2]+1)]), 1], , 1)

	o_min_wilks_list$left_rowids = matrix(o_temp_rowids_matrix[1:o_wilks_value_matrix[1, 4], 1], , 1)
	o_min_wilks_list$right_rowids = matrix(o_temp_rowids_matrix[(o_wilks_value_matrix[1, 4]+1):nrow(o_matrix_rowid), 1], , 1)

	#: clear memory
	rm(o_wilks_value_matrix)

	#: return value
	return(o_min_wilks_list)
}

#: calculate f statistic value and check if it can be divided or cut
#: input: min_wilks_list
#: output: 0 -> no, 1 -> yes
#: NOTE: F statistic -> F(o_df_numerator, o_df_dominator)
f_cal_chk_f = function(min_wilks_list)
{
	n_row_left = nrow(min_wilks_list$left_rowids)
	n_row_right = nrow(min_wilks_list$right_rowids)
	if ((n_row_left + n_row_right) <= (rSCA.env$n_sample_y_cols + 1))
	{
		return(0)
	}
	if (is.na(as.numeric(min_wilks_list$min_wilks_value)))
	{
		#: Exception -> wilks_value is NaN in the following case, e.g.
		#: top matrix: 2 samples
		#: [a1,b1,c1,d1,e1]
		#: [a1,b1,c1,d1,e1]
		#: bot matrix: 5 samples
		#: [a2,b2,c2,d2,e2]
		#: [a2,b2,c2,d2,e2]
		#: [a2,b2,c2,d2,e2]
		#: [a2,b2,c2,d2,e2]
		#: [a2,b2,c2,d2,e2]
		return(1)
	}
	if (as.numeric(min_wilks_list$min_wilks_value) == 0)
	{
		return(1)
	}
	if (as.numeric(min_wilks_list$min_wilks_value) < 0 || as.numeric(min_wilks_list$min_wilks_value) > 1)
	{
		return(1)
	}
	o_df_numerator = rSCA.env$n_sample_y_cols
	o_df_dominator = n_row_left + n_row_right - rSCA.env$n_sample_y_cols - 1
	
	#: calcuate f value
	o_f_value = ((1 - min_wilks_list$min_wilks_value) / min_wilks_list$min_wilks_value) * (o_df_dominator / o_df_numerator)
	#: find the f value to compare
	o_f_criterion = qf((1-rSCA.env$n_alpha), o_df_numerator, o_df_dominator)
	
	#: flag: 0=can not cut, 1=cut
	o_check_flag = 0
	if (o_f_value >= o_f_criterion) o_check_flag = 1
	return(o_check_flag)
}

#: add or update the processing result into out tree
#: input: flag -> 12:cut, 10:leaf, 21:merge
#:	    parent_id -> 0:root, >0:others
#:	    merge_id -> 0:no merge, >0:to be merge [another is parent_id]
#:	    min_wilks_list -> calculated minimum wilks list, [see Function:f_min_wilks]
#: output: o_updatedindex_vector = c(id,id), at most 2 nodes
#:	     merge[21]-> c(n_cursor_tree, 0): new node
#:	     leaf[10]-> c(parent_id, 0)
#:	     cut[12]-> c(left_id, right_id)
#: NOTE: the assignment operator to rSCA.env$o_output_tree needs use '<-'
f_addorupdate_output_tree = function(flag, parent_id, merge_id, min_wilks_list)
{
	#: define the updated index vector
	o_updatedindex_vector = c(0, 0)

	#: output tree list structure
	#: list(id=1, col_index=1, value=1, left=1, right=1, rowids_matrix=matrix())

	if (flag == 21)
	{
		#: merge

		#: set the cursor for output tree
		n_cursor_tree = length(rSCA.env$o_output_tree) + 1

		#: min_wilks_list(min_wilks_value=1, col_id=1, x_value=1, left_rowids=matrix(), right_rowids=matrix())
		o_temp_list = list(id=n_cursor_tree, col_index=min_wilks_list$col_id, value=min_wilks_list$x_value, left=0, right=0, rowids_matrix=rbind(min_wilks_list$left_rowids, min_wilks_list$right_rowids))
		rSCA.env$o_output_tree[[n_cursor_tree]] <- o_temp_list

		rSCA.env$o_output_tree[[parent_id]]$left <- n_cursor_tree
		rSCA.env$o_output_tree[[parent_id]]$right <- n_cursor_tree

		rSCA.env$o_output_tree[[merge_id]]$left <- n_cursor_tree
		rSCA.env$o_output_tree[[merge_id]]$right <- n_cursor_tree

		o_updatedindex_vector[1] = n_cursor_tree

		rm(o_temp_list)
	}
	if (flag == 10)
	{
		#: leaf
		rSCA.env$o_output_tree[[parent_id]]$left <- -1
		rSCA.env$o_output_tree[[parent_id]]$right <- -1
		#: be processed, but still can be merged
		o_updatedindex_vector[1] = parent_id
	}
	if (flag == 12)
	{
		#: cut
		#: set the cursor for output tree
		n_cursor_tree = length(rSCA.env$o_output_tree) + 1

		#: update cut info in the parent node
		rSCA.env$o_output_tree[[parent_id]]$col_index <- min_wilks_list$col_id
		rSCA.env$o_output_tree[[parent_id]]$value <- min_wilks_list$x_value

		#: add left sub node
		o_temp_list_left = list(id=n_cursor_tree, col_index=0, value=0, left=0, right=0, rowids_matrix=min_wilks_list$left_rowids)
		rSCA.env$o_output_tree[[n_cursor_tree]] <- o_temp_list_left

		#: update parent node -> left
		rSCA.env$o_output_tree[[parent_id]]$left <- n_cursor_tree
		o_updatedindex_vector[1] = n_cursor_tree

		#: add right sub node
		n_cursor_tree = length(rSCA.env$o_output_tree) + 1
		o_temp_list_right = list(id=n_cursor_tree, col_index=0, value=0, left=0, right=0, rowids_matrix=min_wilks_list$right_rowids)
		rSCA.env$o_output_tree[[n_cursor_tree]] <- o_temp_list_right

		#: update parent node -> right
		rSCA.env$o_output_tree[[parent_id]]$right <- n_cursor_tree
		o_updatedindex_vector[2] = n_cursor_tree

		#: store left and right subnodes into hash table
		f_hash_store(n_cursor_tree - 1, n_cursor_tree, 1)

		rm(o_temp_list_left)
		rm(o_temp_list_right)
	}

	#: return the updated index vector
	return(o_updatedindex_vector)
}

#: set node in output tree as a leaf
#: input: node_id 
#: output: o_updatedindex_vector [see Function:f_addorupdate_output_tree]
f_setnode_asleaf = function(node_id)
{
	f_addorupdate_output_tree(flag=10, parent_id=node_id, merge_id=0, min_wilks_list=list(0))
}

#: cut node
#: input: node_id, min_wilks_list
#: output: o_updatedindex_vector [see Function:f_addorupdate_output_tree]
f_cutnode = function(node_id, min_wilks_list)
{
	f_addorupdate_output_tree(flag=12, parent_id=node_id, merge_id=0, min_wilks_list)
}

#: merge nodes
#: input: node_id, merge_id, min_wilks_list
#: output: o_updatedindex_vector [see Function:f_addorupdate_output_tree]
f_mergenodes = function(node_id, merge_id, min_wilks_list)
{
	f_addorupdate_output_tree(flag=21, parent_id=node_id, merge_id, min_wilks_list)
}

#: check if node is a leaf
#: input: node id
#: output: 1:yes, 0:no
f_checkif_leaf = function(n_nodeid)
{

	#: NOTE: the minimum sample size for each node is = cols of y + 1
	#: see the df of F(p, n1+n2-p-1)
	#: that means a node is a leaf when it's sample size <= cols of y + 1
	if (nrow(rSCA.env$o_output_tree[[n_nodeid]]$rowids_matrix) <= (rSCA.env$n_sample_y_cols + 1))
	{
		#: set as a leaf
		o_updatedindex_vector = f_setnode_asleaf(n_nodeid)
		return(1)
	}
	if (rSCA.env$o_output_tree[[n_nodeid]]$left == -1 && rSCA.env$o_output_tree[[n_nodeid]]$right == -1)
	{
		return(1)
	}
	return(0)
}

#: process node
#: input: node_id
#: NOTE: first process left sub node, then process right sub node
f_processnode = function(node_id)
{
	#: if a leaf, just exit this function
	if (f_checkif_leaf(node_id) == 1)
	{
		return(0)
	}
	#: add this node into stack
	rSCA.env$o_nodeid_stack_cut[rSCA.env$n_nodeid_statck_cut_cursor] <- node_id
	rSCA.env$n_nodeid_statck_cut_cursor <- rSCA.env$n_nodeid_statck_cut_cursor + 1

	#: loop counter
	n_loop_counter = 0

	#: initiate the loop
	rSCA.env$n_flag_cut <- 1
	while(rSCA.env$n_flag_cut == 1)
	{
		n_loop_counter = n_loop_counter + 1
		if (rSCA.env$b_debug)
			cat("Loop:\t\t\t[", n_loop_counter, "]\r\n------------------------------------------\r\n", file = rSCA.env$s_logfilepath, sep = " ", append = TRUE)

		#: do while only if cutting occured
		
		#: clear cut flag
		rSCA.env$n_flag_cut <- 0

		#: do cut till the cut stack is empty
		while(rSCA.env$n_nodeid_statck_cut_cursor > 1)
		{
			#: read a node from stack
			n_nodeid_cut_temp = rSCA.env$o_nodeid_stack_cut[rSCA.env$n_nodeid_statck_cut_cursor - 1]

			#: if this node has been cutted and set as a leaf, continue the next node
			if (f_checkif_leaf(n_nodeid_cut_temp) >= 1)
			{
				#cat("Cutting node_id:", n_nodeid_cut_temp, " ->> NO NEED!");cat("\n")

				#: store into merge stack
				rSCA.env$o_nodeid_stack_merge[rSCA.env$n_nodeid_statck_merge_cursor] <- n_nodeid_cut_temp
				rSCA.env$n_nodeid_statck_merge_cursor <- rSCA.env$n_nodeid_statck_merge_cursor + 1

				rSCA.env$n_nodeid_statck_cut_cursor <- rSCA.env$n_nodeid_statck_cut_cursor - 1
				next
			}

			#: calculate minimum wilks value
			min_wilks_list = list()
			if (rSCA.env$b_GSS)
				min_wilks_list = f_min_wilks_GSS(rSCA.env$o_output_tree[[n_nodeid_cut_temp]]$rowids_matrix)
			else
				min_wilks_list = f_min_wilks(rSCA.env$o_output_tree[[n_nodeid_cut_temp]]$rowids_matrix)

			#: check if this node can be cut
			n_cut_flag = f_cal_chk_f(min_wilks_list)
			if (n_cut_flag == 1)
			{
				#: if can be cut, then cut it and process it's sub nodes, respectively
				o_updatedindex_vector = f_cutnode(n_nodeid_cut_temp, min_wilks_list)
				rSCA.env$n_cut_times <- rSCA.env$n_cut_times + 1

				#: delete this node (replace) and store the 2 sub nodes into cut stack
				rSCA.env$o_nodeid_stack_cut[rSCA.env$n_nodeid_statck_cut_cursor - 1] <- o_updatedindex_vector[1]
				rSCA.env$o_nodeid_stack_cut[rSCA.env$n_nodeid_statck_cut_cursor] <- o_updatedindex_vector[2]
				rSCA.env$n_nodeid_statck_cut_cursor <- rSCA.env$n_nodeid_statck_cut_cursor + 1

				#: update cut flag
				rSCA.env$n_flag_cut <- 1

				if (rSCA.env$b_debug)
				{
					cat("Cutting Action:\t\t[ ", n_nodeid_cut_temp, " ] -> [ ", o_updatedindex_vector[1], ", ", o_updatedindex_vector[2], " ] >>>>>>> SUCCESS!\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
					cat("Left RowIDs:\t\t[", as.numeric(rSCA.env$o_output_tree[[o_updatedindex_vector[1]]]$rowids_matrix),"]\r\n", file = rSCA.env$s_logfilepath, sep = " ", append = TRUE)
					cat("Right RowIDs:\t\t[", as.numeric(rSCA.env$o_output_tree[[o_updatedindex_vector[2]]]$rowids_matrix), "]\r\n", file = rSCA.env$s_logfilepath, sep = " ", append = TRUE)
					cat("------------------------------------------\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
				}
			}
			else
			{
				#: is a leaf, set it as a leaf and add into merge stack
				o_updatedindex_vector_temp = f_setnode_asleaf(n_nodeid_cut_temp)
				rSCA.env$o_nodeid_stack_merge[rSCA.env$n_nodeid_statck_merge_cursor] <- n_nodeid_cut_temp
				rSCA.env$n_nodeid_statck_merge_cursor <- rSCA.env$n_nodeid_statck_merge_cursor + 1

				#: delete this node from cut stack
				rSCA.env$n_nodeid_statck_cut_cursor <- rSCA.env$n_nodeid_statck_cut_cursor - 1

				#cat("Cutting node_id: ", n_nodeid_cut_temp, " ->> FAILED!");cat("\n")
			}
		}

		rSCA.env$o_nodeid_stack_cut <- c()
		rSCA.env$n_nodeid_statck_cut_cursor <- 1
		if (rSCA.env$b_debug)
		{
			#cat("Current Merging Stack:\t\t[", rSCA.env$o_nodeid_stack_merge, "]\r\n", file = rSCA.env$s_logfilepath, sep = " ", append = TRUE)
			#cat("..............Cutting Finished................\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		}

		#: clear merge flag
		rSCA.env$n_flag_merge <- 0

		#: do merge if the merge stack is not empty
		if(rSCA.env$n_nodeid_statck_merge_cursor > 1)
		{
			#: initiate the loop
			rSCA.env$n_flag_merge <- 1
			while(rSCA.env$n_flag_merge == 1)
			{
				rSCA.env$n_flag_merge <- 0
				
				#: bottom index of the merge stack -> used for the no merging case in one for loop
				n_bot_index_merge_stack = 1				
				while(rSCA.env$n_nodeid_statck_merge_cursor > n_bot_index_merge_stack)
				{
					imerge_a = rSCA.env$n_nodeid_statck_merge_cursor - 1

					#: read a node from stack
					n_nodeid_merge_a = rSCA.env$o_nodeid_stack_merge[imerge_a]

					#: try to merge with other nodes (including the newly merged nodes)
					imerge_b = imerge_a - 1
					while(imerge_b >= n_bot_index_merge_stack)
					{
						n_nodeid_merge_b = rSCA.env$o_nodeid_stack_merge[imerge_b]

						#: check if these two nodes have tried to be merged
						if (f_hash_read(n_nodeid_merge_a, n_nodeid_merge_b) > 0)
						{
							#: if yes, continue other nodes
							#cat("Merging node_ids:[", n_nodeid_merge_a, ",", n_nodeid_merge_b, "] ->> NO NEED!");cat("\n")
							imerge_b = imerge_b - 1
							next
						}

						#: if the total number of each sample (a, b) is lower than (rSCA.env$n_sample_y_cols + 1)
						#: then the wilks value is 0, but now they can not be merged!!!
						if (nrow(rSCA.env$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix) <= (rSCA.env$n_sample_y_cols + 1) || nrow(rSCA.env$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix) <= (rSCA.env$n_sample_y_cols + 1))
						{
							#:update the comparing relationship between these two nodes
							o_temp_vector_r_a = rSCA.env$o_output_tree[[n_nodeid_merge_a]]$compared
							o_temp_vector_r_a[length(o_temp_vector_r_a) + 1] = n_nodeid_merge_b
							rSCA.env$o_output_tree[[n_nodeid_merge_a]]$compared <- o_temp_vector_r_a

							o_temp_vector_r_b = rSCA.env$o_output_tree[[n_nodeid_merge_b]]$compared
							o_temp_vector_r_b[length(o_temp_vector_r_b) + 1] = n_nodeid_merge_a
							rSCA.env$o_output_tree[[n_nodeid_merge_b]]$compared <- o_temp_vector_r_b

							#cat("Merging node_ids:[", n_nodeid_merge_a, ",", n_nodeid_merge_b, "] ->> FAILED!");cat("\n")
							imerge_b = imerge_b - 1
							next
						}

						#: construct the top and bot matrix and calculate wilks value
						
						o_top_matrix_temp = matrix(0, nrow(rSCA.env$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix), rSCA.env$n_sample_y_cols)
						o_top_matrix_temp[1:nrow(o_top_matrix_temp), ] = data.matrix(rSCA.env$o_sample_data_y[c(rSCA.env$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix), ])

						o_bot_matrix_temp = matrix(0, nrow(rSCA.env$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix), rSCA.env$n_sample_y_cols)
						o_bot_matrix_temp[1:nrow(o_bot_matrix_temp), ] = data.matrix(rSCA.env$o_sample_data_y[c(rSCA.env$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix), ])

						n_wilks_value = f_wilks_statistic(o_top_matrix_temp, o_bot_matrix_temp)

						#: 2> contruct min wilks list
						o_temp_min_wilks_list = list(min_wilks_value=n_wilks_value, col_id=0, x_value=0, left_rowids=rSCA.env$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix, right_rowids=rSCA.env$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix)

						#: 3> check if they can be merged
						n_merge_flag = f_cal_chk_f(o_temp_min_wilks_list)

						if (n_merge_flag == 0)
						{
							#: can be merged
							#: 1> do merge
							o_new_merged_vector = f_mergenodes(n_nodeid_merge_a, n_nodeid_merge_b, o_temp_min_wilks_list)
							
							#: 2> store the new node into merge stack at [imerge_b]
							rSCA.env$o_nodeid_stack_merge[imerge_b] <- o_new_merged_vector[1]

							#: delete [imerge_a] from merge stack
							rSCA.env$n_nodeid_statck_merge_cursor <- imerge_a

							#: 3> update merge times variable
							rSCA.env$n_flag_merge <- 1
							rSCA.env$n_merge_times <- rSCA.env$n_merge_times + 1
							if (rSCA.env$b_debug)
							{
								cat("Merging Action:\t\t[ ", n_nodeid_merge_a, ",", n_nodeid_merge_b, " ] -> [ ", o_new_merged_vector[1] , " ] >>>>>>> SUCCESS!\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
								cat("Sub RowIDs (1):\t\t[", as.numeric(rSCA.env$o_output_tree[[n_nodeid_merge_a]]$rowids_matrix),"]\r\n", file = rSCA.env$s_logfilepath, sep = " ", append = TRUE)
								cat("Sub RowIDs (2):\t\t[", as.numeric(rSCA.env$o_output_tree[[n_nodeid_merge_b]]$rowids_matrix), "]\r\n", file = rSCA.env$s_logfilepath, sep = " ", append = TRUE)
								cat("------------------------------------------\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
							}

							#: 4> break the FOR loop
							break
						}
						else
						{
							#:update the comparing relationship between these two nodes
							f_hash_store(n_nodeid_merge_a, n_nodeid_merge_b, 1)

							imerge_b = imerge_b - 1

							#cat("Merging node_ids: [", n_nodeid_merge_a, ",", n_nodeid_merge_b, "] ->> FAILED!");cat("\n")
						}
					}
					if (rSCA.env$n_flag_merge == 1)
					{
						#: there is merging action occurring, then restart the whole loop
						break
					}
					else
					{
						#: can not be merged with other nodes
						#: exchange imerge_a with the node pointed by n_bot_index_merge_stack
						o_temp_merge = rSCA.env$o_nodeid_stack_merge[imerge_a]
						rSCA.env$o_nodeid_stack_merge[imerge_a] <- rSCA.env$o_nodeid_stack_merge[n_bot_index_merge_stack]
						rSCA.env$o_nodeid_stack_merge[n_bot_index_merge_stack] <- o_temp_merge
						n_bot_index_merge_stack = n_bot_index_merge_stack + 1
					}
				}
			}
			
			#: copy all node in the merge stack into the cut stack, continue to cut
			rSCA.env$o_nodeid_stack_cut <- rSCA.env$o_nodeid_stack_merge[1:(rSCA.env$n_nodeid_statck_merge_cursor-1)]
			rSCA.env$n_nodeid_statck_cut_cursor <- rSCA.env$n_nodeid_statck_merge_cursor
			if (rSCA.env$b_debug)
			{
				#cat("Current Cutting Stack:\t\t[", rSCA.env$o_nodeid_stack_cut, "]\r\n", file = rSCA.env$s_logfilepath, sep = " ", append = TRUE)
				#cat("..............Merging Finished................\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
			}

			rSCA.env$o_nodeid_stack_merge <- c()
			rSCA.env$n_nodeid_statck_merge_cursor <- 1
		}
	}
}

# ---------------------------------------------------------------
# Initialization functions
# ---------------------------------------------------------------
f_init = function()
{
	if (rSCA.env$b_debug)
		cat("Initializing...\t\t", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
	
	#: print info to screen
	cat("Initializing...\t\t")

	#: matrix to store the sorted results
	rSCA.env$o_sorted_matrix <- matrix(0, rSCA.env$n_sample_size, rSCA.env$n_sample_x_cols)
	rSCA.env$o_sorted_temp_matrix <- matrix(0, rSCA.env$n_sample_size, rSCA.env$n_sample_x_cols)

	#: do sorting
	for (col in 1:rSCA.env$n_sample_x_cols)
	{
		rSCA.env$o_sorted_temp_matrix[,col] <- order(rSCA.env$o_sample_data_x[,col])
	}
	for (col in 1:rSCA.env$n_sample_x_cols)
	{
		for (row in 1:rSCA.env$n_sample_size)
		{
			rSCA.env$o_sorted_matrix[rSCA.env$o_sorted_temp_matrix[row, col], col] <- row
		}
	}

	#: list for output tree
	#: list(id=1, col_index=1, value=1, left=1, right=1, rowids_matrix=matrix())
	rSCA.env$o_output_tree <- list()

	#: initiate output tree
	o_init_tree_list = list(id=1, col_index=0, value=0, left=0, right=0, rowids_matrix=matrix(1:rSCA.env$n_sample_size, , 1))
	rSCA.env$o_output_tree[[1]] <- o_init_tree_list

	#: stack to store the unprocessed node ids in the output tree
	rSCA.env$o_nodeid_stack_cut <- c()
	rSCA.env$n_nodeid_statck_cut_cursor <- 1
	rSCA.env$o_nodeid_stack_merge <- c()
	rSCA.env$n_nodeid_statck_merge_cursor <- 1

	#: cutting and merging flags for loop, 1:do loop, 0:no need
	rSCA.env$n_flag_cut <- 0
	rSCA.env$n_flag_merge <- 0

	#: some statistical infos for the results
	rSCA.env$n_cut_times <- 0
	rSCA.env$n_merge_times <- 0
	rSCA.env$n_leafnodes_count <- 0

	#: hash table, hash function = (a + b) % rSCA.env$n_sample_size
	#: structure: [value, a, b, pointer]
	rSCA.env$o_hashtable_matrix <- matrix(0, rSCA.env$n_sample_size, 4)

	#: hash table list: (value, a, b, pointer)
	rSCA.env$o_hashtable_list <- list()
	rSCA.env$n_hashtable_list_index <- 0

	if (rSCA.env$b_debug)
		cat("SUCCESS!\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		
	#: print info to screen
	cat("SUCCESS!\r\n")
}

# ---------------------------------------------------------------
# Main functions
# ---------------------------------------------------------------
f_main = function()
{
	if (rSCA.env$b_debug)
	{
		cat("Training:\t\tIN PROGRESS!\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
	}
	
	#: print info to screen
	cat("Training...\t\t")
	
	#: cut from the root node
	f_processnode(1)

	#: results matrix structure -> matrix(id, col_id, x_value, left_id, right_id)
	o_results_matrix = matrix(0, length(rSCA.env$o_output_tree), 5)
	
	#: if mapvalue set as "mean" ==> calculate the mean of all samples
	#: else (i.e. set as "interval") ==> find the max and min of all samples to contruct a interval
	n_mapfile_cols = rSCA.env$n_sample_y_cols
	if (rSCA.env$n_mapvalue == "mean") n_mapfile_cols = rSCA.env$n_sample_y_cols
	if (rSCA.env$n_mapvalue == "max") n_mapfile_cols = rSCA.env$n_sample_y_cols
	if (rSCA.env$n_mapvalue == "min") n_mapfile_cols = rSCA.env$n_sample_y_cols
	if (rSCA.env$n_mapvalue == "median") n_mapfile_cols = rSCA.env$n_sample_y_cols
	if (rSCA.env$n_mapvalue == "interval") n_mapfile_cols = rSCA.env$n_sample_y_cols * 2
	if (rSCA.env$n_mapvalue == "radius") n_mapfile_cols = rSCA.env$n_sample_y_cols * 2
	if (rSCA.env$n_mapvalue == "variation") n_mapfile_cols = rSCA.env$n_sample_y_cols * 2
	if (rSCA.env$n_mapvalue == "random") n_mapfile_cols = rSCA.env$n_sample_y_cols * 2
	
	o_y_results_matrix = matrix(0, length(rSCA.env$o_output_tree), n_mapfile_cols)

	for(itree in 1:length(rSCA.env$o_output_tree))
	{
		o_results_matrix[itree,1] = as.numeric(rSCA.env$o_output_tree[[itree]]$id)
		o_results_matrix[itree,2] = as.numeric(rSCA.env$o_output_tree[[itree]]$col_index)
		o_results_matrix[itree,3] = as.numeric(rSCA.env$o_output_tree[[itree]]$value)
		o_results_matrix[itree,4] = as.numeric(rSCA.env$o_output_tree[[itree]]$left)
		o_results_matrix[itree,5] = as.numeric(rSCA.env$o_output_tree[[itree]]$right)

		#: if is leaf, calculate y accordingly: mean or interval
		if (rSCA.env$o_output_tree[[itree]]$left == -1 && rSCA.env$o_output_tree[[itree]]$right == -1)
		{
			rSCA.env$n_leafnodes_count <- rSCA.env$n_leafnodes_count + 1
			n_y_size = nrow(rSCA.env$o_output_tree[[itree]]$rowids_matrix)
			o_y_matrix = matrix(0, n_y_size, rSCA.env$n_sample_y_cols)

			for (iy in 1:n_y_size)
			{
				n_row_id_y = rSCA.env$o_output_tree[[itree]]$rowids_matrix[iy, 1]
				#: NOTE: need to convert to numeric
				o_y_matrix[iy, ] = as.numeric(rSCA.env$o_sample_data_y[n_row_id_y, ])
			}
			
			#: processing mapvalue
			if (rSCA.env$n_mapvalue == "mean")
			{
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), ] = apply(o_y_matrix, MARGIN=c(2), mean)
			}
			if (rSCA.env$n_mapvalue == "max")
			{
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), ] = apply(o_y_matrix, MARGIN=c(2), max)
			}
			if (rSCA.env$n_mapvalue == "min")
			{
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), ] = apply(o_y_matrix, MARGIN=c(2), min)
			}
			if (rSCA.env$n_mapvalue == "median")
			{
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), ] = apply(o_y_matrix, MARGIN=c(2), median)
			}
			if (rSCA.env$n_mapvalue == "interval")
			{
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), 1:rSCA.env$n_sample_y_cols] = apply(o_y_matrix, MARGIN=c(2), min)
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), (rSCA.env$n_sample_y_cols + 1):n_mapfile_cols] = apply(o_y_matrix, MARGIN=c(2), max)
			}
			if (rSCA.env$n_mapvalue == "radius")
			{
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), 1:rSCA.env$n_sample_y_cols] = apply(o_y_matrix, MARGIN=c(2), mean)
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), (rSCA.env$n_sample_y_cols + 1):n_mapfile_cols] = apply(o_y_matrix, MARGIN=c(2), f_cal_radius)
			}
			if (rSCA.env$n_mapvalue == "variation")
			{
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), 1:rSCA.env$n_sample_y_cols] = apply(o_y_matrix, MARGIN=c(2), mean)
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), (rSCA.env$n_sample_y_cols + 1):n_mapfile_cols] = apply(o_y_matrix, MARGIN=c(2), sd)
			}
			if (rSCA.env$n_mapvalue == "random")
			{
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), 1:rSCA.env$n_sample_y_cols] = apply(o_y_matrix, MARGIN=c(2), min)
				o_y_results_matrix[as.numeric(rSCA.env$o_output_tree[[itree]]$id), (rSCA.env$n_sample_y_cols + 1):n_mapfile_cols] = apply(o_y_matrix, MARGIN=c(2), max)
			}
		}
	}
	
	#: add matrix column names
	colnames(o_results_matrix) = c("NID", "xCol", "xVal", "lfNID", "rtNID")
	
	if (rSCA.env$n_mapvalue == "mean" || rSCA.env$n_mapvalue == "max" || rSCA.env$n_mapvalue == "min" || rSCA.env$n_mapvalue == "median")
		colnames(o_y_results_matrix) = colnames(rSCA.env$o_sample_data_y)
	if (rSCA.env$n_mapvalue == "interval" || rSCA.env$n_mapvalue == "radius" || rSCA.env$n_mapvalue == "variation")
		colnames(o_y_results_matrix) = c(colnames(rSCA.env$o_sample_data_y), colnames(rSCA.env$o_sample_data_y))
	
	#: store the tree and map files
	write.table(o_results_matrix, file = rSCA.env$s_tree_filepath, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
	write.table(o_y_results_matrix, file = rSCA.env$s_map_filepath, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
	
	
	if (rSCA.env$b_debug)
	{
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("Training:\t\tCOMPLETE!\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("Total Nodes:\t\t[ ", length(rSCA.env$o_output_tree), " ]\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("Leaf Nodes:\t\t[ ", rSCA.env$n_leafnodes_count, " ]\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("Cutting Actions:\t[ ", rSCA.env$n_cut_times, " ]\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("Merging Actions:\t[ ", rSCA.env$n_merge_times, " ]\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("Tree File:\t\t[ ", rSCA.env$s_tree_file, " ]\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("Map File:\t\t[ ", rSCA.env$s_map_file, " ]\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
	}
	
	#: print info to sreen
	cat("SUCCESS!\r\n")
	cat("Total Nodes:\t\t", length(rSCA.env$o_output_tree), "\r\n", sep = "")
	cat("Leaf Nodes:\t\t", rSCA.env$n_leafnodes_count, "\r\n", sep = "")
	cat("Cutting Actions:\t", rSCA.env$n_cut_times, "\r\n", sep = "")
	cat("Merging Actions:\t", rSCA.env$n_merge_times, "\r\n", sep = "")
	cat("Tree File:\t\t", rSCA.env$s_tree_file, "\r\n", sep = "")
	cat("Map File:\t\t", rSCA.env$s_map_file, "\r\n", sep = "")
	if (rSCA.env$b_debug)
	{
		cat("Log File:\t\t", rSCA.env$s_logfilepath, "\r\n", sep = "")
	}
}

# ---------------------------------------------------------------
# Interface function
# ---------------------------------------------------------------
do_cluster = function()
{
	#: store the start time
	time_stat <- proc.time()

	#: initialize
	f_init()

	#: do main function
	f_main()

	# : calculate the total time used
	time_end <- (proc.time() - time_stat)[[3]]
	Hours <- time_end %/% (60*60)
	Minutes <- (time_end %% 3600) %/% 60
	Seconds <- time_end %% 60
	time_used <- paste(Hours, " h ", Minutes, " m ", Seconds, " s.", sep="")
	if (rSCA.env$b_debug)
	{
		cat("Time Used:\t\t[ ", time_used, " ]\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
		cat("==========================================\r\n", file = rSCA.env$s_logfilepath, sep = "", append = TRUE)
	}
	
	#: print info to screen
	cat("Time Used:\t\t", time_used, "\r\n", sep = "")
}

