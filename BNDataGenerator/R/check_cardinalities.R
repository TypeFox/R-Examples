# Written by Jae-seong Yoo 20141217

check_cardinality = function (arcs_mat, nodename=NULL, cardinality=NULL)
{
	# Check DAG
	check_dag_arcs = as.matrix(arcs_mat)
	if (is_DAG(check_dag_arcs) == FALSE) {
		stop("arcs_mat must a DAG")
	}
	
	
	# Node 개수
	num_of_nodes = dim(arcs_mat)[1];

	
	# nodename가 NULL이면 임의로 node 이름을 부여한다.
	if (is.null(nodename)) {
		nodename = big_letters(num_of_nodes)
	}

	
	# Cardinality가 NULL이면 모두 2로 설정한다.
	# Cardinality는 모두 2보다 커야 한다.
	if (is.null(cardinality)) {
		cardinality = rep(2, num_of_nodes)
	} else if (sum(cardinality < 2) > 0) {
		stop("All cardinality must be at least 2.")
	} else if (num_of_nodes != length(cardinality)) {
		stop("Wrong length of cardinality")
	}

	
	# 각 Node의 Parent Node 개수
	num_of_parent_nodes = apply(arcs_mat, 2, sum);

	list_parent_nodes = list();
	for(m in 1:num_of_nodes) {
		if (length(which(arcs_mat[,m]==1)) == 0) {
			list_parent_nodes[[m]] = NULL;
		} else {
			list_parent_nodes[[m]] = which(arcs_mat[,m]==1);
		}
	}


	# Root node의 개수
	num_of_root_nodes = sum(num_of_parent_nodes == 0);


	# 지정해야할 조건부 확률 개수
	num_of_probs = NULL;
	for (k in 1:num_of_nodes) {
		num_of_probs[k] = (cardinality[k]-1) * prod(cardinality[list_parent_nodes[[k]]])
	}
	
	
	text_of_probs = list();
	for(i in 1:length(num_of_parent_nodes))
	{
		temp_text = NULL;
		present_cardinality = as.matrix(toss_value(1, cardinality[i]));
		
		if (num_of_parent_nodes[i] == 0) {		### it is root node
			for (j in 1:(cardinality[i]-1))
			{
				temp_text = c(	temp_text,
										paste(	"P(",
													nodename[i], " = ", present_cardinality[j,1],
													")",
													sep=""
												)
									)
			}
		} else {		### it is not a root node
			temp_list_of_pn = as.numeric(list_parent_nodes[[i]]);
			
			for (j in 1:(cardinality[i]-1))
			{
				temp_cases = list();
				cases = NULL;
				
				for (k in 1:length(temp_list_of_pn))
				{
					temp_cases[[k]] = toss_value(1, cardinality[temp_list_of_pn[k]])
					if (is.null(cases)) {
						cases = temp_cases[[k]]
						names(cases) = 1;
					} else {
						cases = merge(cases, temp_cases[[k]])
						names(cases) = c(1:dim(cases)[2])
					}
				}
				cases = as.matrix(cases)
				
				
				for (k in 1:dim(cases)[1])
				{
					temp_text_conditional = NULL;
					
					for (m in 1:dim(cases)[2])
					{
						case_value = paste(	nodename[temp_list_of_pn[m]],
														" = ", 
														cases[k, m],
														sep=""
													)
														
						if (m == 1) {
							temp_text_conditional = case_value;
						} else {
							temp_text_conditional = paste( temp_text_conditional, paste(", ", case_value), sep="" );
						}
					}
					
					temp_text = c(	temp_text,
											paste(	"P(",
														nodename[i], " = ", present_cardinality[j,1],
														"|",
														temp_text_conditional,
														")",
														sep=""
													)
										)
				}
			}
		}

		text_of_probs[[i]] = temp_text;
	}

	
	res = list(	cardinality = cardinality,
					nodename = nodename,
					num_of_root_nodes = num_of_root_nodes,
					num_of_probs = num_of_probs,
					num_of_parent_nodes = num_of_parent_nodes,
					list_parent_nodes = list_parent_nodes,
					list_of_probs = text_of_probs
				);

	return(res);
}