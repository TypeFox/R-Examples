# Written by Jae-seong Yoo 20141101

BN_Data_Generator = function (arcs_mat, Probs, n, nodename=NULL, cardinality=NULL)
{
	# Check Sample Size
	if (n <= 0) {
		stop("Sample size 'n' must be greater than 0.")
	}

	# 20141209: sample size가 1000개보다 적으면 데이터가 올바르게 생성되지 않는 버그가 있다.
	# 이를 보완하기 위한 부분.
	if (n < 10000) {
		temp_n = 1000;
	} else {
		temp_n = n;
	}
	#####

	# Node 개수
	num_of_nodes = dim(arcs_mat)[1];
	
	# 결과는 여기에 저장이 된다.
	result_mat = matrix(0, temp_n, num_of_nodes);
	dimnames(result_mat)[[2]] = nodename;
	# result_mat

	# Check Input Probs & cardinality
	checker = check_cardinality(arcs_mat = arcs_mat, nodename = nodename, cardinality = cardinality)
	cardinality = checker$cardinality;
	nodename = checker$nodename;
	list_parent_nodes = checker$list_parent_nodes;
	num_of_probs = checker$num_of_probs;
	num_of_parent_nodes = checker$num_of_parent_nodes;
	num_of_root_nodes = checker$num_of_root_nodes;

	
	# 지정해야할 조건부 확률 개수만큼 input이 맞는지 확인. 만일 false이면 프로그램 종료
	input_prob_len = length(Probs);
	for (i in 1:input_prob_len) {
		if	(	as.numeric(length(Probs[[i]]))
								!=
				as.numeric(num_of_probs[i])
			)
		{
			stop("Input Probs != num_of_probs!");
		}
	}
	
	# Root Node Initialization
	for(i in 1:num_of_root_nodes) {
		p = Probs[[i]];
		mat_values = merge("Value", c(1:cardinality[i]))
		mat_values = paste(mat_values[,1], mat_values[,2], sep="")

		result_mat[,i] = sample(
											mat_values, temp_n,
											prob=c(p, 1-sum(p)), replace=T
											);
	}


	# Generator
	init = num_of_root_nodes + 1;

	mat = NULL
	for (i in init:num_of_nodes) {
		p = Probs[[i]]
		temp_list_of_pn = as.numeric(list_parent_nodes[[i]]);
		num_of_c_cases = prod(cardinality[temp_list_of_pn])

		temp_cases = list();
		cases = NULL;
		for (j in 1:length(temp_list_of_pn))
		{
			temp_cases[[j]] = toss_value(1, cardinality[temp_list_of_pn[j]])
			if (is.null(cases)) {
				cases = temp_cases[[j]]
				names(cases) = 1;
			} else {
				cases = merge(cases, temp_cases[[j]])
				names(cases) = c(1:dim(cases)[2])
			}
		}
		cases = as.matrix(cases)
		mat_values = merge("Value", c(1:cardinality[i]))
		mat_values = sort(paste(mat_values[,1], mat_values[,2], sep=""))
		
		stack = 1
		for(j in 1:dim(cases)[1])
		{	
			mat = t(
						t(
							as.matrix(result_mat[,temp_list_of_pn])
						) ==
							cases[j,]
					);
			mat = (apply(mat, 1, sum) == dim(mat)[2])
			
			if (cardinality[i] == 2)
			{
				temp_p = p[j]
			} else {
				temp_p = p[stack:(stack + cardinality[i]-2)]
			}
			len = length(which(mat))
			
			result_mat[which(mat),i] = sample(
															mat_values, len,
															prob=c(temp_p, 1-sum(temp_p)), replace=T
														);
			
			stack = stack + (cardinality[i]-1)
		}
	}



	# 20141209: sample size가 1000개보다 적으면 데이터가 올바르게 생성되지 않는 버그가 있다.
	# 이를 보완하기 위한 부분.
	if (n < 1000) {
		result_mat = result_mat[sample(c(1:1000), size=n), ]
	}
	#####

	res = list(	data = data.frame(result_mat),
					nodename = nodename,
					num_of_nodes = num_of_nodes,
					num_of_parent_nodes = num_of_parent_nodes,
					list_parent_nodes = list_parent_nodes);
	return(res);
}