check.mass.balance <-
function(equation, chemical_table, id_col='chebi', formula_col='formula', direction_type=c(' <=> ', ' => ', ' <\\?> ')) {
	# Check parameters
	if(formula_col %in% colnames(chemical_table) ==F) {
		stop('Please provide proper column names\n')
	}

	# Split equation
	nonAtom = 'R'

	directionalities = character(length(equation))
	for(i in direction_type) {
  	ind_direction = grepl(i, equation)
  	directionalities[ind_direction] = i
	}

	# Check direction_type
	if(length(grep('^$', directionalities))>0) {
        print(grep('^$', directionalities, value=T))
		stop('Your equation has undefined direction symbol. Please provide proper direction symbol')
	}

	result = logical(length(equation))
	for(i in 1:length(equation)) {
		substrates = unlist(strsplit(equation[i], directionalities[i]))[1]
		substrates = unlist(strsplit(substrates, ' \\+ '))
    substrates = gsub('\\(n\\+1\\) ', '2 ', substrates)
    substrates = gsub('([0-9])n ', '\\1 ', substrates)
    substrates = gsub('n ', '1 ', substrates)
    
		substrates = gsub('\\(.+\\)', '', substrates)
    if(grepl('n-1', equation[i])) {
      substrates = gsub('n-1 ', '2 ', substrates)
      substrates = gsub('n-2 ', '1 ', substrates)
    }

		pattern = '([0-9]+) (.+)'

		sub_coefficient = rep(1, length(substrates))
		ind_sub = grep(pattern, substrates)
		if(length(ind_sub) > 0) {
			sub_coefficient[ind_sub] = sub(pattern, '\\1', substrates[ind_sub])	
			sub_coefficient = as.numeric(sub_coefficient)
		}
	
		substrates = sub(pattern, '\\2', substrates)
		substrates = rep(substrates, sub_coefficient)

		products = unlist(strsplit(equation[i], directionalities[i]))[2]
		products = unlist(strsplit(products, ' \\+ '))
    products = gsub('\\(n\\+1\\) ', '2 ', products)
    products = gsub('([0-9])n ', '\\1 ', products)
    products = gsub('n ', '1 ', products)
    
		products = gsub('\\(.+\\)', '', products)
    if(grepl('n-1', equation[i])) {
      products = gsub('n-1 ', '2 ', products)
      products = gsub('n-2 ', '1 ', products)
    }

		pro_coefficient = rep(1, length(products))
    ind_pro = grep(pattern, products)
    if(length(ind_pro) > 0) {
      pro_coefficient[ind_pro] = sub(pattern, '\\1', products[ind_pro])
    }
	
		products = sub(pattern, '\\2', products)
		products = rep(products, pro_coefficient)

		# Convert participants into chemical formula
		substrates_formula = lapply(substrates, .id2formula, chemical_table, id_col, formula_col)
		substrates_formula = unlist(substrates_formula)
    substrates_formula = gsub('\\)///\\(','', substrates_formula)
		substrates_formula = gsub(' ','', substrates_formula)
		substrates_formula = gsub('^\\(','', substrates_formula)
		substrates_formula = gsub('\\)$','', substrates_formula)
		products_formula = lapply(products, .id2formula, chemical_table, id_col, formula_col)
		products_formula = unlist(products_formula)
		products_formula = gsub('\\)///\\(','', products_formula)
		products_formula = gsub(' ','', products_formula)
		products_formula = gsub('^\\(','', products_formula)
		products_formula = gsub('\\)$','', products_formula)

		sub_formula_matrix = .formula2matrix(substrates_formula)
		pro_formula_matrix = .formula2matrix(products_formula)

		if(nonAtom %in% c(names(sub_formula_matrix), names(products_formula))) {
			result[i] = 'Not available'
		} else if(length(which(is.na(c(substrates_formula,products_formula)))) > 0) {
		  result[i] = 'Not available'
		} else if(length(products) != length(products_formula)) {
		  result[i] = 'Not available'
		} else if(length(substrates) != length(substrates_formula)) {
		  result[i] = 'Not available'
		} else {
			tmp_result = identical(sub_formula_matrix, pro_formula_matrix)
			result[i] = tmp_result
		}
	}
	return(result)
}
