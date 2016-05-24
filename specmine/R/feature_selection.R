###############################################################
#####################FEATURE SELECTION#########################
###############################################################

# feature selection
# dataset: data and metadata structure
# column.class: metadata column class
# method: "rfe" (recursive feature elimination) or "filter" (feature selection using univariate filters)
# functions:  
feature_selection = function(dataset, column.class, method = "rfe", functions, validation = "cv", 
                             repeats = 5, number = 10, subsets = 2^(2:4)){
	if (method == "rfe"){
		result = recursive_feature_elimination(dataset$data, dataset$metadata[,column.class], functions, 
                                           validation, repeats, number, subsets)
	} 
  else if (method == "filter"){
		result = filter_feature_selection(dataset$data, dataset$metadata[,column.class], functions, 
                                      validation, repeats)
	}
	result
}


#Recursive Feature Elimination
#funcs list: lmFuncs, rfFuncs, treebagFuncs, ldaFuncs, nbFuncs, gamFuncs, lrFuncs
recursive_feature_elimination = function(datamat, samples.class, functions = caret::rfFuncs, method = "cv", 
                                         repeats = 5, number = 10, subsets = 2^(2:4)){
	samples.df = data.frame(t(datamat))
	ctrl <- caret::rfeControl(functions = functions,
                   method = method,
                   repeats = repeats,
                   number = number,
                   verbose = FALSE)
	rfe.result <- caret::rfe(samples.df, samples.class,
                 sizes = subsets,
                 rfeControl = ctrl)
	
	rfe.result
}


#Feature Selection Using Univariate Filters
#functions list: lmSBF, rfSBF, treebagSBF, ldaSBF and nbSBF.
filter_feature_selection = function(datamat, samples.class, functions = caret::rfSBF, method = "cv", 
                                    repeats = 5) { 
	samples.df = data.frame(t(datamat))
	filterCtrl = caret::sbfControl(functions = functions, method = method, repeats = repeats)
	filter.result = caret::sbf(samples.df, samples.class, sbfControl = filterCtrl)
	filter.result
}
