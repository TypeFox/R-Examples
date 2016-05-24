setRFoptions <- function(ntree_cat = NULL, ntree_cont = NULL,
	nodesize_cat = NULL, nodesize_cont = NULL,
	maxnodes_cat = NULL, maxnodes_cont = NULL){
	# Records these global options for use by CALIBERrfimpute
	for (opname in c('ntree_cat', 'ntree_cont',
		'nodesize_cat', 'nodesize_cont',
		'maxnodes_cat', 'maxnodes_cont')){
		fullopname <- paste('CALIBERrfimpute', opname, sep='_')
		if (!is.null(get(opname))){
			message(paste(c('Setting option ', fullopname, ' = ',
				get(opname)), collapse=''))
			eval(parse(text=paste('options(', fullopname, '=', get(opname), ')')))
		}
	}
}

