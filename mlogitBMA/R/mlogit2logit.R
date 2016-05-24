mlogit2logit <- function(f, data, choices=NULL, base.choice=1, varying = NULL, sep='.') {
	mnlspec <- mnl.spec(f, data=data, choices=choices, base.choice=base.choice, varying=varying, sep=sep)
	binlog <- do.mlogit2logit(data=data, spec=mnlspec)
	return(binlog)
}


do.mlogit2logit <- function(data, spec) {

	resp.star <- paste(spec$response,'.star', sep='')
	data <- cbind(data, 0)
	colnames(data)[ncol(data)] <- resp.star
	
	nmodes <- length(spec$choices)
	base.mode <- spec$choices[spec$base.choice]
	other.modes <- setdiff(spec$choices, base.mode)
	other.modes.ids <- setdiff(1:nmodes, spec$base.choice)
	
	# add variables of differences (from base alt.) for alternative specific variables
	delta.vars <- altspec.full <- list()
	if (length(spec$varying.names) > 0) {
		for (ivar in 1:length(spec$varying.names)) {
			this.name <- spec$varying.names[ivar]
			altspec.full[[this.name]] <- spec$full.var.names[,this.name]
			names(altspec.full[[this.name]]) <- spec$choices
			delta.vars[[this.name]] <- paste('d', 
				altspec.full[[this.name]][other.modes.ids], sep='_')
			names(delta.vars[[this.name]]) <- spec$choices[other.modes.ids]
			values <- data[,altspec.full[[this.name]][other.modes.ids]] - data[,altspec.full[[this.name]][base.mode]]
			data <- cbind(data, values)
			colnames(data)[(ncol(data)-length(delta.vars[[this.name]])+1):ncol(data)] <- delta.vars[[this.name]]
		}		
	}
	
	f <- paste(resp.star, '~ ')

	idx.base.data <- data[,spec$response]==base.mode
	base.data <- data[idx.base.data,]
	varying.coefs <- colnames(spec$variable.used)[!spec$same.coefs]
	lvarying.coefs <- length(varying.coefs)
	z.index <- z.names <-  c()
	zcols <- list()
	bin.data <- c()
	add.cols <- c()
	# for each alternative add a block to the data
	for (choice in other.modes) {
		#create a dataset
		idx <- data[,spec$response]==choice
		this.data <- data[idx,]
		this.data[,resp.star] <- 1
		mod.base.data <- base.data
		if (length(spec$varying.names) > 0) { # create W columns
			for (ivar in 1:length(spec$varying.names)) {
				this.data <- cbind(this.data, data[idx, delta.vars[[spec$varying.names[ivar]]][as.character(choice)]])
				colnames(this.data)[ncol(this.data)] <- spec$varying.names[ivar]
				mod.base.data <- cbind(mod.base.data, data[idx.base.data, 
									delta.vars[[spec$varying.names[ivar]]][as.character(choice)]])
				colnames(mod.base.data)[ncol(mod.base.data)] <- spec$varying.names[ivar]
			}
		}
		this.data <- rbind(mod.base.data, this.data)
		ladd.cols <- length(add.cols)
		if(ladd.cols > 0) {
			this.data <- cbind(this.data, matrix(0, nrow=nrow(this.data), ncol=ladd.cols))
			colnames(this.data)[(ncol(this.data)-ladd.cols+1):ncol(this.data)] <- add.cols
		}
		#create dummy z-variables
		intname <- paste('z_',choice, sep='') 
		this.data <- cbind(this.data, 1)
		colnames(this.data)[ncol(this.data)] <- intname
		zcols[[intname]] <- c(name=intname, choice=choice, intercept=TRUE)
		z.index <- c(z.index, ncol(this.data))
		z.names <- c(z.names, intname)
		# create interactions for variables whose coefficients differs in each alternative 
		zvars <- c()
		if (lvarying.coefs > 0) {
			zvars <- paste(intname, varying.coefs, sep='.')
			for (ivar in 1:length(zvars)) {
				values <- this.data[,intname] * this.data[,varying.coefs[ivar]]
				this.data <- cbind(this.data, values)
				colnames(this.data)[ncol(this.data)] <- zvars[ivar]
				zcols[[zvars[ivar]]] <- c(name=varying.coefs[ivar], choice=choice, intercept=FALSE)
			}
			sep <- ''
			if (!is.null(bin.data)) {
				bin.data <- cbind(bin.data, matrix(0, ncol=lvarying.coefs, nrow=nrow(bin.data)))
				colnames(bin.data)[(ncol(bin.data)-lvarying.coefs+1):ncol(bin.data)] <- zvars
				sep <- '+'
			} 
			f <- paste(f, paste(zvars, collapse='+'), sep=sep)
		}
		if (!is.null(bin.data)) {
			bin.data <- cbind(bin.data, 0)
			colnames(bin.data)[ncol(bin.data)] <- intname
			if (spec$intercepts[other.modes.ids[1]])
				f <- paste(f, intname, sep='+')
		}
		add.cols <- c(add.cols, intname, zvars)
		bin.data <- rbind(bin.data, this.data)
	}
	constant.coefs <- colnames(spec$variable.used)[spec$same.coefs]
	if (length(constant.coefs) > 0)
		f <- paste(f, paste(constant.coefs, collapse='+'), sep='+')
	if (!spec$intercepts[other.modes.ids[1]])
		f <- paste(f, '-1')
	return(list(data=bin.data, formula=formula(f), nobs=dim(data)[1], z.index=z.index, z.names=z.names,
				zcols=zcols, choices=spec$choices, 
				choice.main.intercept=other.modes.ids[1]))
}
