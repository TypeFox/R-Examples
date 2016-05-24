
subtable <- function(data, cols, freqvar = "Freq", keep.zero=FALSE, allfactor =FALSE, return.type = class(data)){
	data=as.data.frame(data)
	
	
	if(!is.null(freqvar)){
		if(! (freqvar %in% names(data)) ){
			freqvar <- NULL
		} 
	}
	
	
	names(data)[which(names(data)==freqvar)] <- "Freq"
	
	if(allfactor){
		int = which(sapply(data[,cols], function(v) is.integer(v)|is.numeric(v)))
		for( i in int ){
			data[,cols[i]] = as.factor(data[,cols[i]])
		}
	}

	if(!keep.zero){
		if("Freq" %in% names(data)){
			fid = which(names(data)=="Freq")
			#ss = count2(data, cols, weights = data[,fid])
			ss = count(data, cols, wt_var = fid)
		}else{
			#ss = count2(data, cols)
			ss = count(data, cols)
		}
		names(ss)[names(ss)=="freq"] <- "Freq"
		
		ss <- ss[ss$Freq > 0,]
	}else{
		if("Freq" %in% names(data)){
			fid = which(names(data)=="Freq")
			ss = as.data.frame(xtabs(Freq~.,data=data[,c(cols,fid)]))
		}else{
			ss = as.data.frame(table(data[,c(cols)]))
		}
	}
#	if(allfactor){
#			ss = data.frame(sapply(ss[,-ncol(ss)],factor),ss$Freq)
#			names(ss)[ncol(ss)] = "Freq"
#	}
	if("table" %in% return.type )	ss <- xtabs(Freq~.,data=ss)
		
	return(ss)
}
	
# count2 = function (df, vars = NULL, weights = NULL) 
# {
    # if (is.vector(df)) {
        # df <- data.frame(x = df)
    # }
    # if (!is.null(vars)) {
        # vars <- as.quoted(vars)
        # df <- quickdf(eval.quoted(vars, df))
    # }
# if(!is.null(weights)){
    # id <- plyr:::ninteraction(df, drop = TRUE)
    # u_id <- !duplicated(id)
    # labels <- df[u_id, , drop = FALSE]
    # labels <- labels[order(id[u_id]), , drop = FALSE]
    # Freq <- xtabs(weights~id)
	# class(Freq) = "vector"
    # unrowname(data.frame(labels, Freq))
# }else{
	# id <- plyr:::ninteraction(df, drop = TRUE)
    # u_id <- !duplicated(id)
    # labels <- df[u_id, , drop = FALSE]
    # labels <- labels[order(id[u_id]), , drop = FALSE]
    # Freq <- tabulate(id, attr(id, "n"))
    # unrowname(data.frame(labels, Freq))
# }
# }

subtable.table <- function(x,cols){
	x <- as.table(x)
	x2 <- apply(x,cols,sum)
	dim(x2) <- dim(x)[cols]
	return(as.table(x2))
}

# subtable2 <- function(x, cols, freqvar = NULL, keep.zero=FALSE, allfactor =FALSE, return.type = class(x)){
	# require(data.table)
	# tab <- is.table(x)
	# x <- as.data.table(x)
	
	# if(is.null(freqvar) & tab){
		# freqvar <- "Freq"
		# setnames(x,"N","Freq")
	# }
	# if(is.null(freqvar) & "Freq" %in% names(x)){
		# freqvar <- "Freq"
	# }
	
	# if(!is.null(freqvar)){
		# setnames(x,freqvar,"Freq")
	# }
	
	
	# nv <- names(x)[cols]
	# setkeyv(x,nv)
	
	
	
	# if(is.null(freqvar)){
		# x <- x[,nrow(.SD),by=key(x)] 
	# }else{
		# if(!keep.zero){
			# x <- x[Freq>0]
		# }
		# x <- x[,sum(Freq),by=key(x)]
	# }
	
	# nmz <- names(x)
	# nmz[length(nmz)] <- "Freq"
	# setnames(x,nmz)
	
	# if(!keep.zero){
		# x <- x[Freq>0]
	# }
	
	# if(allfactor){
		# inum <- function(v) is.integer(v)|is.numeric(v)
		# int <- which(!sapply(x,inum))
		# for( i in int ){
			# vn <- names(x)[cols[i]]
			# if(vn != "Freq"){
				# x[,cols[i]] <- as.factor(x[,cols[i]])
			# }
		# }
	# }
	
	
	# if("table" %in% return.type ){
		# x <- xtabs(Freq~.,data=x)
		# return(x)
	# }	
	
	# if(!is.null(freqvar)){
		# nmz[length(nmz)] <- freqvar
		# setnames(x,nmz)
	# }
	# return(x)
# }

