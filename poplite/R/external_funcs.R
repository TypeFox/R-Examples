#Went with an S3 method here and for select for the S4 Database class to go with the S3 generics in dplyr
filter_.Database <- function(.data, ...,.dots)
	  {
	    
	    use.expr <- all_dots(.dots, ..., all_named = TRUE)
	    
	    if (length(use.expr) != 1)
	    {
		stop("ERROR: Please supply a single statement which would result in a logical vector")
	    }
	    
	    .get.var.names <- function(x)
	    {
		if (length(x) == 1)
		{
		    return(as.character(x))
		}else{
		    which.calls = sapply(x, class)
		    if (any(which.calls %in% c('call', '(')))
		    {
			return(lapply(x[which.calls %in% c('call', '(')], .get.var.names))
		    }else{
			return(as.character(x[[2]]))
		    }
		    
		}
	    }
	    
	    needed.vars <- unlist(.get.var.names(use.expr[[1]]))
	    
	    #figure out which tables the requested variables are in
	    
	    parse.tables <- get.tables.from.vars(needed.vars)
	    
	    not.prov.tabs <- sapply(parse.tables, is.null)
	    
	    if (any(not.prov.tabs))
	    {
		col.to.tab <- stack(columns(.data))
		
		np.tabs <- lapply(needed.vars[not.prov.tabs], function(x)
				  {
					found.tabs <- as.character(col.to.tab$ind[col.to.tab$values %in% x])
				  })
		
		np.tabs <- append(np.tabs, parse.tables[not.prov.tabs == F])
		
		np.tab.len <- sapply(np.tabs, length)
		
		if(any(np.tab.len == 0))
		{
		    stop(paste("ERROR: Cannot find corresponding table for", paste(needed.vars[np.tab.len == 0],collapse=",")))
		}
		else if (any(np.tab.len > 1))
		{
		    #if all other columns match one specific table, and the multi-mapping col does too, then use that table
		    ##otherwise throw an error
		    
		    unique.tabs <- np.tabs[np.tab.len == 1]
		    
		    conc.tab.list <- lapply(np.tabs[np.tab.len > 1], function(x) x[x %in% unlist(unique.tabs)])
		    conc.list.len <- sapply(conc.tab.list, length)
		    
		    if (all(length(conc.list.len) == 1) && length(unique(unlist(unique.tabs))) == 1)
		    {
			#as everybody agrees on a given table
			np.tabs <- unique.tabs[1]
			
		    }else{
			stop("ERROR: Cannot uniquely map columns to table, try using: tableX.columnY")
		    }
		    
		}
		
		needed.tables <- unique(unlist(np.tabs))
		
	    }else{
		#names are provided for all columns
		needed.tables <- unique(unlist(parse.tables))
	    }
	    
	    my_db_tbl <- join(.data, needed.tables)
	    
	    #carry out the filter method of dplyr on the temporary or otherwise table
	    
	    cur.stat <- paste(deparse(use.expr[[1]]$expr), collapse=" ")
	    
	    for(i in needed.tables)
	    {
		cur.stat <- gsub(paste0(i, "."), "", cur.stat)
	    }
	   
	    #use.statement <- paste("filter(my_db_tbl,", cur.stat, ")")
	    #return(eval(parse(text=use.statement)))
	    return (filter_(my_db_tbl, cur.stat))
    }

select <- function(.data,..., .tables=NULL)
{   
    use.dots <- lazy_dots(...)
    
    if ((missing(.tables) || is.null(.tables) || is.na(.tables))==F && inherits(.data, "Database") == TRUE)
    {
	use.dots <- append(use.dots, list(.tables=.tables))
    }
    
    select_(.data, .dots = use.dots)
}

select_.Database <- function(.data, ..., .dots)
{
    if (missing(.dots) || is.null(.dots) || all(is.na(.dots)))
    {
      .dots <- list()
    }
  
    #check to see if .tables is part of .dots
    if('.tables' %in% names(.dots))
    {
	.tables <- .dots$.tables
	.dots <- .dots[-which(names(.dots) == ".tables")]
    }else{
	.tables <- NULL
    }
    
    use.expr <- all_dots(.dots, ..., all_named = TRUE)
    
    if (is.null(.tables) == F)
    {
	if (is.character(.tables) && length(.tables) >= 1 && all(.tables %in% tables(.data)))
	{
	   use.tables <- .tables
	}
	else
	{
	    stop("ERROR: .tables needs to be a vector of table names use 'tables(.data)' for a listing")
	}
	
	clean.cols <- sapply(use.expr, function(x) deparse(x$expr))
	
	if (length(clean.cols) > 1)
	{
	    stop(".tables can only be used with at most one select statement, use the dot notation: e.g. 'tableX.columnY'")
	}else if (length(clean.cols) == 1)
	{
	    temp.table <- use.tables
	    use.tables <- clean.cols
	    names(use.tables) <- temp.table
	}
	
    }else if (length(use.expr) == 0 && is.null(.tables))
    {
	stop("ERROR: Please either supply desired columns (columns(.data)) or specify valid table(s) in .tables ('tables(.data)')")
    }else{
	#attempt to figure out what the tables are from the specified columns...
	use.col.list <- lapply(use.expr, function(x)
			   {
				expr.only <- x$expr
				if (length(expr.only) == 1)
				{
				    return(as.character(expr.only))
				}else if (length(expr.only) == 3 && expr.only[[1]] == ":")
				{
				    return(sapply(expr.only[2:3], as.character))
				}else{
				    stop("ERROR: Accepted expression are of the form 'column' or 'column1':'columnN'")
				}
			   })
	
	inp.tab.list <- get.tables.from.vars(use.col.list)
	
	#where the nulls are non-input tables
	not.sup.tab <- sapply(inp.tab.list, function(x) all(is.null(x)))
	
	clean.cols <- sapply(1:length(use.expr), function(x)
				 {
				    if (is.null(inp.tab.list[[x]])==F)
				    {
					return(gsub(paste0(inp.tab.list[[x]], "\\."), "", deparse(use.expr[[x]]$expr)))
				    }else{
					return(deparse(use.expr[[x]]$expr))
				    }
				    
				 })
	
	if(any(not.sup.tab))
	{
	    tab.cols <- columns(.data)
	    
	    #if no tables are specified, need to first try to find a table which contains all the columns
	    tab.ord.mat <- sapply(tab.cols, function(x) sapply(use.col.list[not.sup.tab], function(y) {
		
		temp.match.rank <- rank(match(y, x), na.last=NA)
		
		if (length(temp.match.rank) > 0 && all(temp.match.rank == seq_along(y)))
		{
		    return(T)
		}else{
		    return(F)
		}
	    }))
	    
	    if (is.matrix(tab.ord.mat) == F)
	    {
		tab.ord.mat <- matrix(tab.ord.mat, nrow=1, ncol=length(tab.ord.mat), dimnames=list(NULL, names(tab.cols)))
	    }
	    
	    if (all(apply(tab.ord.mat, 1, function(x) sum(x) == 1)))
	    {
		sel.tabs <- apply(tab.ord.mat, 2, which)
		
		which.to.use <- sapply(sel.tabs, length) > 0
		
		##bugfix to address multiple columns from same table
		#unl.tabs <- unlist(sel.tabs[which.to.use])
		#use.tables <- clean.cols[not.sup.tab][unl.tabs]
		#names(use.tables) <- names(unl.tabs)
		
		unl.tabs <- stack(sel.tabs[which.to.use])
		use.tables <- sapply(sel.tabs[which.to.use], function(x) paste(clean.cols[not.sup.tab][x], collapse=","))
		#use.tables <- clean.cols[not.sup.tab][unl.tabs$values]
		#names(use.tables) <- as.character(unl.tabs$ind)
		
	    }else{
		#otherwise the columns would need to be contiguous between multiple tables
		##probably not going to be easy for user to specify, will throw an error for now...
		
		not.contig <- apply(tab.ord.mat, 1, any) == F
		
		if (any(not.contig))
		{
		    stop(paste("ERROR: Column(s):",paste(sapply(use.expr[not.sup.tab][not.contig], function(x) deparse(x$expr)), collapse=","), "are not contiguous in any table.  Try specifying a table: e.g. 'tableX.columnY'"))
		}else{
		    
		    multi.tab <- apply(tab.ord.mat, 1, function(x) sum(x) > 1)
		    
		    stop(paste("ERROR: Column(s):",paste(sapply(use.expr[not.sup.tab][multi.tab], function(x) deparse(x$expr)), collapse=","), "are found in multiple tables try specifying a table: 'tableX.columnY'"))
		}
		
	    }
	    
	    if(any(not.sup.tab==F))
	    {
		#a mix of specified and non-specified
		##so add in the additional tables, if the columns exist...
		
		if (any(sapply(inp.tab.list[not.sup.tab==F], length) != 1))
		{
		    stop("ERROR: Only a single table should be specified per statement.")
		}
		
		use.tabs <- sapply(inp.tab.list[not.sup.tab==F], "[", 1)
		
		valid.cols <- mapply(function(cols, tabs, cur.tab)
		       {
			    strip.cols <- gsub(paste0(cur.tab, "\\."), "", cols)
			    
			    temp.match.rank <- rank(match(strip.cols, tabs), na.last=NA)
		
			    if (length(temp.match.rank) > 0 && all(temp.match.rank == seq_along(cols)))
			    {
				return(T)
			    }else{
				return(F)
			    }
		       }, use.col.list[not.sup.tab==F], tab.cols[use.tabs], cur.tab=use.tabs)
		
		
		if (all(valid.cols) == F)
		{
		    stop(paste("ERROR: Provided column(s):", paste(sapply(use.expr[not.sup.tab==F][valid.cols==F], function(x) deparse(x$expr)), collapse=","), "could not be found in the specified tables"))
		}
		
		new.tab.list <- clean.cols[not.sup.tab==F]
		names(new.tab.list) <- use.tabs
		
		for(tab.name in names(new.tab.list))
		{
		  if(tab.name %in% names(use.tables))
		  {
		    use.tables[tab.name] <- paste(use.tables[tab.name], new.tab.list[tab.name], sep=",")
		  }else{
		    use.tables <- append(use.tables,new.tab.list[tab.name])
		  }
		}
	    }
	    
	}else{
	    unl.tabs <- stack(inp.tab.list)
	    unl.tabs$cols <- sapply(strsplit(as.character(unl.tabs$ind), "\\."), "[", 2)
	    
	    split.tabs <- split(unl.tabs, unl.tabs$values)
	    
	    use.tables <- sapply(split.tabs, function(x) paste(x$cols, collapse=","))
	  
	    #use.tables <- clean.cols
	    #names(use.tables) <- sapply(inp.tab.list, "[", 1)
	}
	
	#return(eval(parse(text=paste("select(join(.data, use.tables), ", paste(clean.cols, collapse=","), ")"))))
	return(select_(join(.data, use.tables), .dots=as.list(clean.cols)))
    }
    
    return(join(.data, use.tables))
}
