
#interface with igraph
tsl.to.graph <- function(tsl)
{
    graph.list <- lapply(tsl@tab.list, function(x)
           {
                cur.edges <- unlist(lapply(x$foreign.keys, "[[", "local.keys"))
                if (is.null(cur.edges))
                {
                    return(NULL)
                }
                else
                {
                    common.col <- intersect(x$db.cols, cur.edges)
                    use.fk <- sapply(x$foreign.keys, function(x) all(x$local.keys %in% common.col))
                    return(names(use.fk)[use.fk == TRUE])
                }
           })
    
    graph.comp.list <- graph.list[sapply(graph.list, is.null)==F]
    
    return(graph.data.frame(stack(graph.comp.list)))
}

get.starting.point <- function(tbsl, use.tables)
{
    tsl.graph <- tsl.to.graph(tbsl)
    sp.mat <- shortest.paths(tsl.graph, mode="all")
    sp.mat[is.infinite(sp.mat)] <- 0
    diag(sp.mat) <- NA
    
    if (use.tables %in% rownames(sp.mat) == F || use.tables %in% colnames(sp.mat) == F)
    {
	browser()
    }
    
    sp.mat <- sp.mat[use.tables, use.tables, drop=F]
    
    is.valid <- apply(sp.mat, 1, function(x) all(na.omit(x) > 0))
    
    valid.mat <- sp.mat[is.valid,,drop=F]
    
    dists <- apply(valid.mat, 1, function(x) sum(na.omit(x)))
    
    #assuming that the furthest table is at one end of the query path and 'should' make a good start point
    return(names(dists)[which.max(dists)])
}

get.shortest.query.path <- function(tbsl, start=NULL, finish=NULL, reverse=TRUE, undirected=TRUE)
{   
    tsl.graph <- tsl.to.graph(tbsl)
    
    if (missing(finish) || is.null(finish) || all(is.na(finish)))
    {
        finish <- V(tsl.graph)
    }
    
    if (undirected)
    {
        use.mode <- "all"
    }else{
        use.mode <- "out"
    }
    
    table.path <- get.shortest.paths(graph=tsl.graph,from=start,to=finish, mode = use.mode, weights = NULL, output="vpath",predecessors = FALSE, inbound.edges = FALSE)
    #do it backwards as that is how we want to merge it with the pd tables
    
    mask.query <- lapply(table.path$vpath, function(x)
			 {
			    if (reverse==TRUE)
			    {
				return(rev(V(tsl.graph)[x]$name))
			    }
			    else
			    {
				return(V(tsl.graph)[x]$name)
			    }
			 })
    
    return(mask.query)
}





#need to add validitity checks to default.search.cols to below function
valid.TableSchemaList <- function(object)
{
    if (is.null(names(object@tab.list)))
    {
        return("The supplied tab.list needs to have names")
    }
    
    valid.list <- sapply(object@tab.list, function(x)
           {
                if (is.null(names(x)) == TRUE || all(names(x) %in% c( "db.cols","db.schema", "db.constr", "dta.func", "should.ignore", "foreign.keys")) == FALSE)
                {
                    return(FALSE)
                }
                else
                {
                    if (length(x$db.schema) == length(x$db.cols))
                    {
                        if (is.null(x$foreign.keys))
                        {
                            return(TRUE)
                        }
                        else if (class(x$foreign.keys) == "list" && is.null(names(x$foreign.keys)) == FALSE)
                        {
                            return(all(sapply(x$foreign.keys, function(x)
                                   {
                                        return(all(names(x) %in% c("local.keys", "ext.keys")))
                                   })))
                        }
                        else
                        {
                            return(FALSE)
                        }
                    }
                    else
                    {
                        return(FALSE)
                    }
                }
           })
    
    if (all(valid.list) == TRUE)
    {
        return(TRUE)
    }
    else
    {
        return(paste("Invalid input for: ", names(valid.list)[valid.list == FALSE]))
    }
}

setClass(Class="TableSchemaList", representation=list(tab.list="list"), prototype=prototype(tab.list=list()), validity=valid.TableSchemaList)

#need to fix me...
TableSchemaList <- function(tab.list=NULL)
{
    return(new("TableSchemaList", tab.list=tab.list))
}

setMethod("show", signature("TableSchemaList"), function(object)
          {
                message(paste("TableSchemaList containing", length(object), "tables"))
		
		if (length(object) > 0)
		{
		    for (i in 1:length(object))
		    {
			num.cols <- length(object@tab.list[[i]]$db.cols)
			col.val <- ifelse(num.cols == 1, "Column", "Columns")
			message(paste("   ", "-", names(object@tab.list)[i], "(", num.cols, col.val,")"))
			if (is.null(object@tab.list[[i]]$foreign.keys) ==F)
			{
			    for (j in names(object@tab.list[[i]]$foreign.keys))
			    {
				all.loc.keys <- object@tab.list[[i]]$foreign.keys[[j]]$local.keys
				
				#check if the loc.keys are still retained in the table or only used for merging purposes
				
				loc.keys <- all.loc.keys[all.loc.keys %in% object@tab.list[[i]]$db.cols]
				
				if (length(loc.keys) > 0)
				{
				     if (length(loc.keys) == 1)
				    {
					message(paste("      - Column", loc.keys, "is derived from table", j))
				    }else{
					message(paste("      - Columns", paste(loc.keys, collapse=","), "are derived from table", j))
				    }
				}
    
			    }
			}
			
			if (object@tab.list[[i]]$db.constr != "")
			{
			    constr.cols <- gsub("\\s+", "", regmatches(object@tab.list[[i]]$db.constr, regexec("\\((.+)\\)", object@tab.list[[i]]$db.constr))[[1]][2])
			    message(paste("      - Uniqueness constraints on:", constr.cols))
			}
		    }
		}
		
		
          })

#may need to use BiocGenerics at some point...

setGeneric("append")
setMethod("append", signature("TableSchemaList", "TableSchemaList"), function(x, values, after=length(x))
	  {
		 lengx <- length(x)
		    
		    if (!after) 
			return(new("TableSchemaList", tab.list=c(values@tab.list, x@tab.list)))
		    else if (after >= lengx) 
			return(new("TableSchemaList", tab.list=c(x@tab.list, values@tab.list)))
		    else return(new("TableSchemaList", tab.list=c(x@tab.list[1L:after], values@tab.list, x@tab.list[(after + 1L):lengx])))
		
	  })

#setGeneric("length", def=function(x), standardGeneric("length"))
setMethod("length", signature("TableSchemaList"), function(x)
	  {
		return(length(x@tab.list))
	  })

subset.TableSchemaList <- function(x, table.name, ...)
{
    if (all(table.name %in% names(x@tab.list)) == FALSE)
            {
                stop("ERROR: Only valid names can be used for subsetting")
            }
            
            return(new("TableSchemaList", tab.list=x@tab.list[table.name]))
}

makeSchemaFromFunction <- function(dta.func,name,...)
{
    if (missing(dta.func) || is.null(dta.func) || is.function(dta.func) == F)
    {
	stop("ERROR: Please supply a function for 'func'")
    }
    
    pass.objs <- list(...)
    
    if (length(pass.objs) == 0)
    {
	stop("ERROR: Please supply object(s) to apply 'func' to")
    }
    
    makeSchemaFromData(do.call(dta.func, pass.objs), name=name, dta.func=dta.func)
}

makeSchemaFromData <- function(tab.df, name=NULL, dta.func=NULL)
{
  if (missing(name) || is.null(name) || is.na(name))
  {
      stop("ERROR: Please supply a name for the table")
  }
  
  if (class(tab.df) != "data.frame")
    {
	stop("ERROR: tab.df needs to be a data.frame")
    }
  
  if (is.null(names(tab.df)))
  {
      stop("ERROR: tab.df needs to have column names")
  }
  
  if (valid.db.names(names(tab.df)) == F)
  {
      stop("ERROR: The names of the supplied data.frame need to be modified for the database see correct.df.names")
  }
  
  if (missing(dta.func) || is.null(dta.func))
    {
	use.func <- eval(parse(text=paste("function(",name,")", name)))
    }else{
	use.func <- dta.func
    }
  
  cur.list <- list(db.cols=character(0), db.schema=character(0), db.constr="", dta.func=use.func , should.ignore=T, foreign.keys=NULL)
  
    cur.list$db.constr <- ""
    cur.list$db.cols <- paste0(name, "_ind")
    cur.list$db.schema <- "INTEGER PRIMARY KEY AUTOINCREMENT"
 
  
  cur.list$db.cols <- append(cur.list$db.cols, names(tab.df))
  
  cur.list$db.schema <- append(cur.list$db.schema, determine.db.types(tab.df))
  
  tab.list <- list(cur.list)
  names(tab.list) <- name
  
  return(new("TableSchemaList", tab.list=tab.list))
}

correct.df.names <- function(tab.df)
{
    names(tab.df) <- make.db.names.default(names(tab.df))
    return(tab.df)
}

valid.db.names <- function(inp.names)
{
    conv.names <- make.db.names.default(inp.names)

    if (length(setdiff(inp.names, conv.names)) > 0)
    {
	return(FALSE)
    }
    else
    {
	return(TRUE)
    }
}

character.to.type <- function(val.class)
{
    char.val <- as.character(val.class)
    
    exist.nas <- sum(is.na(char.val))
    is.numeric.type <- suppressWarnings(as.numeric(char.val))
    
    if (sum(is.na(is.numeric.type)) > exist.nas)
    {
        return("TEXT")
    }
    else
    {
        if(any(grepl("\\d*\\.\\d*",char.val)))
        {
            return("NUMERIC")
        }
        else
        {
            return("INTEGER")
        }
    }
}

basic.integer <- function(x)
{
    return("INTEGER")
}


basic.text <- function(x)
{
    return("TEXT")
}

determine.db.types <- function(dta)
{
    return(sapply(names(dta), function(x)
           {
                return(switch(class(dta[,x]), character=basic.text, factor=character.to.type,
                    numeric=character.to.type, integer=basic.integer, basic.text)(dta[,x]))
           }))
}


return.element <- function(use.obj, name)
{
    return(sapply(use.obj@tab.list, "[[", name))
}

setGeneric("columns", def=function(obj,...) standardGeneric("columns"))
setMethod("columns", signature("TableSchemaList"), function(obj)
	  {
	    
	    ret.list <- lapply(schemaNames(obj), function(x)
		   {
			colNames(obj, x)
		   })
	    
	    names(ret.list) <- schemaNames(obj)
	    
	    return(ret.list)
	  })


get.tables.from.vars <- function(col.list)
{
    tb.list <- lapply(col.list, function(x)
	   {
		split.x <- strsplit(x, "\\.")[[1]]
	    
		if (length(split.x) == 1)
		{
		    return(NULL)
		    
		}else if(length(split.x) == 2){
		    
		    return(split.x[1])
		    
		}else{
		    stop("ERROR: Unexpected format of supplied variable")
		}
	   })
    
    return(tb.list)
}


setGeneric("foreignExtKeyCols", def=function(obj, ...) standardGeneric("foreignExtKeyCols"))
setMethod("foreignExtKeyCols", signature("TableSchemaList"), function(obj, table.name)
          {
	    lapply(return.element(obj, "foreign.keys")[[table.name]], function(x)
		   {
			return(x$ext.keys)
		   })
          })

setGeneric("foreignLocalKeyCols", def=function(obj, ...) standardGeneric("foreignLocalKeyCols"))
setMethod("foreignLocalKeyCols", signature("TableSchemaList"), function(obj, table.name, join.table=NULL)
          {
	    ret.list <- lapply(return.element(obj, "foreign.keys")[[table.name]], function(x)
		   {
			return(x$local.keys)
		   })
	    
	    if (is.null(join.table))
	    {
		return(ret.list)
	    }else if (length(join.table) == 1){
		return(ret.list[[join.table]])
	    }else{
		return(ret.list[join.table])
	    }
            
          })

setGeneric("foreignExtKeySchema", def=function(obj, ...) standardGeneric("foreignExtKeySchema"))
setMethod("foreignExtKeySchema", signature("TableSchemaList"), function(obj, table.name)
          {
	    
	    ext.keys <- foreignExtKeyCols(obj, table.name)
	    
	    ret.list <- lapply(names(ext.keys), function(x,e.k)
		   {
			colSchema(obj, x, mode="normal")[match(e.k[[x]], colNames(obj, x, mode="normal"))]
		   }, ext.keys)
	    
	    names(ret.list) <- names(ext.keys)
	    return(ret.list)
	    
          })

setGeneric("bindDataFunction", def=function(obj, ...) standardGeneric("bindDataFunction"))
setMethod("bindDataFunction", signature("TableSchemaList"), function(obj, table.name, bind.vals, mode=c("normal", "merge"))
        {
            table.mode <- match.arg(mode)
            
	    cur.func <- return.element(subset(obj, table.name), "dta.func")[[1]]
	    
	    cur.forms <- formals(cur.func)
	    
	    if (all(names(cur.forms) %in% names(bind.vals)) == F)
	    {
		stop(paste("ERROR: Cannot find variables with names:", paste(names(cur.forms), collapse=",")))
	    }
	    
            vcf.dta <- do.call(cur.func, bind.vals[names(cur.forms)])
            
            cur.cols <- colNames(obj, table.name, mode=table.mode)
            cur.schema <- colSchema(obj, table.name, mode=table.mode)
            
            #don't need to supply columns which are autoincremented, they will be automatically added to the data.frame
            auto.col <- cur.cols[cur.schema == "INTEGER PRIMARY KEY AUTOINCREMENT"]
            #stopifnot(length(auto.col) == 1)
            
            diff.cols <- setdiff(cur.cols, colnames(vcf.dta))
            
            if (length(diff.cols) == 0)
            {
                return(vcf.dta[,cur.cols])
            }
            else if (length(auto.col) == 1 && length(diff.cols) == 1 && diff.cols == auto.col)
            {
                temp.vcf.dta <- cbind(vcf.dta, NA_integer_)
                names(temp.vcf.dta) <- c(names(vcf.dta), auto.col)
                return(temp.vcf.dta[,cur.cols])
            }
            else
            {
                nf.cols <- setdiff(cur.cols, colnames(vcf.dta))
                stop(paste("ERROR: Cannot find column(s)", paste(nf.cols, collapse=",")))
            }
            
        })

setGeneric("shouldIgnore", def=function(obj, ...) standardGeneric("shouldIgnore"))
setMethod("shouldIgnore", signature("TableSchemaList"), function(obj, table.name)
        {
            return(return.element(subset(obj, table.name), "should.ignore"))
        })

setGeneric("shouldMerge", def=function(obj, ...) standardGeneric("shouldMerge"))
setMethod("shouldMerge", signature("TableSchemaList"), function(obj, table.name=NULL)
        {
            if (missing(table.name) || is.null(table.name))
            {
                sub.obj <- obj
            }
            else
            {
                sub.obj <- subset(obj, table.name)
            }
            
	    #this is a fix on 8-18-2014 to allow a relationship between two tables to be specified using the same column(s)
	    
	    return(any(sapply(return.element(sub.obj, "foreign.keys"), function(x)
		   {
			if (is.null(x))
			{
			    return(F)
			}
			else if (length(intersect(x$local.keys, x$ext.keys)) == length(union(x$local.keys, x$ext.keys)))
			{
			    return(F)
			}
			else
			{
			    return(T)
			}
		   })))
	    
        })

setGeneric("schemaNames", def=function(obj, ...) standardGeneric("schemaNames"))
setMethod("schemaNames", signature("TableSchemaList"), function(obj)
        {
            names(obj@tab.list)
        })

setGeneric("tableConstr", def=function(obj, ...) standardGeneric("tableConstr"))
setMethod("tableConstr", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
        {
            table.mode <- match.arg(mode)
            
            ret.el <- return.element(subset(obj, table.name), "db.constr")
            
            if (is.null(ret.el) || is.na(ret.el) || table.mode == "merge")
            {
                return("")
            }
            else
            {
                return(ret.el)
            }
        })

setGeneric("tableName", def=function(obj, ...) standardGeneric("tableName"))
setMethod("tableName", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                table.mode <- match.arg(mode)
            
                if (table.name %in% names(obj@tab.list) == FALSE)
                {
                    stop("ERROR: Invalid table.name supplied")
                }
                
                return(switch(table.mode, normal=table.name, merge=paste(table.name, "temp", sep="_")))
          })

#need to make colNames and colSchema consistent and deal with the direct keys
setGeneric("colNames", def=function(obj, ...) standardGeneric("colNames"))
setMethod("colNames", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                .get.names.schema(obj, table.name, mode, type="cols")
          })

setGeneric("colSchema", def=function(obj, ...) standardGeneric("colSchema"))
setMethod("colSchema", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                .get.names.schema(obj, table.name, mode, type="schema")
          })

.get.names.schema <- function(obj, table.name, mode=c("normal", "merge"), type=c("cols", "schema"))
{
    table.mode <- match.arg(mode)
    type <- match.arg(type)
    sub.obj <- subset(obj, table.name)
    
    base.schema <- as.character(return.element(sub.obj, "db.schema"))
    base.cols <- as.character(return.element(sub.obj, "db.cols"))
    
    if (table.mode == "normal")
    {
	return(switch(type, cols=base.cols, schema=base.schema))
    }
    else
    {
	foreign.schema <- foreignExtKeySchema(obj, table.name)
	foreign.cols <- foreignExtKeyCols(obj, table.name)
	local.cols <- foreignLocalKeyCols(obj, table.name)
	
	is.direct.key <- isDirectKey(obj, table.name)
	
	direct.keys <- as.character(unlist(foreign.cols[is.direct.key]))
	unl.f.c <- as.character(unlist(foreign.cols))
	unl.l.c <- as.character(unlist(local.cols))
	
	rm.cols <- setdiff(unl.l.c, direct.keys)
	
	#if there are any direct keys that are derived from other tables 'local.keys' then remove as well
	
	derived.direct <- direct.keys %in% unlist(local.cols[is.direct.key==F])
	
	if (any(derived.direct))
	{
	    rm.cols <- append(rm.cols, direct.keys[derived.direct])
	}
	
	keep.foreign.cols <- (unl.f.c %in% rm.cols == F) & (unl.f.c %in% base.cols == F)
	
	#keep.foreign.cols <- unlist(mapply(function(x,y){
	#			    #also remove the columns that will be present in the final table but not part of the initial table but keep the direct keys
	#			    rm.cols <- setdiff(y, direct.keys)
	#			    should.keep <- (x %in% rm.cols == F) & (x %in% base.cols == F)
	#			    return(should.keep)
	#		       }, foreign.cols, local.cols))
	
	rm.base.cols <- base.cols %in% rm.cols#setdiff(unlist(local.cols), direct.keys)
	rm.base.cols <- rm.base.cols | base.schema == "INTEGER PRIMARY KEY AUTOINCREMENT"
	
	if (type == "cols")
	{
	    return(as.character(c(unl.f.c[keep.foreign.cols], base.cols[rm.base.cols == F])))
	}
	else
	{
	    return(as.character(c(unlist(foreign.schema)[keep.foreign.cols], base.schema[rm.base.cols == F])))
	}
    }
    
    
}

setGeneric("createTable", def=function(obj, ...) standardGeneric("createTable"))
setMethod("createTable", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                table.mode <- match.arg(mode)
                if (shouldMerge(obj, table.name) == TRUE && table.mode == "merge")
                {
                    temp.str <- "TEMPORARY"
                }
                else if (table.mode == "normal")
                {
                    temp.str <- ""
                }
                else
                {
                    stop("ERROR: Cannot generate statement with mode set to 'merge' and a NULL foreign.key element")
                }
                
                use.cols <- colNames(obj, table.name, mode=table.mode)
                use.schema <- colSchema(obj, table.name, mode=table.mode)
                
                tab.constr <- tableConstr(obj, table.name, mode=table.mode)
                tab.constr <- ifelse(tab.constr == "", tab.constr, paste0(",", tab.constr))
                
                return(paste("CREATE",temp.str,"TABLE", tableName(obj, table.name, table.mode), "(", paste(paste(use.cols, use.schema), collapse=","), tab.constr, ")"))
          })

setGeneric("directKeys", def=function(obj, ...) standardGeneric("directKeys"))
setMethod("directKeys", signature("TableSchemaList"), function(obj, table.name)
	  {
		all.keys <- return.element(obj, "foreign.keys")[[table.name]]
		
		is.direct.keys <- sapply(all.keys, function(x) length(intersect(x$local.keys, x$ext.keys)) == length(union(x$local.keys, x$ext.keys)))
		
		return(unique(as.character(unlist(all.keys[is.direct.keys]))))
	  })

setGeneric("isDirectKey", def=function(obj, ...) standardGeneric("isDirectKey"))
setMethod("isDirectKey", signature("TableSchemaList"), function(obj, table.name)
	  {
		all.keys <- return.element(obj, "foreign.keys")[[table.name]]
		
		return(sapply(all.keys, function(x) length(intersect(x$local.keys, x$ext.keys)) == length(union(x$local.keys, x$ext.keys))))
	  })

setGeneric("mergeStatement", def=function(obj, ...) standardGeneric("mergeStatement"))
setMethod("mergeStatement", signature("TableSchemaList"), function(obj, table.name)
          {
	    
                #currently, probably the temporary table
                cur.db <- tableName(obj, table.name, mode="merge")
                #table trying to create
                target.db <- tableName(obj, table.name, mode="normal")
                
		target.cols <- colNames(obj, table.name, mode="normal")
                target.schema <- colSchema(obj, table.name, mode="normal")
                #remove the autoincrement column first
                target.cols <- target.cols[target.schema != "INTEGER PRIMARY KEY AUTOINCREMENT"]
		
                #create the join statement using the foreign.keys slot
                
                fk <- return.element(obj, "foreign.keys")[[table.name]]
                
                if (is.null(fk))
                {
                    stop("ERROR: Cannot generate statement if the foreign key element is NULL")
                }
                
		#don't merge on the tables that have direct keys
		is.direct.keys <- sapply(fk, function(x) length(intersect(x$local.keys, x$ext.keys)) == length(union(x$local.keys, x$ext.keys)))
		
		fk <- fk[is.direct.keys == F] 
		
                keys <- sapply(fk, function(y)
                       {
                            return(paste(y$ext.keys, collapse=","))
                       })
		
		join.statement <- paste(paste("JOIN", names(keys), "USING", paste0("(", keys,")")), collapse=" ")
		
		#in case columns besides the keys are duplicated, make sure the select statement refers to the appropriate table
                
		names(target.cols) <- rep(cur.db, length(target.cols))
		
		for(i in names(fk))
		{
		    join.tab.cols <- target.cols %in% fk[[i]]$local.keys
		    
		    if (any(join.tab.cols))
		    {
			names(target.cols)[join.tab.cols] <- i
		    }
		}
		
		plain.targs <- paste(target.cols, collapse=",")
		
		paste.targs <- paste(paste(names(target.cols), target.cols, sep="."), collapse=",")
		
                if (shouldIgnore(obj, table.name))
                {
                    ignore.str <- "OR IGNORE"
                }
                else
                {
                    ignore.str <- ""
                }
                
                return(paste("INSERT",ignore.str,"INTO", target.db, "(", plain.targs,") SELECT", paste.targs,"FROM", cur.db , join.statement))
          })

setGeneric("insertStatement", def=function(obj, ...) standardGeneric("insertStatement"))
setMethod("insertStatement", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                table.mode <- match.arg(mode)
                
                if (shouldMerge(obj, table.name) == FALSE && table.mode == "merge")
                {
                    stop("ERROR: Cannot run statement with mode 'merge' and a NULL foreign.key element")
                }
                
                use.cols <- colNames(obj, table.name, mode=table.mode)
                
                if (shouldIgnore(obj, table.name))
                {
                    ignore.str <- "OR IGNORE"
                }
                else
                {
                    ignore.str <- ""
                }
                
                return(paste("INSERT",ignore.str,"INTO", tableName(obj, table.name, table.mode), "VALUES (", paste(paste0(":", use.cols), collapse=","), ")"))
          })

.get.model.side <- function(form, side=c("left", "right"), num.components=3)
{    
    side <-  match.arg(side)
    
    char.form <- as.character(form)
    
    if (length(char.form) == num.components)
    {
        use.ind <- switch(side, right=num.components, left=num.components-1)
        
        split.lhs <- strsplit(char.form[use.ind], "\\s+\\+\\s+")[[1]]
        return(split.lhs)
    }
    else
    {
        stop("ERROR: Unexpected parsing for input formula")
    }
}



rhs <- function(form)
{
    .get.model.side(form, "right")
}

lhs <- function(form)
{
    .get.model.side(form, "left")
}

#check if a variable is a reference to a primary key ie: .table
#check to be sure that the values are in 'to'
.resolve.rhs.fk <- function(obj, cur.rhs, cur.tab)
{
    if (any(grepl("^\\.", cur.rhs, perl=T)))
    {
	ref.pos <- grepl("^\\.", cur.rhs, perl=T)
	ref.cols <- cur.rhs[ref.pos]
	ref.tables <- sub("\\.", "", ref.cols)
	
	if (all(ref.tables %in% tables(obj)) == F)
	{
	    stop(paste("ERROR: invalid tables specified:", paste(setdiff(ref.tables, tables(obj)), collapse=",")))
	}
	
	valid.refs <- as.character(sapply(ref.tables, function(x) obj@tab.list[[x]]$db.cols[obj@tab.list[[x]]$db.schema == "INTEGER PRIMARY KEY AUTOINCREMENT"]))
	
	cur.rhs[ref.pos] <- valid.refs
    }
    else if (all(cur.rhs %in% obj@tab.list[[cur.tab]]$db.cols) == F)
    {
	stop("ERROR: All values on the rhs of the formula need to be in the 'to' table's columns")
    }
    
    return(cur.rhs)
}

setGeneric("constraint<-", def=function(obj, table.name,..., value) standardGeneric("constraint<-"))
setReplaceMethod("constraint", signature("TableSchemaList"), function(obj, table.name,value, should.ignore=T, constr.name=NULL)
		 {
		    
		    if (missing(table.name) || is.null(table.name) || length(table.name) != 1  || is.na(table.name) || all(table.name %in% tables(obj))==F)
		    {
			stop("ERROR: please supply a single valid table for 'table.name'")
		    }
		    
		    if (missing(constr.name) || is.null(constr.name))
		    {
			constr.name <- paste(table.name, "idx", sep="_")
		    }
		    
		    if (length(should.ignore) != 1 || is.logical(should.ignore) == F)
		    {
			stop("ERROR: should.ignore should be a single logial value")
		    }
		    
		    if (is.null(value))
		    {
			
			obj@tab.list[[table.name]]$db.constr <- ""
			
		    }else{
			cur.rhs <- .get.model.side(value, "right", num.components=2)
		    
			cur.rhs <- .resolve.rhs.fk(obj, cur.rhs, table.name)
			
			obj@tab.list[[table.name]]$db.constr <- paste("CONSTRAINT",constr.name,"UNIQUE (",paste(cur.rhs, collapse=","),")")
		    }
		    
                    obj@tab.list[[table.name]]$should.ignore <- should.ignore
		    
		    validObject(obj)
                    return(obj)
		    
		 })

setGeneric("relationship<-", def=function(obj, ..., value) standardGeneric("relationship<-"))
setReplaceMethod("relationship", signature("TableSchemaList"), function(obj, value, from, to)
                 {
                    #check to be sure 'from' comes before 'to'
                    
                    from.pos <- which(names(obj@tab.list) == from)
                    to.pos <-  which(names(obj@tab.list) == to)
                    
                    if (to.pos <= from.pos)
                    {
                        stop("ERROR: The 'from' table needs to be specified before the 'to' table")
                    }
		    
                    cur.lhs <- lhs(value)
                    cur.rhs <- rhs(value)
                    
                    if (length(cur.lhs) == 1 && cur.lhs == ".")
                    {
                        #check to see if there is an auto-incremented primary key for from
                        which.prim <- which(obj@tab.list[[from]]$db.schema == "INTEGER PRIMARY KEY AUTOINCREMENT")
                        if (length(which.prim) == 1)
                        {
                            cur.lhs <- obj@tab.list[[from]]$db.cols[which.prim]
                        }
                        else
                        {
                            stop("ERROR: Can't specify '.' in the formula when a autoincremented primary key is not available")
                        }
                    }
                    else
                    {
                        #check to be sure that the values are in from
                        if (all(cur.lhs %in% obj@tab.list[[from]]$db.cols) == F)
                        {
                            stop("ERROR: All values on the lhs of the formula need to be in the 'from' table's columns")
                        }
                    }
		    
		    cur.rhs <- .resolve.rhs.fk(obj, cur.rhs, to)
                    
                    #additionally, they all need to be in from's as well
                    
                    if (all(cur.rhs %in% obj@tab.list[[from]]$db.cols) == F)
                    {
                        stop("ERROR: All values on the rhs of the formula need to be in the 'from' table's columns")
                    }
                
                    #if everything looks good, then add in the relationships to the to table
                    
                    use.fk <- list(list(local.keys=cur.lhs, ext.keys=cur.rhs))
                    names(use.fk) <- from
                    
                    if (is.null(obj@tab.list[[to]]$foreign.keys))
                    {
                        obj@tab.list[[to]]$foreign.keys <- use.fk
                    }
                    else
                    {
                        obj@tab.list[[to]]$foreign.keys <- append(obj@tab.list[[to]]$foreign.keys, use.fk)
                    }
                    
                    #also modify db.cols and db.schema to reflect the new keys/relationships if a non one-to-one relationship is specified
                    
		    if (length(intersect(cur.rhs, cur.lhs)) != length(union(cur.rhs, cur.lhs)))
		    {
			#if any of the previous keys are one-to-one keys then keep them as part of the table
			
			direct.keys <- directKeys(obj, to)
			
			cur.rhs <- setdiff(cur.rhs, direct.keys)
			
			which.rhs <- which(obj@tab.list[[to]]$db.cols %in% cur.rhs)
			obj@tab.list[[to]]$db.cols <- obj@tab.list[[to]]$db.cols[-which.rhs]
			obj@tab.list[[to]]$db.schema <- obj@tab.list[[to]]$db.schema[-which.rhs] 
			
			obj@tab.list[[to]]$db.cols <- append(obj@tab.list[[to]]$db.cols, cur.lhs)
			
			#add in what the schema is in from, cleaning a little for integer autoincrement
			
			curm <- match(cur.lhs, obj@tab.list[[from]]$db.cols)
			
			obj@tab.list[[to]]$db.schema <- append(obj@tab.list[[to]]$db.schema, sapply(strsplit(obj@tab.list[[from]]$db.schema[curm], "\\s+"), "[[", 1))
		    }
		    
                    validObject(obj)
                    return(obj)
                 })


read.database.tables <- function(db.name, num.rows=10)
{
    temp.con <- dbConnect(SQLite(), db.name)
    
    tab.names <- dbListTables(temp.con)
    
    tab.list <- lapply(tab.names, function(x) dbGetQuery(temp.con, paste0('SELECT * FROM ', x, ' LIMIT ', num.rows)))
    
    names(tab.list) <- tab.names
    
    tab.list <- tab.list[-which(names(tab.list) == "sqlite_sequence")]
    
    dbDisconnect(temp.con)
    
    return(tab.list)
}