setClass(Class="Database", representation=list(tbsl="TableSchemaList", db.file="character"))

setMethod("show", signature("Database"), function(object)
	  {
	    message("Database Instance")
	    message(paste0("File: '", dbFile(object), "'"))
	    show(schema(object))
	  })

setGeneric("schema", def=function(obj,...) standardGeneric("schema"))
setMethod("schema", signature("Database"), function(obj)
	  {
	    return(obj@tbsl)
	  })

setGeneric("dbFile", def=function(obj,...) standardGeneric("dbFile"))
setMethod("dbFile", signature("Database"), function(obj)
	  {
	    return(obj@db.file)
	  })

setGeneric("tables", def=function(obj,...) standardGeneric("tables"))
setMethod("tables", signature("TableSchemaList"), function(obj)
	  {
	    return(schemaNames(obj))
	  })
setMethod("tables", signature("Database"), function(obj)
	  {
	    return(tables(schema(obj)))
	  })

setMethod("columns", signature("Database"), function(obj)
	  {
	    cur.schema <- schema(obj)
	    
	    return(columns(cur.schema))
	  })

Database <- function(tbsl, db.file)
{
    if (class(tbsl) != "TableSchemaList")
    {
	stop("ERROR: tbsl needs to be an instance of class TableSchemaList")
    }
    
    if ((is.character(db.file) && length(db.file) == 1)==F)
    {
	stop("ERROR: db.file needs to be a single path to a file")
    }
    
    #just want to make sure there is an available DB file...somewhat wasteful
    if (file.exists(db.file) == F)
    {
	temp.con <- dbConnect(SQLite(), db.file)
	dbDisconnect(temp.con)
    }
    
    return(new("Database", tbsl=tbsl, db.file=db.file))
}

#where cur.table is the index in use.path which is a character vector of tables and obj is a tableschemalist
get.join.keys <- function(cur.table, use.path, obj, ancil.tables)
{
	#maybe it simply needs to be if they share a direct key...
		    
	common.key <- intersect(directKeys(schema(obj), use.path[cur.table]), directKeys(schema(obj), use.path[cur.table+1]))
	
	if(length(common.key) > 0){
	 
		add.keys <- common.key
	}else{
		add.keys <- NULL
	}
	
	if (is.null(ancil.tables) == F && use.path[cur.table] %in% names(ancil.tables))
	{
		temp.ancil.keys <- lapply(ancil.tables[[use.path[cur.table]]], function(x) get.join.keys(1, c(use.path[cur.table+1], x) , obj, NULL))
		
		add.keys <- append(add.keys, unname(unlist(temp.ancil.keys)))
	}
	
	#finally the direct keys from one table to the next
	for.join <- foreignLocalKeyCols(schema(obj), use.path[cur.table], use.path[cur.table+1])
	
	if (is.null(for.join))
	{
		back.join <- foreignLocalKeyCols(schema(obj), use.path[cur.table+1], use.path[cur.table])
	 
		if (is.null(back.join))
		{
		    stop("ERROR: Cannot determine join structure")
		}
		else
		{
		    #the as.character has to do with inner_join and company have different semantics based on named versus unnamed
		    return(unique(as.character(append(back.join, add.keys))))
		}
	}
	else
	{
		return(unique(as.character(append(for.join, add.keys))))
	}
}

#still under construction, need to deal with multiple tables and possibly outer joins and such
setGeneric("join", def=function(obj, ...) standardGeneric("join"))
setMethod("join", signature("Database"), function(obj, needed.tables)
	  {
	    if (is.character(needed.tables) == F || (all(needed.tables %in% tables(obj))==F && all(names(needed.tables %in% tables(obj))) == F))
	    {
			stop("ERROR: needed.tables needs to be a character vector corresponding to table names")
	    }
	    
	    if (length(needed.tables) > 1)
	    {
			#use the TBSL object to determine how to join the tables and create a temporary table
			if (is.null(names(needed.tables)))
			{
				start.node <- get.starting.point(schema(obj), needed.tables)
			}else{
				start.node <- get.starting.point(schema(obj), names(needed.tables))
			}
			
			table.path <- get.shortest.query.path(schema(obj), start=start.node, finish=NULL, reverse=F, undirected=T)
			
			valid.path <- sapply(table.path, function(x) all(needed.tables %in% x) || (is.null(names(needed.tables)) == F && all(names(needed.tables) %in% x)))
			
			ancil.tables <- NULL
			
			if (all(valid.path == F))
			{
				
				
				tsl.graph <- tsl.to.graph(schema(obj))
				num.triang<- adjacent.triangles(tsl.graph)
				
				#if (any(num.triang > 0))
				#{
				#	browser()
				#}
				
				if (all(num.triang > 0) && length(num.triang) == 3)
				{
					min.tree <- minimum.spanning.tree(tsl.graph)
					edge.mat <- get.edges(min.tree, E(min.tree))
					
					use.path <- V(min.tree)[edge.mat[1,]]$name
					ancil.tables <- list(V(min.tree)[edge.mat[2,2]]$name)
					names(ancil.tables) <- V(min.tree)[edge.mat[2,1]]$name
					
				}else{
					#need to add in tables that are not part of the main path put should be added in per use.tables
					#also
					longest.table.path <- table.path[[which.max(sapply(table.path, length))]]
					
					if (is.null(names(needed.tables)))
					{
						lo.tables <- setdiff(needed.tables, longest.table.path)
					}else{
						lo.tables <- setdiff(names(needed.tables), longest.table.path)
					}
					
					#figure out if all the lo.tables can be joined directly to one or more tables on the longest.table.path
					
					temp.ancil.tables <- lapply(lo.tables, function(x)
									   {
											temp.sp <- get.shortest.query.path(schema(obj), start=x, finish=longest.table.path, reverse=F, undirected=T)
											temp.sp.lens <- sapply(temp.sp, length)
											
											if (any(temp.sp.lens == 2))
											{
												return(sapply(temp.sp[temp.sp.lens == 2], "[", 2))
											}else{
												stop(paste("ERROR: Cannot determine how to join table:",x,"query cannot be carried out"))
											}
									   })
					names(temp.ancil.tables) <- lo.tables
					
					stacked.ancil <- stack(temp.ancil.tables)
					ancil.tables <- split(as.character(stacked.ancil$ind), stacked.ancil$values)
					use.path <- longest.table.path
				}
				
			 
			}else{
				min.valid.path <- which.min(sapply(table.path[valid.path], length))
				use.path <- table.path[valid.path][[min.valid.path]]
				ancil.tables <- NULL
			}
		
		#the joining needs to take into account not just the direct keys from one table to the next but also the necessary
		#keys if one table has already been merged to another as well as any keys that are shared between the tables that
		#were derived from a downstream table
		
		join.cols <- lapply(1:(length(use.path)-1), function(x) {
		    
		    get.join.keys(x, use.path, obj, ancil.tables)
		})
		
		
		if (is.null(ancil.tables) == F)
		{
			#for each of the ancillary tables, join them in a piecewise fashion to their respective use.path tables
			ancil.join.cols <- mapply(function(a.tabs, tab){
				
				temp.keys <- lapply(1:length(a.tabs), function(x) get.join.keys(x, c(tab, a.tabs), obj, NULL))
				
				#temp.keys <- lapply(a.tabs, function(x)
				#	  {
				#		get.join.keys(1, c(x, tab), obj, NULL)
				#	  })
				names(temp.keys) <- a.tabs
				return(temp.keys)
			}, ancil.tables, names(ancil.tables), SIMPLIFY=F)
			
		}else{
			ancil.join.cols <- NULL
		}
		
		#now using dplyr::inner_join(x,y,by=NULL)
		
		src.db <- src_sqlite(dbFile(obj), create = F)
		
		if (is.null(names(needed.tables)))
		{
		    all.tab <- tbl(src.db, use.path[1])
			i <- 1
		    if (is.null(ancil.tables) == F && use.path[i] %in% names(ancil.tables))
			{
				for(j in ancil.tables[[use.path[i]]])
				{
					all.tab <- inner_join(all.tab, tbl(src.db,j), by=ancil.join.cols[[use.path[i]]][[j]])
				}
			}
		    use.path <- use.path[-1]
		    rm(i)
		    
		    for(i in seq_along(use.path))
		    {
				if (is.null(ancil.tables) == F && use.path[i] %in% names(ancil.tables))
				{
					for(j in ancil.tables[[use.path[i]]])
					{
						all.tab <- inner_join(all.tab, tbl(src.db,j), by=ancil.join.cols[[use.path[i]]][[j]])
					}
				}
				
				all.tab <- inner_join(all.tab, tbl(src.db,use.path[i]), by=join.cols[[i]])
		    }
		    
		}else{
		    
		    .get.select.cols <- function(tab, tab.exp, nec.cols, src.db)
		    {
				nec.cols <- nec.cols[!is.na(nec.cols)]
				#browser()
				temp.tab <- tbl(src.db, tab)
				#try it once to see what columns the evaluation brings back
				#temp.tab <- eval(parse(text=paste("select(temp.tab, ", tab.exp , ")")))
				temp.tab <- select_(temp.tab, .dots=as.list(unlist(strsplit(setNames(tab.exp, NULL), ","))))
				
				#if not all the columns necessary for joining are present, then add them and execute again
				if (all(nec.cols %in% colnames(temp.tab) ==F))
				{
				    diff.cols <- setdiff(nec.cols, colnames(temp.tab))
				    temp.tab <- tbl(src.db, tab)
				    #temp.tab <- eval(parse(text=paste("select(temp.tab, ", paste(diff.cols, collapse=",") , ",",tab.exp, ")")))
				    temp.tab <- select_(temp.tab, .dots=as.list(unlist(strsplit(setNames(c(diff.cols, tab.exp), NULL), ","))))
				    
				}
				
				return(temp.tab)
		    }
		    
		    #There can be tables needed simply to complete the query, not to retrieve columns from
		    if (use.path[1] %in% names(needed.tables))
		    {
				if (is.null(ancil.tables) == F && use.path[1] %in% names(ancil.tables))
				{
					
					#also make sure that all the necessary columns for joining to the ancilary tables are present
					all.tab <- .get.select.cols(use.path[1], needed.tables[use.path[1]], c(join.cols[[1]], unlist(ancil.join.cols[[use.path[1]]],use.names=F)), src.db)
					
					for(j in ancil.tables[[use.path[1]]])
					{
						if (j %in% names(needed.tables))
						{
							new.tab <- .get.select.cols(j, needed.tables[j], ancil.join.cols[[use.path[1]]][[j]], src.db)
						}else{
							new.tab <- tbl(src.db, j)
						}
						
						all.tab <- inner_join(all.tab, new.tab, by=ancil.join.cols[[use.path[1]]][[j]])
					}
					
				}else{
					all.tab <- .get.select.cols(use.path[1], needed.tables[use.path[1]], join.cols[[1]], src.db) 
				}
				
		    }else{
				all.tab <- tbl(src.db, use.path[1])
			
				if (is.null(ancil.tables) == F && use.path[1] %in% names(ancil.tables))
				{
					for(j in ancil.tables[[use.path[1]]])
					{
						if (j %in% names(needed.tables))
						{
							new.tab <- .get.select.cols(j, needed.tables[j], ancil.join.cols[[use.path[1]]][[j]], src.db)
						}else{
							new.tab <- tbl(src.db, j)
						}
						
						all.tab <- inner_join(all.tab, new.tab, by=ancil.join.cols[[use.path[1]]][[j]])
					}
				}
		    }
		    
		    use.path <- use.path[-1]
		    
		    for(i in seq_along(use.path))
		    {
				if (use.path[i] %in% names(needed.tables))
				{
				    if (is.null(ancil.tables) == F && use.path[i] %in% names(ancil.tables))
					{
						new.tab <- .get.select.cols(use.path[i], needed.tables[use.path[i]],
											c(join.cols[[i]],
											ifelse(i == length(use.path), NA, join.cols[[i+1]]),
											unlist(ancil.join.cols[[use.path[i]]],use.names=F)),
											src.db)
					
						for(j in ancil.tables[[use.path[i]]])
						{
							if (j %in% names(needed.tables))
							{
								new.ancil.tab <- .get.select.cols(j, needed.tables[j], ancil.join.cols[[use.path[i]]][[j]], src.db)
							}else{
								new.ancil.tab <- tbl(src.db, j)
							}
							
							new.tab <- inner_join(new.tab, new.ancil, by=ancil.join.cols[[use.path[i]]][[j]])
						}
					}else
					{
						new.tab <- .get.select.cols(use.path[i],
											   needed.tables[use.path[i]],
											   c(join.cols[[i]],
												ifelse(i == length(use.path), NA, join.cols[[i+1]])),
											   src.db) 
					}
				}else{
					
					new.tab <- tbl(src.db, use.path[i])
					
					if (is.null(ancil.tables) == F && use.path[i] %in% names(ancil.tables))
					{
						for(j in ancil.tables[[use.path[i]]])
						{
							if (j %in% names(needed.tables))
							{
								new.ancil.tab <- .get.select.cols(j, needed.tables[j], ancil.join.cols[[use.path[i]]][[j]], src.db)
							}else{
								new.ancil.tab <- tbl(src.db, j)
							}
							
							new.tab <- inner_join(new.tab, new.ancil.tab, by=ancil.join.cols[[use.path[i]]][[j]])
						}
					}
				}
			
				all.tab <- inner_join(all.tab, new.tab, by=join.cols[[i]])
		    }
		    
		    #make sure the final table only includes the requested columns
		    #all.tab <- eval(parse(text=paste("select(all.tab, ", paste(needed.tables, collapse=","), ")")))
		    #browser()
		    all.tab <- select_(all.tab, .dots=as.list(unlist(strsplit(setNames(needed.tables, NULL), ","))))
		}
		
	    }else{
		
		my_db <- src_sqlite(dbFile(obj), create = F)
		
		if (is.null(names(needed.tables))==F)
		{
		    all.tab <- tbl(my_db, names(needed.tables))
		    #all.tab <- eval(parse(text=paste("select(all.tab, ", needed.tables, ")")))
		    all.tab <- select_(all.tab, .dots=as.list(unlist(strsplit(setNames(needed.tables, NULL), ","))))
		    
		}else{
		    all.tab <- tbl(my_db, needed.tables)
		}
	    }
	    
	    return(all.tab)
	  
	  })


setGeneric("populate", def=function(obj, ...) standardGeneric("populate"))
setMethod("populate", signature("Database"), function(obj, ..., use.tables=NULL, should.debug=FALSE)
	  {
	    db.con <- dbConnect(SQLite(), dbFile(obj))
	    
	    .populate(schema(obj), db.con, ins.vals=list(...), use.tables=use.tables, should.debug=should.debug)
	    
	    invisible(dbDisconnect(db.con))
	  })

.populate <- function(obj, db.con,ins.vals=NULL, use.tables=NULL, should.debug=FALSE)
{
    db.schema <- obj
    
    if (class(db.con) != "SQLiteConnection")
    {
        stop("ERROR: db.con needs to be of class SQLiteConnection")
    }
    
    if (missing(ins.vals) || is.null(ins.vals))
    {
        stop("ERROR: ins.vals cannot be missing or NULL")
    }
	
    if (missing(use.tables) || is.null(use.tables) || is.na(use.tables))
    {
        use.tables <- schemaNames(db.schema)
    }
    else if (all(use.tables %in% schemaNames(db.schema)) == FALSE)
    {
        stop("ERROR: Invalid values for use.tables")
    }
    
    #schemaNames should be arranged in the order of population
    for(i in use.tables)
    {
        message(paste("Starting", i))
        #if table doesn't exist, then create it
        if (dbExistsTable(db.con, tableName(db.schema, i, mode="normal")) == FALSE)
        {
            if (should.debug) message("Creating database table")
            if (should.debug) message(createTable(db.schema, i, mode="normal"))
            dbGetQuery(db.con, createTable(db.schema, i, mode="normal"))
        }
        
        #then merge with existing databases as necessary

        if (shouldMerge(db.schema, i))
        {
            if (should.debug) message("Creating temporary table for merging")
            
            if (dbExistsTable(db.con, tableName(db.schema, i, mode="merge")))
            {
                stop("ERROR: Temporary tables should not exist prior to this loop")
            }
            
            if (should.debug) message(createTable(db.schema, i, mode="merge"))
            dbGetQuery(db.con, createTable(db.schema, i, mode="merge"))
            
            if (should.debug) message("Adding to temporary table")
            if (should.debug) message(insertStatement(db.schema, i, mode="merge"))
            #first add the data to temporary database
	   
	    dbBegin(db.con)
            dbGetPreparedQuery(db.con, insertStatement(db.schema, i, mode="merge"), bind.data = bindDataFunction(db.schema, i, ins.vals, mode="merge"))
            dbCommit(db.con)
            
            #merge from temporary into main table
            if (should.debug) message("Merging with existing table(s)")
            if (should.debug) message(mergeStatement(db.schema, i))
            dbGetQuery(db.con, mergeStatement(db.schema, i))
            
            #then also drop intermediate tables
            if (should.debug) message("Removing temporary table")
            if (should.debug) message(paste("DROP TABLE", tableName(db.schema, i, mode="merge")))
            dbGetQuery(db.con, paste("DROP TABLE", tableName(db.schema, i, mode="merge")))
        }else
        {
            if (should.debug) message("Adding to database table")
            if (should.debug) message(insertStatement(db.schema, i))
            #add the data to database
            dbBegin(db.con)
            dbGetPreparedQuery(db.con, insertStatement(db.schema, i), bind.data = bindDataFunction(db.schema, i, ins.vals))
            dbCommit(db.con)
        }
        
    }
}