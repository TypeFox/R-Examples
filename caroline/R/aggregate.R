 
 .pasteUnique <- function(x, ...){
    paste(unique(x), ...)
 }
 
 #=c('sum','mean','var','sd','max','min','length','concat','none')

.aggregateDFrow <- function(x, funcs, clmns, distinct=FALSE, collapse=',', ...){

  #funcs <- match.arg(funcs, several.ok=TRUE)
   
  if(length(clmns)!= length(funcs))
    stop('clmns length does not equal funcs length')
  dfl <- list()
  for(i in seq(along=clmns)){
    args <- list(x[,clmns[i]])
    if(funcs[i]=='paste'){
      args <-  append(list(collapse=collapse), args)
      if(distinct)
        funcs[i] <- '.pasteUnique'
    }
    if(!funcs[i] %in% c('length','paste'))
      args <- append(args, list(...))
    dfl[[paste(clmns[i],funcs[i],sep='_')]] <- do.call(funcs[i], args)
  }
  df <- as.data.frame(dfl, stringsAsFactors=FALSE);

  return(df)  
}


.rowlist2df <- function(x)
 do.call(rbind.data.frame, x)


bestBy <- function(df, by, best, clmns=names(df), inverse=FALSE, sql=FALSE){
  if(class(best) != 'character')
    stop('best column name must be of class caracter')
   
  if(!best %in% clmns){
    warning("'best' not included in 'clmns', adding it for you")
    clmns <- c(clmns, best)
  }
  byinclmns <- by %in% clmns
  if(!byinclmns)
      clmns <- c(clmns, by)
  
  if(!sql){    
    gb <- groupBy(df, by, clmns=clmns, aggregation=NULL)
    out <- .rowlist2df(lapply(gb, function(g)   g[order(g[,best], decreasing=inverse)[1],]))
    out <- out[order(out[,best], decreasing=inverse),clmns]			
  }else{
    require(RSQLite)
    
    dots.in.names <- any(grepl('\\.', names(df)))
    if(dots.in.names){
      ## protect SQL from R's allowance of '.'s in variable names
      by <- gsub('\\.','_', by)  
      names(df) <- gsub('\\.','_', names(df))  
      dotpos.clmns <- .gedots(clmns)     #remember
      clmns <- gsub('\\.','_', clmns)
    }
    
    tmpfile <- tempfile() 
    con <- dbConnect(dbDriver("SQLite"), dbname = tmpfile)
    dbWriteTable(con, 'tab', df)
    sql <- paste("SELECT * FROM tab AS tab1 WHERE tab1.oid IN 
                                   (SELECT tab2.oid FROM tab AS tab2
         				WHERE tab1.",by," = tab2.",by,"
         				ORDER BY tab2.", best," ", c('ASC','DESC')[inverse+1],"
     				LIMIT 1)", sep='')
    out <- dbGetQuery(con, sql)
    out <- out[order(out[,best], decreasing = inverse), clmns]
    
    ## convert underscores back to '.'s 
    if(dots.in.names)
      names(out) <- .redots(names(out), dotpos.clmns) 
    
    rownames(out)<- out[,by]
    if(!byinclmns)
	out[,by] <- NULL
  }
  return(out)
}


## the next two functions compensate for sql turning '.'s in names to '_' 
## it finds the positions, let's R -> SQL do it's conversion ...
.gedots <- function(v)
  lapply(v,function(x) gregexpr('\\.',x)[[1]])
  
# ... and then remembers where the original '.'s were and puts them back
.redots <- function(v, pos.list){
  for(i in seq(along=v)){
    if(any(pos.list[[i]] < 0)){
      next
    }else{
     for(pos in pos.list[[i]]){
       begin <- substr(v[i],1,pos-1)
       end <- substr(v[i],pos+1, nchar(v[i]))
       v[i] <-  paste(begin,end, sep='.')
     }
    } 
   }  
  return(v)
}

groupBy <- function(df, by, aggregation, clmns=names(df), collapse=',', distinct=FALSE, sql=FALSE, full.names=FALSE, ...){

  #aggregation <- match.arg(aggregation, several.ok=TRUE)
  
  #if(any(aggregation!='none'))
  if(!any(is.null(aggregation)))
    if(length(aggregation) != length(clmns))
      if(length(aggregation) == 1){
        warning(paste("automatically extending your only specified aggregation to all ",length(clmns)," clmns"))
        aggregation <- rep(aggregation, length(clmns))
      }else{
        stop("length of 'aggregation' does not equal length of 'clmns'")
      }
  
  if(is.character(by)){
    if(!by %in% names(df)){
      stop("character argument 'by' not in names of dataframe 'df'")
    }else{
       if(!sql)
         by <- df[,by]
    }
  }else{  #not character
    if(sql)
      stop('sql version requires character field name for grouping factor')
  }
  
  if(any(sapply(as.character(by), nchar)> 255))
    stop("levels of by must have no more than 256 characters")

  if(!sql){

    if(any(is.null(aggregation))){
      return(by(df, by, function(x) x[,clmns]))
    }else{
      #aggregation <- sub('count','length', aggregation)
      out <- .rowlist2df(
               by(df, by, function(x) .aggregateDFrow(x, funcs=aggregation, clmns=clmns, collapse=collapse, ...)))
    }
  }else{
    if(any(is.null(aggregation)))
      stop('cannot return group by with out agregation for sql version')
    if(any(c('sd','var') %in% aggregation))
      stop('SQLite does not yet support STDDEV or VARIANCE SQL calls')
   
    supported.sql.funcs <- c('sum','mean','max','min','length','paste')
    if(!all(aggregation %in% supported.sql.funcs))
      stop("In SQL mode, please stick to the supported functions: ", paste(supported.sql.funcs, collapse=','))
    
    dots.in.names <- any(grepl('\\.', names(df)))
    if(dots.in.names){
        ## protect SQL from R's allowance of '.'s in variable names
	by <- gsub('\\.','_', by)  
	names(df) <- gsub('\\.','_', names(df))  
        dotpos.clmns <- .gedots(clmns) #remember positions for reinitroduction
        clmns <- gsub('\\.','_', clmns)
    }
    
    require(RSQLite)
    tmpfile <- tempfile() 
    con <- dbConnect(dbDriver("SQLite"), dbname = tmpfile)
    ## sort by 'by' column: important for ordering match to 'unique' row renaming below
    df <- df[order(df[,by]),]
    dbWriteTable(con, 'tab', df)
    aggregation <- sub('mean','avg', aggregation)
    aggregation <- sub('length','count', aggregation)
    #aggregation <- sub('sd','stddev', aggregation)  #not yet supported by SQLite
    #aggregation <- sub('var','variance', aggregation) #not yet supported by SQLite
    sql.concat <- 'group_concat('
    sql.collapse <- paste(", '",collapse,"'",sep='')

    if(distinct){
      sql.concat <- paste(sql.concat, 'DISTINCT ', sep='')
      if(collapse == ','){
        sql.collapse <- ''
      }else{
        stop("collapse must be ',' when distinct, sql are both TRUE")
      }
    }
	
    select <- paste(paste(aggregation,"(",clmns,") AS ", clmns,'_',aggregation, sep=''), collapse=', ')  #
    select <- gsub('paste\\(([^\\)]+)',paste(sql.concat, "\\1", sql.collapse ,sep=''), select) #array_to_string(array_agg(field), '; ') #postgresql
    
    ## add the 'by' variable
    select <- paste(by,"AS sortbyid,", select)
    
    sql <- paste("SELECT", select, "FROM tab GROUP BY", by)
    out <- dbGetQuery(con, sql)
    
    rownames(out) <- out$sortbyid
    out$sortbyid <- NULL
    
    ## convert underscores back to '.'s 
    if(dots.in.names){
      clmns      <- .redots(clmns     , dotpos.clmns)  
      names(out) <- .redots(names(out), dotpos.clmns) 
    }
  }

  if(!full.names){
    if(length(unique(clmns)) != ncol(out))
      warning('Selected columns are not unique, cannot use original names. Using full names instead')
    else
      names(out) <- clmns
  }
  
  
  return(out)
}


#.dbGroupBy <- function(df, by, clmns, aggregation){

#}


regroup <- function(df, old, new, clmns=names(df), funcs=rep('sum',length(clmns)), combine=TRUE){
  
  if(length(old) != length(new))
    stop('old and new must be the same length')

  if(combine)
    groupings <- row.names(df)
  else
    groupings <- rep(NA, nrow(df))
  matches <- match(old, row.names(df))
  groupings[matches[!is.na(matches)]] <- new[!is.na(matches)]

  regroup <- groupBy(df, by=as.factor(groupings), funcs, clmns=clmns)
  return(regroup)
}

  
addFactLevs <- function(x, new.levs=NA){

  if(class(x) != 'data.frame')
    stop('x must be a data.frame')

  isfacts <- sapply(x, is.factor)
  for(clmn in names(x)[isfacts])
    x[,clmn] <- factor(x[,clmn], levels=c(levels(x[,clmn]), new.levs))
  x
}



sstable <- function(x, idx.clmns, ct.clmns=NULL, na.label='NA'){#,exclude=exclude, ...){

  if(class(x) != 'data.frame')
    stop('x must be a data.frame')

  if(any(sapply(x[,idx.clmns], class)!= 'factor')){
    warning('some or all of your index columns are not factors. coercing them...')
    for(idx in idx.clmns)
      x[,idx] <- as.factor(x[,idx])
  }
    
  
  # should I make this run more like the other 'table' functions?
  if(is.null(ct.clmns)){
    ## plain old table method
    tab <- table(x[,idx.clmns], exclude='ifany')
    rownames(tab)[is.na(rownames(tab))] <- na.label
    
  }else{
    ## this whole section need reworking into a much simpler agg2tab function perhaps    
    ncc <- length(ct.clmns)
    nic <- length(idx.clmns)
  
    if(nic > 2)
      stop("can't have more than 2 index columns")
    if(ncc > 1 & nic > 1)
      stop('cannot have both more than one index and count columns')
       
    y <- tab2df(x[,idx.clmns])
    y <- addFactLevs(y, 'NA')
    y[is.na(y)] <- na.label
      
    if(nic == 2 & ncc == 1){
      tmp <- as.table( by(x[,ct.clmns], as.list(y), sum))
    }else{
      tmp <- as.table(sapply(ct.clmns, function(cc)
                             by(x[,cc], as.list(y), sum)
                             )
                      )
    }

    tab <- tab2df(tmp)
    tab [is.na(tab)] <- 0

    if(dim(tmp)[1]==2 & NROW(tab)==1 & ncc==2){
      colnames(tab) <- ct.clmns
      rownames(tab) <- unique(y)
    }
    
  }
  if(NCOL(tab) == 1 & NROW(tab) == 1){
    ret <- data.frame(tab[1,1], row.names=rownames(tab))
    colnames(ret) <- colnames(tab)
    return(ret)
  }
  
  if(NCOL(tab) > 1){
    tab <- tab2df(tab)    
    tab$sum <- apply(tab, 1, sum)

    tab <- subset(tab, sum > 0)
    
    return(tab[rev(order(tab$sum)),])
  }else{
    return(tab[rev(order(tab[,1])),])
  }
}


leghead <- function(x, n=7, tabulate=FALSE, colors=TRUE, na.name='NA',na.col='white', other.col='gray', na.last=TRUE){

  if(as.logical(tabulate))
    x <- sstable(x, tabulate)
  
  if(n > 3 & n < nrow(x)){
    
    o.rn <- rownames(x)    
    is.n <- sapply(x, is.numeric)
    n <- n + sum(!is.n)  # so that we can utilize (an) extra color(s)
    others.summed <- apply(x[(n+1):nrow(x), is.n],2,sum)    
    x <- rbind(head(x,n), c(others.summed, 'NA'[!is.n]))
    rownames(x) <- c(o.rn[1:n],'other')
    
  }

  if('color' %in% colnames(x)){
    warn.col.clmn <- 'dataframe already has a color column:'
    isna <- rownames(x)== 'NA'
    if(sum(isna) > 0){
      warn.col.clmn <- paste(warn.col.clmn , 'changing NA levels to', na.col, '.')
      x$color[isna] <- na.col
    }else{
      warn.col.clmn <- paste(warn.col.clmn , 'leaving it as is.')
    }
    warning(warn.col.clmn)
    colors <- FALSE
  }
  
  if(all(colors != FALSE)){
    
    if(length(colors)>1){

      matches <- match(names(colors),rownames(x))
      if(all(is.na(matches)))
        stop('color vector names do not match the top n table row names')
      x$color <- colors[matches]
      
    }else{
      isna <- rownames(x)=='NA'
      x$color[!isna] <- rainbow(nrow(x) - sum(isna))
      x$color[ isna] <- na.col
    }
  }
  if('color' %in% colnames(x))
    x$color[rownames(x)=='other'] <- other.col  

  if(na.name != 'NA')
    rownames(x)[match('NA', rownames(x))] <- na.name
    
  if(na.last){
    isna <- rownames(x)==na.name
    if(sum(isna) > 0)
      x <- x[c(rownames(x)[!isna], na.name),]
  }
  x
}




.vle2df <- function(vl,i){
  ## vector list element to dataframe (preserves list element names)

  if(class(vl[[i]])=='data.frame'){
    df <- vl[[i]]
    names(df) <- paste(names(df),'.',i,sep='')    
  }else{
    df <- as.data.frame(vl[[i]])
    colnames(df) <- i
    rownames(df) <- names(vl[[i]])
  }
  df$rownames <- rownames(df)  #necessary because the built in b='row.names' merge is really slow (if not completely broken)
  df
}


nerge <- function(l, ...){
  ## named data.frame or vector merge

  if(!all(sapply(l, function(k) class(k) == 'data.frame' | is.vector(k) | is.factor(k))))
     stop('list elements must be either of class data.frame or of type vector (or factor)')
  if(length(l) < 2)
    stop('list l must have at least 2 elements')
  if(is.null(names(l))){
    warning("each merge element in the list 'l' should have a name. making some up.")
    names(l) <- letters[1:length(l)]
  }
  if(!all(sapply(l, function(j) !is.null(rownames(j)) | !is.null(names(j)))))
    stop('all list elements must have named components (dataframes must have row names)')
    
  df <- .vle2df(l,names(l)[1])
  for(e in 2:length(l)){
    df <- merge(df, .vle2df(l, names(l)[e]), by='rownames', ...)
    rownames(df) <- df$rownames
  }
  df$rownames <- NULL


  ## removing appended colnames if unnecessary
  orig.names <- sub(paste('\\.(',paste(names(l),collapse='|'),')$',sep=''),'', names(df))
  if(length(orig.names)== length(unique(orig.names)))
    names(df) <- orig.names
  
  df
}


