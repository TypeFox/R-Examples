#
#    _____      _ ____  ____
#   / ___/_____(_) __ \/ __ )
#   \__ \/ ___/ / / / / __  |
#  ___/ / /__/ / /_/ / /_/ / 
# /____/\___/_/_____/_____/  
#
#
#
# BEGIN_COPYRIGHT
#
# This file is part of SciDB.
# Copyright (C) 2008-2014 SciDB, Inc.
#
# SciDB is free software: you can redistribute it and/or modify
# it under the terms of the AFFERO GNU General Public License as published by
# the Free Software Foundation.
#
# SciDB is distributed "AS-IS" AND WITHOUT ANY WARRANTY OF ANY KIND,
# INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
# NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE. See
# the AFFERO GNU General Public License for the complete license terms.
#
# You should have received a copy of the AFFERO GNU General Public License
# along with SciDB.  If not, see <http://www.gnu.org/licenses/agpl-3.0.html>
#
# END_COPYRIGHT
#

# The functions and methods defined below are based closely on native SciDB
# functions, some of which have weak or limited analogs in R. The functions
# defined below work with objects of class scidb (arrays), scidbdf (data
# frames). They can be efficiently nested by explicitly setting eval=FALSE on
# inner functions, deferring computation until eval=TRUE.

# SciDB rename wrapper
# Note! that the default garbage collection option here is to *not* remove.
rename = function(A, name=A@name, gc)
{
  if(!(inherits(A,"scidb") || inherits(A,"scidbdf"))) stop("`A` must be a scidb object.")
  if(missing(gc)) gc = FALSE
  query = sprintf("rename(%s,%s)",A@name, name)
  iquery(query)
  scidb(name,gc=gc)
}

remove_old_versions = function( stored_array )
{
  if(!(inherits(stored_array,"scidb") || inherits(stored_array,"scidbdf"))) stop("`stored_array` must be a scidb object.")
  versions_query = sprintf("aggregate(versions(%s), max(version_id) as max_version)", stored_array@name)
  versions = iqdf(versions_query, n=Inf)
  max_version = versions$max_version
  if (is.na(max_version))
  {
    stop("The array has no versions to remove")
  }
  remove_query = sprintf("remove_versions(%s, %i)", stored_array@name, max_version)
  iquery(remove_query)
}

unpack_scidb = function(x, `eval`=FALSE)
{
  dimname = make.unique_(c(dimensions(x),scidb_attributes(x)), "i")
  query = sprintf("unpack(%s, %s)",x@name, dimname)
  ans = .scidbeval(query,eval,depend=list(x))
  class(ans) = "scidbdf"
  ans
}

attribute_rename = function(x, old, `new`, `eval`=FALSE)
{
  atr = scidb_attributes(x)
  if(missing(old)) old=x@attributes
# Positional attributes
  if(is.numeric(old))
  {
    old = atr[old]
  }
  query = sprintf("attribute_rename(%s,%s)",x@name,
    paste(paste(old,new,sep=","),collapse=","))
  .scidbeval(query,eval,depend=list(x))
}

dimension_rename = function(x, old, `new`, `eval`=FALSE)
{
  if(!(is.scidb(x) || is.scidbdf(x))) stop("Requires a scidb or scidbdf object")
  if(missing(old)) old = dimensions(x)
  dnames = dimensions(x)
  if(!is.numeric(old))
  {
    old = which(dnames %in% old)
  }
  idx = old
  if(length(idx)!=length(new)) stop("Invalid old dimension name specified")
  dnames[idx] = `new`
  query = sprintf("cast(%s, %s%s)", x@name, build_attr_schema(x),
             build_dim_schema(x, newnames=dnames))
  .scidbeval(query,eval,depend=list(x))
}

slice = function(x, d, n, `eval`=FALSE)
{
  if(!(is.scidb(x) || is.scidbdf(x))) stop("Requires a scidb or scidbdf object")
  N = length(dimensions(x))
  i = d
  if(is.character(d))
  {
    i = which(dimensions(x) %in% d)
  }
  if(length(i)==0 || i>N)
  {
    stop("Invalid dimension specified")
  }
  if(missing(n)) n = scidb_coordinate_bounds(x)$start[i]
  query = sprintf("slice(%s, %s)",x@name,paste(paste(x@dimensions[i],noE(n),sep=","),collapse=","))
  .scidbeval(query,eval,depend=list(x))
}

# SciDB substitute wrapper. Default behavior strips nulls in a clever way.
replaceNA = function(x, value, `attribute`, `eval`=FALSE, ...)
{
  if(!(is.scidb(x) || is.scidbdf(x))) stop("Requires a scidb or scidbdf object")
  if(!any(scidb_nullable(x))) return(x)
  if(missing(attribute))
  {
    attribute = 1:length(scidb_attributes(x))
  }
  if(!is.numeric(attribute))
  {
    attribute = which(scidb_attributes(x) %in% attribute)
  }
  if(length(attribute)<1) stop("Invalid attribute(s)")
  if(missing(value))
  {
    ba = paste("build(",lapply(attribute,function(a) build_attr_schema(x,I=a,nullable=FALSE,newnames="____")),"[i=0:0,1,0]",",",.scidb_default_subst[scidb_types(x)],")",sep="")
  } else
  {
    ba = paste("build(",lapply(attribute,function(a) build_attr_schema(x,I=a,nullable=FALSE),newnames="____"),"[i=0:0,1,0]",",",value,")",sep="")
  }
  query = x@name
  for(j in 1:length(ba))
  {
    query = sprintf("substitute(%s,%s,%s)",query,ba[j],scidb_attributes(x)[j])
  }
  .scidbeval(query, `eval`, depend=list(x))
}

subarray = function(x, limits, between=FALSE, `eval`=FALSE)
{
  if(!(class(x) %in% c("scidb","scidbdf"))) stop("Invalid SciDB object")
  if(missing(limits)) limits=paste(rep("null",2*length(dimensions(x))),collapse=",")
  else if(is.character(limits))
  {
# Assume user has supplied a schema string
    limits = paste(between_coordinate_bounds(limits),collapse=",")
  } else if(is.scidb(limits) || is.scidbdf(limits))
  {
# User has supplied an array
    limits = paste(between_coordinate_bounds(schema(limits)),collapse=",")
  } else
  {
# Assume a vector of limits
    if(length(limits)!=2*length(dimensions(x))) stop("Mismatched bounds length")
    limits = paste(noE(limits),collapse=",")
  }
  limits = gsub("\\*",.scidb_DIM_MAX,limits)
  if(between)
    query = sprintf("between(%s,%s)",x@name,limits)
  else
    query = sprintf("subarray(%s,%s)",x@name,limits)
  .scidbeval(query, `eval`, depend=list(x))
}

cast = function(x, schema, `eval`=FALSE)
{
  if(!(class(x) %in% c("scidb","scidbdf"))) stop("Invalid SciDB object")
  if(missing(schema)) stop("Missing cast schema")
  if(is.scidb(schema) || is.scidbdf(schema)) schema = schema(schema) # wow!
  query = sprintf("cast(%s,%s)",x@name,schema)
  .scidbeval(query,eval,depend=list(x))
}

# SciDB build wrapper, intended to act something like the R 'array' function.
build = function(data, dim, names, type,
                 start, name, chunksize, overlap, gc=TRUE, `eval`=FALSE)
{
  if(missing(type))
  {
    type = typeof(data)
    if(is.character(data))
    {
      if(length(grep("\\(",data))>0) type="double"
      else
      {
        type = "string"
        data = sprintf("'%s'",data)
      }
    }
  }
# Special case:
  if(is.scidb(dim) || is.scidbdf(dim))
  {
    schema = sprintf("%s%s",build_attr_schema(dim,I=1),build_dim_schema(dim))
    query = sprintf("build(%s,%s)",schema,data)
    ans = .scidbeval(query,eval)
# We know that the output of build is not sparse
    attr(ans,"sparse") = FALSE
    return(ans)
  }
  if(missing(start)) start = rep(0,length(dim))
  if(missing(overlap)) overlap = rep(0,length(dim))
  if(missing(chunksize))
  {
    chunksize = rep(ceiling(1e6^(1/length(dim))),length(dim))
  }
  if(length(start)!=length(dim)) stop("Mismatched dimension/start lengths")
  if(length(chunksize)!=length(dim)) stop("Mismatched dimension/chunksize lengths")
  if(length(overlap)!=length(dim)) stop("Mismatched dimension/overlap lengths")
  if(missing(names))
  {
    names = c("val", letters[9:(8+length(dim))])
  }
# No scientific notation please
  chunksize = noE(chunksize)
  overlap = noE(overlap)
  dim = noE(dim + (start - 1))
  start = noE(start)
  schema = paste("<",names[1],":",type,">",sep="")
  schema = paste(schema, paste("[",paste(paste(paste(
        paste(names[-1],start,sep="="), dim, sep=":"),
        chunksize, overlap, sep=","), collapse=","),"]",sep=""), sep="")
  query = sprintf("build(%s,%s)",schema,data)
  if(missing(name)) return(.scidbeval(query,eval))
  ans = .scidbeval(query,eval,name)
# We know that the output of build is not sparse
  attr(ans,"sparse") = FALSE
  ans
}

# Count the number of non-empty cells
count = function(x)
{
  if(!(class(x) %in% c("scidb","scidbdf"))) stop("Invalid SciDB object")
  iquery(sprintf("aggregate(%s, count(*) as count)",x@name),return=TRUE)$count
}

# Filter the attributes of the scidb, scidbdf object to contain
# only those specified in expr.
# X:    a scidb, scidbdf object
# attributes: a character vector describing the list of attributes to project onto
# eval: a boolean value. If TRUE, the query is executed returning a scidb array.
#       If FALSE, a promise object describing the query is returned.
project = function(X,attributes,`eval`=FALSE)
{
  xname = X
  if(is.logical(attributes))
    attributes = X@attributes[which(attributes)]
  if(is.numeric(attributes))
    attributes = X@attributes[attributes]
  if(class(X) %in% c("scidbdf","scidb")) xname = X@name
  query = sprintf("project(%s,%s)", xname,paste(attributes,collapse=","))
  .scidbeval(query,eval,depend=list(X))
}

# This is the SciDB filter operation, not the R timeseries one.
# X is either a scidb, scidbdf object.
# expr is a valid SciDB expression (character)
# eval=TRUE means run the query and return a scidb object.
# eval=FALSE means return a promise object representing the query.
`filter_scidb` = function(X,expr,`eval`=FALSE)
{
  xname = X
  if(class(X) %in% c("scidbdf","scidb")) xname = X@name
# Check for special filter cases and adjust expr
  if(length(scidb_attributes(X))==2 && nchar(expr)==1)
  {
# Case 1: single comparison and two attributes
    expr = paste(scidb_attributes(X), collapse=expr)
  }
  query = sprintf("filter(%s,%s)", xname,expr)
  .scidbeval(query,eval,depend=list(X))
}


`index_lookup` = function(X, I, attr, new_attr, `eval`=FALSE)
{
  if(missing(attr)) attr = X@attributes[[1]]
  if(missing(new_attr)) new_attr=paste(attr,"index",sep="_")
  al = scidb_alias(X,I)
  xname = X
  if(class(X) %in% c("scidb","scidbdf")) xname=X@name
  iname = I
  if(class(I) %in% c("scidb","scidbdf"))
  {
    if(length(scidb_attributes(I))>1) I = project(I,1)
    if(scidb_nullable(I)) I = replaceNA(I)
    iname=I@name
  }
  query = sprintf("index_lookup(%s as %s, %s as %s, %s.%s, %s)",xname, al[1], iname, al[2], al[1], attr, new_attr)
  .scidbeval(query,eval,depend=list(X,I))
}

# Sort of like cbind for data frames.
bind = function(X, name, FUN, `eval`=FALSE)
{
  aname = X
  if(class(X) %in% c("scidb","scidbdf")) aname=X@name
# Auto-generate names like X_n:
  if(missing(name))
  {
    name = rep("X",length(FUN))
  }
  name = make.unique_(c(scidb_attributes(X),dimensions(X)), name)
  if(length(name)!=length(FUN)) stop("name and FUN must be character vectors of identical length")
  expr = paste(paste(name,FUN,sep=","),collapse=",")
  query = sprintf("apply(%s, %s)",aname, expr)
  .scidbeval(query,eval,depend=list(X))
}

unique_scidb = function(x, incomparables=FALSE, sort=TRUE, ...)
{
  mc = list(...)
  `eval` = ifelse(is.null(mc$eval), FALSE, mc$eval)
  if(incomparables!=FALSE) warning("The incomparables option is not available yet.")
  if(any(x@attributes %in% "i"))
  {
    new_attrs = x@attributes
    new_attrs = new_attrs[x@attributes %in% "i"] = make.unique_(x@attributes,"i")
    x = attribute_rename(x,x@attributes,new_attrs)
  }
  if(sort)
  {
# XXX XXX There is a problem here if there is an attribute named 'n' (see sort function
# below)...this must be fixed.
    rs = sprintf("%s[n=0:%s,%s,0]",build_attr_schema(x,I=1),.scidb_DIM_MAX,noE(min(1e6,prod(dim(x)))))
    if(length(x@attributes)>1)
    {
      query = sprintf("uniq(redimension(sort(project(%s,%s)),%s))",x@name,x@attributes[[1]],rs)
    }
    else
    {
      query = sprintf("uniq(redimension(sort(%s),%s))",x@name,rs)
    }
  } else
  {
    query = sprintf("uniq(%s)",x@name)
  }
  .scidbeval(query,eval,depend=list(x),`data.frame`=TRUE)
}

sort_scidb = function(X, decreasing = FALSE, ...)
{
  mc = list(...)
  if(!is.null(mc$na.last))
    warning("na.last option not supported by SciDB sort. Missing values are treated as less than other values by SciDB sort.")
  dflag = ifelse(decreasing, 'desc', 'asc')
# Check for ridiculous SciDB name conflict problem
  if(any(X@attributes %in% "n"))
  {
    new_attrs = X@attributes
    new_attrs = new_attrs[X@attributes %in% "n"] = make.unique_(X@attributes,"n")
    X = attribute_rename(X,X@attributes,new_attrs)
  }
  if(is.null(mc$attributes))
  {
    if(length(X@attributes)>1) warning("Array contains more than one attribute, sorting on all of them.\nUse the attributes= option to restrict the sort.")
    mc$attributes=X@attributes
  }
  `eval` = ifelse(is.null(mc$eval), FALSE, mc$eval)
  a = paste(paste(mc$attributes, dflag, sep=" "),collapse=",")
  if(!is.null(mc$chunk_size)) a = paste(a, mc$chunk_size, sep=",")

#  rs = sprintf("%s[n=0:%s,%s,0]",build_attr_schema(X,I=1),.scidb_DIM_MAX,noE(min(1e6,prod(dim(X)))))
#  query = sprintf("redimension(sort(%s,%s),%s)", X@name,a,rs)
  query = sprintf("sort(%s,%s)", X@name,a)
  .scidbeval(query,eval,depend=list(X))
}

# S3 methods
`merge.scidb` = function(x,y,by=intersect(dimensions(x),dimensions(y)),...) merge_scidb(x,y,by,...)
`merge.scidbdf` = function(x,y,by=intersect(dimensions(x),dimensions(y)),...) merge_scidb(x,y,by,...)
`sort.scidb` = function(...,na.last=TRUE,decreasing=FALSE) sort_scidb(...,decreasing)
`sort.scidbdf` = function(x,decreasing=FALSE,...) sort_scidb(x,decreasing,...)
`unique.scidb` = function(x,incomparables=FALSE,...) unique_scidb(x,incomparables,...)
`unique.scidbdf` = function(x,incomparables=FALSE,...) unique_scidb(x,incomparables,...)
`subset.scidb` = function(x,subset,...) filter_scidb(x,expr=subset,...)
`subset.scidbdf` = function(x,subset,...) filter_scidb(x,expr=subset,...)
