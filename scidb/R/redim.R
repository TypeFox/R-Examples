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

# Array redimension, repart, reshape operations.
reshape_scidb = function(data, schema, shape, dimnames, start, chunks, `eval`=FALSE)
{
  if(!missing(schema))
  {
    if(is.scidb(schema)||is.scidbdf(schema)) schema=schema(schema) # <- that's nutty notation Malkovich!
    query = sprintf("reshape(%s,%s)",data@name,schema)
    return(.scidbeval(query,eval,depend=list(data)))
  }
  if(missing(shape)) stop("Missing dimension shape")
  N = length(shape)
  if(missing(dimnames))
  {
    dimnames=letters[9:(9+N-1)]
  }
  if(missing(chunks))
  {
    chunks = ceiling(1e6^(1/N))
  }
  if(missing(start)) start = rep(0,N)
  shape = shape - 1 + start
  D = build_dim_schema(data, newstart=start, newnames=dimnames, newend=shape, newchunk=chunks)
  query = sprintf("reshape(%s,%s%s)",data@name,build_attr_schema(data),D)
  .scidbeval(query,eval,depend=list(data))
}

repart = function(x, schema, upper, chunk, overlap, `eval`=FALSE)
{
  if(!missing(schema))
  {
    query = sprintf("repart(%s, %s)", x@name, schema)
    return(.scidbeval(query,eval,depend=list(x)))
  }
  if(missing(upper)) upper = scidb_coordinate_end(x)
  if(missing(chunk)) chunk = scidb_coordinate_chunksize(x)
  if(missing(overlap)) overlap = scidb_coordinate_overlap(x)
  a = build_attr_schema(x)
  schema = sprintf("%s%s", a, build_dim_schema(x,newend=upper,newchunk=chunk,newoverlap=overlap))
  query = sprintf("repart(%s, %s)", x@name, schema)
  .scidbeval(query,eval,depend=list(x))
}

# SciDB redimension wrapper
#
# Either supply schema or dim. dim is a list of new dimensions made up from the
# attributes and existing dimensions. FUN is a limited scidb aggregation
# expression.
redimension = function(x, schema, dim, FUN)
{
  if(!(class(x) %in% c("scidb","scidbdf"))) stop("Invalid SciDB object")
# NB SciDB NULL is not allowed along a coordinate axis prior to SciDB 12.11,
# which could lead to a run time error here.
  if(missing(schema)) schema = NULL
  if(missing(dim)) dim = NULL
  s = schema
  if(is.null(s) && is.null(dim) ||
    (!is.null(s) && !is.null(dim)))
  {
    stop("Exactly one of schema or dim must be specified")
  }
  if((class(s) %in% c("scidb","scidbdf"))) s = schema(s)
  dnames = c()
  if(!is.null(dim))
  {
    d = unlist(dim)
    dnames = vector("list",length(dim))
    ia = which(scidb_attributes(x) %in% d)
    if(is.numeric(d)) id = d
    else id = which(dimensions(x) %in% d)
    if(length(ia)<1 && length(id)<1) stop("Invalid dimensions")
    as = build_attr_schema(x, I=-ia)
    if(length(id>0))
    {
      ds = build_dim_schema(x, I=id, bracket=FALSE)
    } else
    {
      ds = c()
    }
    if(length(ia)>0)
    {
# We'll be converting attributes to dimensions here.
# First, we make sure that they are all int64. If not, we add a new
# auxiliary attribute with index_lookup and dimension along that instead.
      reindexed = FALSE
      xold = x
      for(nid in x@attributes[ia])
      {
        idx = which(x@attributes %in% nid)
        if(scidb_types(x)[idx] != "int64")
        {
          reindexed = TRUE
          didx = which(d %in% nid)
          newat = sprintf("%s_index",nid)
          newat = make.unique_(x@attributes, newat)
          dnames[didx] = unique(xold[,nid])
          names(dnames)[didx] = newat
          x = index_lookup(x, unique(xold[,nid]), nid, newat)
          d[didx] = newat
        }
      }
      if(reindexed)
      {
        ia = which(x@attributes %in% d)
        as = build_attr_schema(x, I=-ia) # remove _index attributes
      }
# Note! we need to keep around the original nid attributes for the
# index_lookups. This can screw with some aggregation functions.
# In such cases, use aggregate for now. Eventually, maybe split
# the nids out into a new array?

# Add the new dimension(s)
      a = x@attributes[ia]
      x@attributes = x@attributes[-ia]
      f = paste(paste("min(",a,"), max(",a,")",sep=""),collapse=",")
# Explicitly bounding this dimension is nice, but not necessary. Can
# just use a '*' bound instead (cheaper)--see commented line... XXX
      m = matrix(aggregate(x, FUN=f, unpack=FALSE)[],ncol=2,byrow=TRUE)
#      m = cbind(rep("0",length(a)), rep("*",length(a)))
      p = prod(as.numeric(scidb_coordinate_chunksize(x)[-id]))
      chunk = ceiling((1e6/p)^(1/length(ia)))
      new = apply(m,1,paste,collapse=":")
      new = paste(a,new,sep="=")
      new = paste(new, noE(chunk), "0", sep=",")
      new = paste(new,collapse=",")
      ds = ifelse(length(ds)>0,paste(ds,new,sep=","),new)
    }
    s = sprintf("%s[%s]",as,ds)
  }
# Check to see if the new schema is identical to the original schema.
# If so, don't bother with redimension, and return the input
  if(isTRUE(compare_schema(x,s)))
  {
    return(x)
  }
  if(!missing(FUN))
  {
    fn = NULL
    if(is.function(FUN))
    {
# convert to aggregation expression that can be parsed by aparser
# special count case (usually what is desired)
      if(isTRUE(all.equal(FUN,count)) || isTRUE(all.equal(FUN,length)))
        fn = "count(*) as count"
      else
        fn = paste(sprintf("%s(%s) as %s",.scidbfun(FUN),scidb_attributes(s),scidb_attributes(s)), collapse=",")
    }
    if(is.character(FUN)) fn = FUN
    if(is.null(fn))
      stop("`FUN` requires a function or SciDB aggregation expression")
    if(is.null(schema))
    {
      s = sprintf("%s%s",aparser(x,fn), build_dim_schema(s))
    }
    s = sprintf("%s, %s", s, fn)
  }
  query = sprintf("redimension(%s,%s)",x@name,s)
  ans = .scidbeval(query,`eval`=FALSE,depend=list(x))
  if(!is.null(dnames))
  {
    if(class(ans) %in% "scidbdf") rownames(ans)=dnames[[1]]
    else {
      if(any(dimensions(ans) %in% names(dnames)))
      {
        dnew = vector("list",length(dimensions(ans)))
        names(dnew) = dimensions(ans)
        dnames = dnames[unlist(lapply(dnames,function(x)!is.null(x)))]
        dnew[names(dnames)] = dnames
        dimnames(ans) = dnew
      }
    }
  }
  ans
}
