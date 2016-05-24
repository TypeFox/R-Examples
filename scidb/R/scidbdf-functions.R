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

# Element-wise operations
Ops.scidbdf = function(e1,e2) {
  switch(.Generic,
    '<' = .compare(e1,e2,"<"),
    '<=' = .compare(e1,e2,"<="),
    '>' = .compare(e1,e2,">"),
    '>=' = .compare(e1,e2,">="),
    '==' = .compare(e1,e2,"="),
    '!=' = .compare(e1,e2,"<>"),
    default = stop("Unsupported binary operation.")
  )
}

cbind.scidbdf = function(x, y)
{
  if(missing(y))
  {
    newdim=make.unique_(x@attributes, "j")
    nd = sprintf("%s[%s,%s=0:0,1,0]",build_attr_schema(x) , build_dim_schema(x,bracket=FALSE),newdim)
    return(redimension(bind(x,newdim,0), nd))
  }
  if(!is.scidb(y) && !is.scidbdf(y)) stop("cbind requires either a single argument or two SciDB arrays")
  i = intersect(dimensions(x),dimensions(y))
  if(length(i)<1) stop("Non-conformable arrays") # XXX Should really try harder
  merge(x,y,by=i)
}

rbind.scidbdf = function(x, y)
{
  c(x,y)
}

colnames.scidbdf = function(x)
{
  x@attributes
}

rownames.scidbdf = function(x)
{
  if(is.null(x@gc$dimnames)) return(NULL)
  x@gc$dimnames[[1]]
}

`rownames<-.scidbdf` = function(x, value)
{
  x
}

`dimnames<-.scidbdf` = function(x, value)
{
  y = do.call("names<-.scidbdf",args=list(x=x,value=value[[2]]))
  do.call("row.names<-.scidbdf",args=list(x=y,value=value[[1]])) # order matters here
}

row.names.scidbdf = function(x)
{
  if(is.null(x@gc$dimnames)) return(NULL)
  x@gc$dimnames[[1]]
}

`row.names<-.scidbdf` = function(x, value)
{
  y = x
  class(y) = "scidb"
  z = do.call("dimnames<-.scidb", args=list(x=y, value=list(value)))
  z@gc$depend=c(z@gc$depend, x)
  class(z) = "scidbdf"
  z
}

names.scidbdf = function(x)
{
  x@attributes
}

`names<-.scidbdf` = function(x,value)
{
  old = x@attributes
  if(is.null(value)) return(x)
  if(all(old==value)) return(x)
  if(length(value)!=length(old)) stop(paste("Incorrect number of names (should be",length(old),")"))
  arg = paste(paste(old,value,sep=","),collapse=",")
  query = sprintf("attribute_rename(%s,%s)",x@name,arg)
  ans = .scidbeval(query, eval=FALSE, gc=TRUE, depend=list(x))
  row.names(ans) = rownames(x)  # Preserve row names
  ans
}

dimnames.scidbdf = function(x)
{
  list(rownames.scidbdf(x), x@attributes)
}

`$.scidbdf` = function(x, ...)
{
  M = match.call()
  M[1] = call("[.scidbdf")
  M[2] = x
  M[4] = as.character(M[3])
  M[3] = expression(NULL)
  eval(M)
}

# data.frame subsetting wrapper.
# x: A Scidbdf array object
# ...: list of dimensions
# iterative: return a data.frame iterator
# n: if iterative, how many rows to return
# 
`[.scidbdf` = function(x, ..., iterative=FALSE, n=Inf, row.names)
{
  M = match.call()
  if(missing(row.names)) row.names=1
  else row.names = M[["row.names"]]
  drop = ifelse(is.null(M$drop),TRUE,M$drop)
  redim = ifelse(is.null(M$redim),TRUE,M$redim)
# Passing along a NULL argument is harder than it should be...
  M = M[3:length(M)]
  if(!is.null(names(M))) M = M[!(names(M) %in% c("drop","iterative","n","row.names","redim"))]
# i shall contain a list of requested index values
  E = parent.frame()
  i = lapply(1:length(M), function(j) tryCatch(eval(M[j][[1]],E),error=function(e)c()))
# User wants this materialized to R...
  if(all(sapply(i,is.null)))
    if(iterative)
    {
      ans = iquery(sprintf("%s",x@name),`return`=TRUE,iterative=TRUE,n=n,excludecol=1,colClasses=scidbdfcc(x))
      if(!is.null(dimnames(x)[[1]])) warning("row labels will not be displayed")
      return(ans)
    }
    else
    {
      if(!is.null(dimnames(x)[[1]]))
      {
        row.names = iquery(dimnames(x)[[1]]@name, `return`=TRUE,binary=TRUE)[,2]
        ans = iquery(sprintf("%s",x@name),`return`=TRUE,binary=TRUE,buffer=nrow(x),row.names=row.names)[,-1]
      } else
      {
        ans = iquery(sprintf("%s",x@name),`return`=TRUE,binary=TRUE,buffer=nrow(x),row.names=row.names)
      }
      return(ans)
    }
# Not materializing, return a SciDB array
  if(length(i)!=length(dim(x))) stop("Dimension mismatch")
  scidbdf_subset(x,i,drop,redim)
}

`dim.scidbdf` = function(x)
{
  c(as.numeric(scidb_coordinate_bounds(x)$length), length(scidb_attributes(x)))
}

`dim<-.scidbdf` = function(x, value)
{
  reshape(x,shape=value)
}


str.scidbdf = function(object, ...)
{
  .scidbstr(object)
}

ncol.scidbdf = function(x) length(scidb_attributes(x))
nrow.scidbdf = function(x) 
  {
    dim(x)[1]
  }
# This is consistent with regular data frames:
length.scidbdf = function(x) ncol(x)




# 'si' sequential numeric index range, for example c(1,2,3,4,5)
# 'bi' special between index range, that is a function that returns upper/lower limits
# 'ui' not specified range (everything, by R convention)
# 'ci' lookup-style range, a non-sequential numeric or labeled set, for example
#      c(3,3,1,5,3)   or  c('a1','a3')
scidbdf_subset = function(x, i, drop=FALSE, redim=TRUE)
{
  attribute_range = i[[2]]
  if(is.logical(attribute_range)) attribute_range = which(attribute_range)
  if(is.null(attribute_range)) attribute_range = x@attributes
  if(is.numeric(attribute_range))
  {
    attribute_range = x@attributes[attribute_range]
  }

  y = project(x, attribute_range)
  row.names(y) = rownames(x)  # Preserve row names
  class(y) = "scidb"
  ans = dimfilter(y, list(i[[1]]), `eval`=FALSE, drop=drop, redim=redim)
  if(!(length(dim(ans)==1) && drop)) class(ans) = "scidbdf"
  ans@gc$depend = c(ans@gc$depend, x)
  ans
}

betweenbound = function(x, m, n)
{
  ans = sprintf("between(%s, %s, %s)", x@name, noE(m), noE(n))
# Reset just the upper dimension index, use of redimension here is overkill
# but WE NEED IT HERE in case users over or undershoot dimension bounds.
  schema = sprintf("%s%s",build_attr_schema(x), build_dim_schema(x,newstart=m,newend=n))
  ans = sprintf("redimension(%s,%s)", ans, schema)
# seemingly more efficient but not general, sadly:
# XXX FIX ME...put an if/then clause here to detect efficient cases... Argggh.
#  s = subarray(x, c(m,n))
#  ans = sprintf("repart(subarray(%s, %s, %s),%s)", x@name, noE(m), noE(n),schema)
  ans
}
