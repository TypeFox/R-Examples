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

# A utility function that returns TRUE if entries in the numeric vector j
# are sequential and in increasing order.
checkseq = function(j)
{
  if(!is.numeric(j)) return(FALSE)
  !any(diff(j)!=1)
}

# Returns a function that evaluates to a list of bounds
between = function(a,b)
{
  if(missing(b))
  {
    if(length(a)==2)
    {
      b = a[[2]]
      a = a[[1]]
    } else stop("between requires two arguments or a single argument with two elements")
  }
  l = lapply(list(a,b), function(x) gsub(".*Inf","null",noE(x)))
  function() l
}

`$.scidb` = function(x, ...)
{
  M = match.call()
  a = as.character(M[3])
  project(x, a)
}

# dimfilter: The workhorse array subsetting function
# INPUT
# x: A scidb object
# i: a list of index expressions
# eval: (logical) if TRUE, return a new SciDB array, otherwise just a promise
# drop: (logical) if TRUE, delete dimensions of an array with only one level
# redim: (logical) default TRUE, if false the supress bounding redimension
# OUTPUT
# a scidb object
#
# dimfilter distinguishes between four kinds of indexing operations:
# 'si' sequential numeric index range, for example c(1,2,3,4,5)
# 'bi' special between index range, that is a function or a list
#      that returns upper/lower limits
# 'ui' not specified range (everything, by R convention)
# 'ci' other, for example c(3,1,2,5) or c(1,1)
#
dimfilter = function(x, i, eval, drop, redim)
{
  if(missing(redim)) redim=TRUE
# Partition the indices into class:
# Identify sequential, numeric indices
  si = sapply(i, checkseq)
# Identify explicit between-type indices (functions, lists)
  bi = sapply(i, function(x) inherits(x,"function"))
# Unspecified range
  ui = sapply(i,is.null)
# Identify everything else
  ci = !(si | bi | ui)

  r = lapply(1:length(bi), function(j)
    {
      if(bi[j])
      {
# Just use the provided between-style range
        unlist(i[j][[1]]())
      }
      else if(si[j])
      {
# sequential numeric or lookup-type range
        noE(c(min(i[j][[1]]),max(i[j][[1]])))
      }
      else
       {
# Unspecified range or special index (ci case), which we handle later.
         c('null','null')
       }
    })
  ranges = r
  r = unlist(lapply(r,noE))
  everything = all(r %in% "null")
  ro = r[seq(from=1,to=length(r),by=2)]
  re = r[seq(from=2,to=length(r),by=2)]
  r = paste(c(ro,re),collapse=",")
  q = sprintf("between(%s,%s)",x@name,r)

  new_dimnames = c()
  new_depend = x
  if(!everything)
  {
    newstart = unlist(lapply(ranges,function(z)z[1]))
    newstart[newstart=="null"] = NA
    ina = is.na(as.numeric(newstart))
    if(any(ina))
    {
      newstart[ina] = scidb_coordinate_start(x)[ina]
    }
    newend = unlist(lapply(ranges,function(z)z[2]))
    newend[newend=="null"] = NA
    ina = is.na(newend)
    if(any(ina))
    {
      newend[ina] = scidb_coordinate_end(x)[ina]
    }
# Propagate dimension labels
    if(!is.null(dimnames(x)))
    {
      new_dimnames =  vector(mode="list",length=length(dim(x)))
      for(j in 1:length(new_dimnames))
      {
        if(!is.null(dimnames(x)[[j]]))
        {
          new_dimnames[[j]] = scidb(sprintf("project(join(redimension(%s,%s%s) as x, %s as y), y.%s)",q,build_attr_schema(x),build_dim_schema(x,I=j,newend="*"), dimnames(x)[[j]]@name, scidb_attributes(dimnames(x)[[j]])[[1]]))
# This dependency is already covered by x, which we already depend on:
#          new_depend = c(new_depend, dimnames(x)[[j]])
        }
      }
    }
# redim will be handled by the special_index function in the ci case:
    if(redim && !any(ci))
    {
      q = sprintf("subarray(redimension(%s, %s%s),%s)", q, build_attr_schema(x),
          build_dim_schema(x,newend=newend,newstart=newstart), paste(rep('null',2*length(dim(x))),collapse=","))
# We need to redimension any dimname arrays conformably. The only thing that
# could change is the starting coordinate.
      if(!is.null(new_dimnames))
      { for(j in 1:length(new_dimnames))
        { if(!is.null(new_dimnames[[j]]))
          { if(newstart[j] != scidb_coordinate_start(new_dimnames[[j]]))
            {
              new_dimnames[[j]] = scidb(sprintf("subarray(redimension(%s, %s%s),null,null)", new_dimnames[[j]]@name, build_attr_schema(new_dimnames[[j]]), build_dim_schema(new_dimnames[[j]],newstart=newstart[j],newend="*")))
            } else
            {
              new_dimnames[[j]] = scidb(sprintf("subarray(%s,null,null)", new_dimnames[[j]]@name ))
            }
          } } }
    }
  } else new_dimnames = dimnames(x)
# Return a new scidb array reference
  ans = .scidbeval(q,eval=FALSE,gc=TRUE,`data.frame`=FALSE,depend=new_depend)
#  dimnames(ans) = new_dimnames # no--it does too much housekeeping for us,
# we know that we already have a conformable schema so just update the dimname
# directly:
  ans@gc$dimnames = new_dimnames
  if(any(ci)) 
  {
    return(special_index(ans, ans@name, i, ci, eval, drop, redim))
  }
# Drop singleton dimensions if instructed to
  if(drop)
  {
    ans = drop_dim(ans)
  }
  if(`eval`)
  {
    ans = scidbeval(ans)
  }
  ans
}

# Helper function to drop dimensions of length 1 (cf R's drop)
drop_dim = function(ans)
{
  i = as.numeric(scidb_coordinate_bounds(ans)$length) == 1
  if(all(i))
  {
    i[1] = FALSE
  }
  if(any(i))
  {
    j = which(i)
    A = build_attr_schema(ans)
    D = build_dim_schema(ans,I=-j)
    query = sprintf("redimension(%s, %s%s)",ans@name,A,D)
    dn = dimnames(ans)
    ans = .scidbeval(query,`eval`=FALSE,depend=list(ans))
    if(!is.null(dn)) ans@gc$dimnames = dn[!i]
  }
  ans
}



# Materialize the single-attribute scidb array x as an R array.
materialize = function(x, drop=FALSE)
{

# If x has multiple attributes, warn.
  if(length(x@attributes)>1)
  {
    warnonce("unpack")
    return(iquery(x, return=TRUE,n=Inf))
  }
  type = names(.scidbtypes[.scidbtypes==scidb_types(x)])
# Check for types that are not fully supported yet.
  xstart = as.numeric(scidb_coordinate_start(x))
  attr = .get_attribute(x)
  if(length(type)<1)
  {
    u = unpack(x)[]
    ans = tryCatch(
      {
        array(dim=dim(x))
      },error = function(e)
      {
        n = length(dim(x))
        array(dim=apply(u[,1:n,drop=FALSE],2,function(x){max(x)+1}))
      })
    i = as.matrix(u[,1:length(dim(x))])
    for(j in 1:length(dim(x))) i[,j] = i[,j] + 1 - xstart[j]
    ans[i] = u[,ncol(u)]
    return(ans)
  }


# Set array index origin to zero. We need the zero origin here to reconstruct
# array indices in R.
  d     = dim(x)
  ndim  = length(dimensions(x))
  N     = paste(rep("null",2*ndim),collapse=",")
  query = sprintf("subarray(project(%s,%s),%s)",x@name,attr,N)

# Unpack into a staging data frame
  data  = scidb_unpack_to_dataframe(query)
  nelem = nrow(data)
  if(is.null(nelem)) nelem = 0
  p     = prod(d)
# Adjust indexing origin
  data[,1:ndim] = data[,1:ndim] + 1

# Handle coordinate labels
  labels = NULL
  if(!is.null(dimnames(x)))
  {
    labels = lapply(1:ndim, function(j)
    {
      if(is.null(dimnames(x)[[j]]))
        return(tryCatch(seq(from=as.numeric(scidb_coordinate_start(x)[j]),length.out=dim(x)[j]), error=invisible))
      if(!is.scidb(dimnames(x)[[j]])) return(NULL)
      dn = iquery(dimnames(x)[[j]]@name, `return`=TRUE, binary=TRUE, n=Inf)
      if(is.null(dn)) return("")
      if(nrow(dn)==0) return("")
      if(nrow(dn)==dim(x)[j])
      {
        nm = dn[,2]
      } else
      {
        dn[,1] = dn[,1]  - as.numeric(scidb_coordinate_start(x)[j]) + 1
        nm =  rep("",dim(x)[j])
        nm[dn[,1]] = dn[,2]
      }
      nm
    })
  } else
  {
    labels = lapply(1:ndim, function(j)
    {
      tryCatch(seq(from=as.numeric(scidb_coordinate_start(x)[j]),length.out=dim(x)[j]), error=invisible)
    })
  }

# Check for sparse matrix or sparse vector case. The tryCatch guards
# against unsupported types in R's sparse Matrix package and returns the
# raw data frames in bad cases.
  if(ndim==2 && nelem < p)
  {
    ans = tryCatch(
          {
            if(any(is.infinite(d)))
            {
              warnonce("toobig")
# de-adjust indexing origin
              data[,1:ndim] = data[,1:ndim] - 1
              return(data)
            }
            if(is.null(data)) t = Matrix::Matrix(0.0,nrow=dim(x)[1],ncol=dim(x)[2])
            else t = Matrix::sparseMatrix(i=data[,1],j=data[,2],x=data[,3],dims=d)
            dimnames(t) = labels
            t
          }, error=function(e) {warnonce("nonum"); data[,1:ndim] = data[,1:ndim] - 1; data})
    return(ans)
  } else if(ndim==1 && nelem < p)
  {
    ans = tryCatch(
          {
            if(any(is.infinite(dim(x))))
            {
              warnonce("toobig")
              data[,1:ndim] = data[,1:ndim] - 1
              return(data)
            }
            if(is.null(data)) t = Matrix::sparseVector(0.0,1,length=dim(x))
            else t = Matrix::sparseVector(i=data[,1],x=data[,2],length=p)
            t
          }, error=function(e) {warnonce("nonum");data[,1:ndim] = data[,1:ndim] - 1; data})
    return(ans)
  } else if(nelem < p)
  {
# Don't know how to represent this in R! (R only knows sparse vectors or arrays)
    warning("Note: R does not natively support sparse n-d objects for n>2. Returning data as a data frame.")
    data[,1:ndim] = data[,1:ndim] - 1
    return(data)
  }
# OK, we have a dense array of some kind
  if(length(d)==1) # a vector
  {
    ans = data[,2]
    names(ans) = labels[[1]]
    return(ans)
  }

  ans = array(NA, dim=d)  # A matrix or n-d array
  ans[as.matrix(data[,1:ndim])] = data[,ndim+1]
  dimnames(ans) = labels
  ans
}
