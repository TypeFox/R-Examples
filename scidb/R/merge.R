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

# Auxillary merge functions for each special case follow. The main function
# appears athe bottom of this file.

# Special case 1: full set cross product
merge_scidb_cross = function(x,y)
{
# New attribute schema for y that won't conflict with x:
  newas = build_attr_schema(y,newnames=make.unique_(x@attributes,y@attributes))
# Impose a reasonable chunk size for dense arrays
  chunky = scidb_coordinate_chunksize(y)
  chunkx   = scidb_coordinate_chunksize(x)
  chunk_elements = prod(c(as.numeric(chunky),as.numeric(chunkx)))
# Only compute these counts if we need to
  pdx = prod(dim(x))
  if(is.scidbdf(x)) pdx = dim(x)[1]
  pdy = prod(dim(y))
  if(is.scidbdf(y)) pdy = dim(y)[1]
  if(chunk_elements>1e6 && pdx==count(x) && pdy==count(y))
    {
      NC = length(chunkx) + length(chunky)
      NS = 1e6^(1/NC)
      chunky = rep(noE(NS), length(chunky))
      chunkx = rep(noE(NS), length(chunkx))
      x = repart(x,sprintf("%s%s",build_attr_schema(x), build_dim_schema(x,newchunk=chunkx)))
      y = repart(y,sprintf("%s%s",build_attr_schema(y), build_dim_schema(y,newchunk=chunkx)))
    }
  newds = build_dim_schema(y,newnames=make.unique_(x@dimensions,y@dimensions))
  y = cast(y,sprintf("%s%s",newas,newds))
  query = sprintf("cross_join(%s, %s)",x@name,y@name)
  return(.scidbeval(query,FALSE,depend=list(x,y)))
}

merge_scidb_on_attributes = function(x,y,by.x,by.y)
{
  `eval`=FALSE
  by.x = by.x[[1]]  # Limitation: only one attribute for now
  by.y = by.y[[1]]  # Ditto
  al = scidb_alias(x,y)  # new alias names
  lkup = unique(project(x,by.x),attributes=by.x)
  XI = index_lookup(x,lkup,by.x,`eval`=FALSE)
  YI = index_lookup(y,lkup,by.y,`eval`=FALSE)

  new_dim_name = make.unique_(c(dimensions(x),dimensions(y)),"row")
  a = XI@attributes %in% paste(by.x,"index",sep="_")
  n = XI@attributes[a]
  redim = paste(paste(n,"=-1:",.scidb_DIM_MAX,",100000,0",sep=""), collapse=",")
  S = build_attr_schema(x, I=!(x@attributes %in% by.x))
  D = sprintf("[%s,%s]",redim,build_dim_schema(x,bracket=FALSE))
  q1 = sprintf("redimension(substitute(%s,build(<_i_:int64>[_j_=0:0,1,0],-1),%s),%s%s)",XI@name,n,S,D)

  a = YI@attributes %in% paste(by.y,"index",sep="_")
  n = YI@attributes[a]
  redim = paste(paste(n,"=-1:",.scidb_DIM_MAX,",100000,0",sep=""), collapse=",")
  S = build_attr_schema(y)
  D = sprintf("[%s,%s]",redim,build_dim_schema(y,bracket=FALSE))
  D2 = sprintf("[%s,_%s]",redim,build_dim_schema(y,bracket=FALSE))
  q2 = sprintf("cast(redimension(substitute(%s,build(<_i_:int64>[_j_=0:0,1,0],-1),%s),%s%s),%s%s)",YI@name,n,S,D,S,D2)
  query = sprintf("unpack(cross_join(%s as %s, %s as %s, %s.%s_index, %s.%s_index),%s)",q1,al[1],q2,al[2],al[1],by.x,al[2],by.y,new_dim_name)
  return(.scidbeval(query,eval,depend=list(x,y)))
}

# SciDB join, cross_join, and merge wrapper internal function to support merge
# on various classes (scidb, scidbdf). This is an internal function to support
# R's merge on various SciDB objects.
#
# x and y are SciDB array references of any kind (scidb, scidbdf)
# `by` is either a single character indicating a dimension name common to both
#      arrays to join on, or a two-element list of character vectors of array
#      dimensions to join on.
# `fillin` is an optional argument specifying a value used to fill attributes
#          as required by merge, it defaults to null.
# `all` is an optional argument that, if TRUE, indicates outer join. It only
#       applies in limited settings. The default is inner join.
#
`merge_scidb` = function(x,y,`by`,...)
{
  mc = list(...)
  al = scidb_alias(x,y)
  by.x = by.y = NULL
  `all` = FALSE
  scidbmerge = FALSE
  fillin = "(null)"
  if(!is.null(mc$all)) `all` = mc$all
  if(!is.null(mc$by.x)) by.x = mc$by.x
  if(!is.null(mc$by.y)) by.y = mc$by.y
  if(!is.null(mc$merge)) scidbmerge = mc$merge
  if(!is.null(mc$fillin)) fillin = sprintf("(%s)",mc$fillin)
  `eval` = FALSE
  xname = x@name
  yname = y@name

# Check input
  if(sum(!is.null(by.x), !is.null(by.y))==1)
  {
    stop("Either both or none of by.x and by.y must be specified.")
  }
  if((!is.null(by.x) && !is.null(by.y)))
  {
    `by` = NULL
  }

# Check for full cross case.
  if((is.null(`by`) && is.null(by.x) && is.null(by.y)) ||
      length(`by`)==0 && is.null(by.x) && is.null(by.y))
  {
    if(scidbmerge) stop("SciDB merge not supported in this context")
    return(merge_scidb_cross(x,y))
  }

# Convert identically specified by into separate by.x by.y
  if(length(by)>0)
  {
    by.x = `by`
    by.y = `by`
  }

# Check for numeric `by` specification (dimension index)
  if(is.numeric(by.x)) by.x = dimensions(x)[by.x]
  if(is.numeric(by.y)) by.y = dimensions(y)[by.y]

# Check for special join on attributes case (limited applicability)
# In particular:
# - join on only one attribute per array
# - only inner join
  if(all(by.x %in% x@attributes) && all(by.y %in% y@attributes))
  {
    if(scidbmerge) stop("SciDB merge not supported in this context")
    return(merge_scidb_on_attributes(x,y,by.x,by.y))
  }

# OK, we've ruled out cross and attribute join special cases. We have left
# either the normal SciDB join/merge or cross_join on a subset of dimensions.

# New attribute schema for y that won't conflict with x:
  newas = build_attr_schema(y,newnames=make.unique_(x@attributes,y@attributes))
# Check for join case (easy case)
  if((length(by.x) == length(by.y)) && all(dimensions(x) %in% by.x) && all(dimensions(y) %in% by.y))
  {
# Check for valid starting coordinates (they must be identical)
    if(!isTRUE(all.equal(scidb_coordinate_start(x),scidb_coordinate_start(y))))
    {
#      stop("Mis-matched starting coordinates") # used to error out, now try to redim
# try inserting a redimension
       xless = scidb_coordinate_start(x) < scidb_coordinate_start(y)
       yless = scidb_coordinate_start(y) < scidb_coordinate_start(x)
       if(any(xless))
       {
         newstart = scidb_coordinate_start(y)
         newstart[xless] = scidb_coordinate_start(x)
         y = redimension(y, schema=sprintf("%s%s", build_attr_schema(y), build_dim_schema(y,newstart=newstart)))
       }
       if(any(yless))
       {
         newstart = scidb_coordinate_start(x)
         newstart[yless] = scidb_coordinate_start(y)
         x = redimension(x, schema=sprintf("%s%s", build_attr_schema(x), build_dim_schema(x,newstart=newstart)))
       }
    }
# If the chunk sizes are identical, we're OK (join does not care about the
# upper array bounds). Otherwise we need redimension.
    if(!isTRUE(all.equal(scidb_coordinate_chunksize(x), scidb_coordinate_chunksize(y))))
    {
      newds = build_dim_schema(y,newnames=dimensions(x))
      castschema = sprintf("%s%s", newas, newds)
      reschema = sprintf("%s%s", newas,build_dim_schema(x))
      z = redimension(cast(y,castschema),reschema)
    } else
    {
      castschema = sprintf("%s%s", newas, build_dim_schema(y))
      z = cast(y, castschema)
    }
    if(all)
    {
# Experimental outer join XXX XXX
      if(scidbmerge) stop("at most one of `all` and `merge` may be set TRUE")
      x = make_nullable(x)
      z = make_nullable(z)
# Form a null-valued version of each array in the alternate array coordinate system
      xnames = make.unique_(c(dimensions(z),scidb_attributes(z)),scidb_attributes(x))
      vals = paste(scidb_types(x), rep(fillin,length(scidb_types(x))))
      xnull = make_nullable(attribute_rename(project(bind(z,xnames,vals),xnames),xnames,scidb_attributes(x)))
      znames = make.unique_(c(dimensions(x),scidb_attributes(x)),scidb_attributes(z))
      vals = paste(scidb_types(z), rep(fillin,length(scidb_types(z))))
      znull = make_nullable(attribute_rename(project(bind(x,znames,vals),znames),znames,scidb_attributes(z)))
# Merge each array with its nullified counterpart, then join:
      query = sprintf("join(merge(%s,%s),merge(%s,%s))",x@name,xnull@name,z@name,znull@name)
    }
    else
      if(scidbmerge)
      {
        query = sprintf("merge(%s,%s)",x@name,z@name)
      } else
      {
        query = sprintf("join(%s,%s)",x@name,z@name)
      }
    return(.scidbeval(query,eval,depend=list(x,y)))
  }


# Finally, the cross-join case (trickiest)
  if(scidbmerge) stop("cross-merge not yet supported")
# Cast and redimension y conformably with x along join dimensions:
  idx.x = which(dimensions(x) %in% by.x)
  msk.y = dimensions(y) %in% by.y
  newds = lapply(1:length(dimensions(y)),
    function(j) {
# It's possible to get a SciDB name conflict here (issue #41).
      y_dim = dimensions(y)[j]
      y_new = make.unique_(dimensions(x),y_dim)
      if(!msk.y[j])
      {
        build_dim_schema(y,I=j,bracket=FALSE,newnames=y_new)
      } else
      {
        ind = which(by.y %in% y_dim) # by index
        build_dim_schema(x,I=idx.x[ind],newnames=y_dim,bracket=FALSE,newend=scidb_coordinate_end(y)[j])
      }
    })
  newds = newds[!unlist(lapply(newds,is.null))]
  newds = sprintf("[%s]",paste(newds,collapse=","))

# If the chunk sizes are identical, we're OK (join does not care about the
# upper array bounds). Otherwise we need redimension.
  if(isTRUE(compare_schema(x,y,ignore_attributes=TRUE,ignore_types=TRUE,ignore_nullable=TRUE,s1_dimension_index=idx.x, s2_dimension_index=which(msk.y), ignore_end=TRUE)))
  {
    castschema = sprintf("%s%s",newas,newds)
    z = cast(y, castschema)
  } else
  {
    reschema = sprintf("%s%s", newas,newds)
    castschema = sprintf("%s%s",newas,newds)
    z = redimension(cast(subarray(y,limits=reschema,between=TRUE),castschema),reschema)
  }

# Join on dimensions.
  query = sprintf("cross_join(%s as %s, %s as %s", xname, al[1], z@name, al[2])
  k = min(length(by.x),length(by.y))
  by.x = by.x[1:k]
  by.y = by.y[1:k]
  cterms = unique(paste(c(al[1],al[2]), as.vector(rbind(by.x,by.y)), sep="."))
  cterms = paste(cterms,collapse=",")
  query  = paste(query,",",cterms,")",sep="")
  ans = .scidbeval(query,eval,depend=list(x,y))
  ans
}
