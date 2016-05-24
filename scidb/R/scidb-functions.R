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

# Various functions that support S3 methods for the SciDB class

# cheeky solve
solve.scidb = function(a, b, ...)
{
  A = scidbeval(crossprod(a))
  s = svd(A)
  x = scidbeval(t(t(b) %*% a))
  x = t(s$u) %*% x
  x = s$d^(-1) * x
  x = s$v %*% x
  scidbeval(x)
}

cbind.scidb = function(x)
{
  if(length(dim(x))!=1) return(x)
  newdim=make.unique_(c(dimensions(x),x@attributes), "j")
  nd = sprintf("%s[%s,%s=0:0,1,0]",build_attr_schema(x) , build_dim_schema(x,bracket=FALSE),newdim)
  redimension(x, nd)
}

log.scidb = function(x, base=exp(1))
{
  log_scidb(x,base) 
}

colnames.scidb = function(x)
{
  if(is.null(x@gc$dimnames)) return(NULL)
  if(length(x@gc$dimnames)<2) return(NULL)
  x@gc$dimnames[[2]]
}

`colnames<-.scidb` = function(x, value)
{
  x
}

rownames.scidb = function(x)
{
  if(is.null(x@gc$dimnames)) return(NULL)
  x@gc$dimnames[[1]]
}

`rownames<-.scidb` = function(x, value)
{
  x
}

names.scidb = function(x)
{
  if(is.null(dim(x))) rownames(x)
  colnames(x)
}

`names<-.scidb` = function(x, value)
{
  colnames(x) <- value
}

dimnames.scidb = function(x)
{
  x@gc$dimnames
}

# SciDB labeled coordinates assignment
`dimnames<-.scidb` = function(x, value)
{
  if(is.null(value)) 
  {
    x@gc$dimnames = value
    return(x)
  }
  if(!is.list(value))
    stop("dimnames requires a list")
  if(length(value)!=length(dim(x)))
    stop("incorrect number of dimensions specified")
# Make sure that specified indices are all SciDB arrays.
# If not, try to make them so.
  value = lapply(1:length(value), function(j)
    {
      v = value[[j]]
      if(is.null(v)) return(v)
      if(is.scidbdf(v)) class(v) = "scidb"
      if(is.scidb(v))
      {
        if(length(dim(v))>1) stop("Dimension label arrays must be one dimensional")
        if(length(v@attributes)>1) v = project(v,1)
# Make sure that the label array has '*' upper dimension bound
        check = scidb_coordinate_start(x)[j] == scidb_coordinate_start(v)[1] &&
                scidb_coordinate_chunksize(x)[j] == scidb_coordinate_chunksize(v)[1] &&
                (is.na(dim(v)) || (dim(v) == as.numeric(.scidb_DIM_MAX)))
        if(!check)
        {
# uh oh
            v = redimension(v, sprintf("%s%s", build_attr_schema(v), build_dim_schema(v, newend="*")))
            v = reshape(v, sprintf("%s%s", build_attr_schema(v), build_dim_schema(v, newstart=scidb_coordinate_start(x)[j])))
        }
# join out holes. don't use dimnames! this is crummy.
        xj = apply(x,j,count)
        v = project(merge(v, xj, by.x=dimensions(v), by.y=dimensions(xj)),1)
        return(v);
      }
      scidbeval(unbound(as.scidb(v,
                start=scidb_coordinate_start(x)[j],
                chunkSize=scidb_coordinate_chunksize(x)[j])))
    })

  check = unlist(lapply(1:length(value), function(j)
    {
      is.null(value[[j]]) || (scidb_coordinate_start(value[[j]])[1] == scidb_coordinate_start(x)[j])
    }))
  if(!all(check))
    stop("Label array starting indices don't match data array--please check")
  check = unlist(lapply(1:length(value), function(j)
    {
      is.null(value[[j]]) || (nrow(value[[j]]) == dim(x)[j])
    }))
  check[is.na(check)] = FALSE
  x@gc$dimnames = value
  x@gc$depend =c (x@gc$depend, unlist(value))
  x
}

summary.scidb = function(x)
{
  warning("Not available.")
  invisible()
}

`[<-.scidb` = function(x, ..., value)
{
  m = match.call()
  a = scidb_attributes(x)[1]
  if(!is.null(m$attr)) a = m$attr
  ai = which(scidb_attributes(x) %in% a)
  i = list(...)
  if(!all(unlist(lapply(i,checkseq)))) stop("Assignment is limited to contiguous blocks for now.")
  d = sapply(i, length) # dimensions
  s = sapply(i, function(u) u[1]) # origin
# Note, ordered by rows thanks to aperm
  if(!is.array(value)) value = array(value)
  v = as.scidb(as.vector(aperm(value)), nullable=scidb_nullable(x)[ai], attr=a, reshape=FALSE)
  reschema = sprintf("%s%s", build_attr_schema(v),
              build_dim_schema(x, newstart=s, newlen=d))
  merge(x, redimension(reshape(v, schema=reschema), schema=schema(x)), merge=TRUE)
}

# Array subsetting wrapper.
# x: A Scidb array object
# ...: list of dimensions
# 
# Returns a materialized R array if length(list(...))==0.
# Or, a scidb array promise.
`[.scidb` = function(x, ...)
{
  M = match.call()
  if(is.null(M$drop)) drop=TRUE  # if TRUE follow R drop convention
  else drop=M$drop
  if(is.null(M$eval)) eval=FALSE  # if TRUE eval (not really needed anymore)
  else eval=M$eval
  if(is.null(M$redim)) redim=TRUE  # if FALSE don't reset coordinates
  else redim=M$redim
# inverse of redim (same functionality) for user convenience
  if(!is.null(M$between)) redim=!M$between
  if(!is.null(M$redim) && !is.null(M$between))
  {
    warning("The `redim` and `between` options are antonyms but they provide the same control over the query. When both are specified the `redim` option is ignored and subsetting only pays attention to the between option.")
  }
  M = M[3:length(M)]
  if(!is.null(names(M))) M = M[!(names(M) %in% c("drop","eval","redim","between"))]
# i shall contain a list of requested index values
  E = parent.frame()
  i = lapply(1:length(M), function(j) tryCatch(eval(M[j][[1]],E),error=function(e)c()))
# User wants this materialized to R...
  if(all(sapply(i,is.null)))
    return(materialize(x,drop=drop))
# Not materializing, return a SciDB array
  if(length(i)!=length(dim(x))) stop("Dimension mismatch")
  dimfilter(x,i,eval,drop=drop,redim)
}

`dim<-.scidb` = function(x, value)
{
  reshape(x, shape=value)
}

str.scidb = function(object, ...)
{
  .scidbstr(object)
}

print.scidb = function(x, ...)
{
  cat("\nSciDB array ",x@name)
  cat("\n")
  print(head(x))
  if(is.null(dim(x))) j = as.numeric(scidb_coordinate_bounds(x)$length) - 6
  else j = dim(x)[1]-6
  if(j>2) cat("and ",j,"more rows not displayed...\n")
  if(length(dim(x))>0) {
    k = dim(x)[2] - 6
    if(k>2) cat("and ",k,"more columns not displayed...\n")
  }
}

ncol.scidb = function(x) dim(x)[2]
nrow.scidb = function(x) dim(x)[1]
dim.scidb = function(x) as.numeric(scidb_coordinate_bounds(x)$length)
length.scidb = function(x) prod(as.numeric(scidb_coordinate_bounds(x)$length))

# Vector, Matrix, matrix, or data.frame only.
# XXX Future: Add n-d array support here (TODO)
as.scidb = function(X,
                    name=tmpnam(),
                    chunksize,
                    overlap,
                    start,
                    chunkSize = chunksize,
                    gc=TRUE, ...)
{
  if(inherits(X,"data.frame"))
  {
    if(missing(chunksize))
      return(df2scidb(X,name=name,gc=gc,start=start,...))
    else
      return(df2scidb(X,name=name,chunkSize=as.numeric(chunksize[[1]]),gc=gc,start=start,...))
  }
  force_type = NULL
  if(is.factor(X))
  {
    if("scidb_factor" %in% class(X))
    {
      levels(X) = levels_scidb(X)
      X = as.integer(as.vector(X))  # XXX !!! XXX THIS LIMITS INDEX VALUES TO 31 bits!!! !!!
      force_type = "int64"
    }
    else
      X = as.vector(X)
  }
  args = list(...)
  nullable = TRUE
  if(!is.null(args$nullable)) nullable = as.logical(args$nullable) # control nullability
  attr_name = "val"
  if(!is.null(args$attr)) attr_name = as.character(args$attr)      # attribute name
  do_reshape = TRUE
  if(!is.null(args$reshape)) do_reshape = as.logical(args$reshape) # control reshape
  if(missing(chunksize))
  {
# Note nrow, ncol might be NULL here if X is not a matrix. That's OK, we'll
# deal with that case later.
    chunkSize=c(min(1000L,nrow(X)),min(1000L,ncol(X)))
  }
  chunkSize = as.numeric(chunkSize)
  if(length(chunkSize)==1) chunkSize = c(chunkSize, chunkSize)
  if(!missing(overlap)) warning("Sorry, overlap is not yet supported by the as.scidb function. Consider using the reparition function for now.")
  overlap = c(0,0)
  if(missing(start)) start=c(0,0)
  start     = as.numeric(start)
  if(length(start)==1) start=c(start,start)
  if(inherits(X,"dgCMatrix"))
  {
# Sparse matrix case
    return(.Matrix2scidb(X,name=name,rowChunkSize=chunkSize[[1]],colChunkSize=chunkSize[[2]],start=start,gc=gc,...))
  }
  D = dim(X)
  start = as.integer(start)
  overlap = as.integer(overlap)
  type = .scidbtypes[[typeof(X)]]
  if(is.null(type)) {
    stop(paste("Unupported data type. The package presently supports: ",
       paste(.scidbtypes,collapse=" "),".",sep=""))
   }
  if(is.null(D)) {
# X is a vector
    if(!is.vector(X)) stop ("X must be a matrix or a vector")
    chunkSize = min(chunkSize[[1]],length(X))
    X = as.matrix(X)
    schema = sprintf(
      "< %s : %s null>  [i=%.0f:%.0f,%.0f,%.0f]", attr_name, type, start[[1]],
      nrow(X)-1+start[[1]], min(nrow(X),chunkSize), overlap[[1]])
    load_schema = schema
    if(!is.null(force_type))
    {
      schema = sprintf(
        "< %s : %s null>  [i=%.0f:%.0f,%.0f,%.0f]", attr_name, force_type, start[[1]],
        nrow(X)-1+start[[1]], min(nrow(X),chunkSize), overlap[[1]])
    }
  } else {
# X is a matrix
    schema = sprintf(
      "< %s : %s  null>  [i=%.0f:%.0f,%.0f,%.0f, j=%.0f:%.0f,%.0f,%.0f]", attr_name, type, start[[1]],
      nrow(X)-1+start[[1]], chunkSize[[1]], overlap[[1]], start[[2]], ncol(X)-1+start[[2]],
      chunkSize[[2]], overlap[[2]])
    load_schema = sprintf("<%s:%s null>[__row=1:%.0f,1000000,0]",attr_name, type,  length(X))
  }
  if(!is.matrix(X)) stop ("X must be a matrix or a vector")

  DEBUG = FALSE
  if(!is.null(options("scidb.debug")[[1]]) && TRUE==options("scidb.debug")[[1]]) DEBUG=TRUE
  td1 = proc.time()
# Obtain a session from shim for the upload process
  session = getSession()
  on.exit( GET("/release_session",list(id=session)) ,add=TRUE)

# Upload the data
# XXX I couldn't get RCurl to work using the fileUpload(contents=x), with 'x'
# a raw vector. But we need RCurl to support SSL. As a work-around, we save
# the object. This extra local copy sucks and must be improved !!! XXX
# XXX The bug is in RCurl's curl.c addFormElement function.
# XXX
  fn = tempfile()
  bytes = writeBin(.Call("scidb_raw",as.vector(t(X)),PACKAGE="scidb"),con=fn)
  url = URI("upload_file",list(id=session))
# RCurl madness! postForm with digest auth errors out, even when it
# succeeds! XXX
  hdr = digest_auth("POST",url)
  ans = tryCatch(postForm(uri = url, uploadedfile = fileUpload(filename=fn), .opts = c(curlopts(),curlOptions(httpheader=hdr,'ssl.verifypeer'=0,'ssl.verifyhost'=as.integer(options("scidb.verifyhost"))))), error=invisible)
  unlink(fn)
  ans = ans[[1]]
  ans = gsub("\r","",ans)
  ans = gsub("\n","",ans)
  if(DEBUG)
  {
    cat("Data upload time",(proc.time()-td1)[3],"\n")
  }

# Load query
  if(do_reshape)
  {
    if(!is.null(force_type))
    {
      query = sprintf("store( attribute_rename(project(apply(reshape(input(%s,'%s', 0, '(%s null)'),%s), ___, %s(%s)), ___), ___, %s), %s)",load_schema,ans,type,schema, force_type, attr_name, attr_name, name)
    }
    else {
      query = sprintf("store(reshape(input(%s,'%s', 0, '(%s null)'),%s),%s)",load_schema,ans,type,schema,name)
    }
  }
  else
  {
    if(!is.null(force_type))
    {
      query = sprintf("store(attribute_rename(project(appy(input(%s,'%s', 0, '(%s null)'),___, %s(%s)), ___), ___, %s), %s)",load_schema,ans,type,name,force_type,attr_name,attr_name)
    } else
    {
      query = sprintf("store(input(%s,'%s', 0, '(%s null)'),%s)",load_schema,ans,type,name)
    }
  }
  iquery(query)
  ans = scidb(name,gc=gc)
  if(!nullable) ans = replaceNA(ans)
  ans
}
