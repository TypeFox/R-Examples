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

# Generic function declarations
setGeneric("%<%", def=function(x,y){NULL})
setGeneric("%>%", def=function(x,y){NULL})
setGeneric("%<=%", def=function(x,y){NULL})
setGeneric("%>=%", def=function(x,y){NULL})
setGeneric("%==%", def=function(x,y){NULL})
setGeneric("%!=%", def=function(x,y){NULL})

setOldClass("glm")
setGeneric("glm")
setMethod("glm", signature(formula="ANY", family="ANY", data="scidbdf"), glm_scidb)

setOldClass("crossprod")
setGeneric("crossprod")

setOldClass("crossprod")
setGeneric("crossprod")

setOldClass("na.locf")
setGeneric("na.locf")

setOldClass("hist")
setGeneric("hist")

setOldClass("rank")
setGeneric("rank")

setGeneric("sum")
setGeneric("mean")
setGeneric("median")
setGeneric("min")
setGeneric("max")
setGeneric("sd")
setGeneric("var")
setGeneric("diag")
setGeneric("head")
setGeneric("tail")
setGeneric("is.scidb", function(x) standardGeneric("is.scidb"))
setGeneric("print")#, function(x) standardGeneric("print"))
setGeneric("image")

setOldClass("aggregate")
setGeneric("aggregate")

setOldClass("sweep")
setGeneric("sweep")

setOldClass("apply")
setGeneric("apply")

setMethod("unpack",signature(x="scidb"),unpack_scidb)

setOldClass("reshape")
setGeneric("reshape", function(data,...) stats::reshape(data,...))

setOldClass("dist")
setGeneric("dist")

setOldClass("kmeans")
setGeneric("kmeans")

setOldClass("svd")
setGeneric("svd")

setOldClass("glm.fit")
setGeneric("glm.fit")

setOldClass("t")
setGeneric("t")

setOldClass("lag")
setGeneric("lag")
setGeneric("regrid", def=function(x,grid,expr){NULL})
setGeneric("xgrid", def=function(x,grid){NULL})
setGeneric("Filter")

# Non-traditional masking binary comparison operators
setMethod("%<%",signature(x="scidb", y="ANY"),
  function(x,y)
  {
    .compare(x,y,"<",traditional=FALSE)
  },
  valueClass="scidb"
)
setMethod("%>%",signature(x="scidb", y="ANY"),
  function(x,y)
  {
    .compare(x,y,">",traditional=FALSE)
  },
  valueClass="scidb"
)
setMethod("%<=%",signature(x="scidb", y="ANY"),
  function(x,y)
  {
    .compare(x,y,"<=",traditional=FALSE)
  },
  valueClass="scidb"
)
setMethod("%>=%",signature(x="scidb", y="ANY"),
  function(x,y)
  {
    .compare(x,y,">=",traditional=FALSE)
  },
  valueClass="scidb"
)
setMethod("%==%",signature(x="scidb", y="ANY"),
  function(x,y)
  {
    .compare(x,y,"==",traditional=FALSE)
  },
  valueClass="scidb"
)
setMethod("%!=%",signature(x="scidb", y="ANY"),
  function(x,y)
  {
    .compare(x,y,"!=",traditional=FALSE)
  },
  valueClass="scidb"
)



setMethod("%*%",signature(x="scidb", y="scidbdf"),
  function(x,y)
  {
    scidbmultiply(x,cbind(y))
  },
  valueClass="scidb"
)

setMethod("%*%",signature(x="scidb", y="scidb"),
  function(x,y)
  {
    scidbmultiply(x,y)
  },
  valueClass="scidb"
)

setMethod("%*%",signature(x="matrix", y="scidb"),
  function(x,y)
  {
    as.scidb(x,gc=TRUE) %*% y
  },
  valueClass="character"
)

setMethod("%*%",signature(x="scidb", y="matrix"),
  function(x,y)
  {
    x %*% as.scidb(y,gc=TRUE)
  },
  valueClass="character"
)

setMethod("%*%",signature(x="scidb", y="numeric"),
  function(x,y)
  {
    x %*% as.scidb(matrix(y,nrow=ncol(x)), gc=TRUE)
  },
  valueClass="character"
)

setMethod("%*%",signature(x="numeric", y="scidb"),
  function(x,y)
  {
    as.scidb(matrix(x,ncol=nrow(y)), gc=TRUE) %*% y
  },
  valueClass="character"
)

setMethod("crossprod",signature(x="scidb", y="missing"),
  function(x)
  {
    t(x) %*% x
  },
  valueClass="scidb"
)

setMethod("crossprod",signature(x="scidb", y="scidb"),
  function(x,y)
  {
    t(x) %*% y
  },
  valueClass="scidb"
)

setGeneric("tcrossprod")
setMethod("tcrossprod",signature(x="scidb", y="scidb"),
  function(x,y)
  {
    x %*% t(y)
  },
  valueClass="scidb"
)

setMethod("tcrossprod",signature(x="scidb", y="missing"),
  function(x)
  {
    x %*% t(x)
  },
  valueClass="scidb"
)


scidb_grand = function(x, op)
{
  ag = paste(sprintf("%s(%s) as %s", op, x@attributes, x@attributes),collapse=",")
  query = sprintf("aggregate(%s, %s)", x@name, ag)
  if(length(x@attributes)==1) return(iquery(query, `return`=TRUE)[,2])
  iquery(query,`return`=TRUE)
}

# The remaining functions return data to R:
setMethod("sum", signature(x="scidb"),
function(x)
{
  scidb_grand(x, "sum")
})

setMethod("mean", signature(x="scidb"),
function(x)
{
  scidb_grand(x, "avg")
})

setMethod("median", signature(x="scidb"),
function(x)
{
  scidb_grand(x, "median")
})

setMethod("min", signature(x="scidb"),
function(x)
{
  scidb_grand(x, "min")
})

setMethod("max", signature(x="scidb"),
function(x)
{
  scidb_grand(x, "max")
})

setMethod("sd", signature(x="scidb"),
function(x)
{
  scidb_grand(x, "stdev")
})

setMethod("var", signature(x="scidb"),
function(x)
{
  scidb_grand(x, "var")
})

setMethod("diag", signature(x="scidb"),
function(x)
{
  D = dim(x)
  atr = .get_attribute(x)
  dims = dimensions(x)
  bounds = scidb_coordinate_bounds(x)
  len = as.numeric(bounds$length)
  chunk = scidb_coordinate_chunksize(x)
  overlap = scidb_coordinate_overlap(x)
  xtypes = scidb_types(x)
  if(length(D)>2) stop("diag requires a matrix or vector")
# Two cases
# Case 1: Given a matrix, return its diagonal as a vector.
  if(length(D)==2)
  {
    bschema = sprintf("<%s:double>",atr)
    bschema = sprintf("%s%s",bschema,build_dim_schema(x,newstart=c(0,0)))
    mask = sprintf("build(<%s:double>[%s=0:%s,1000000,0],1)",atr,dims[1],noE(min(len)-1))
    mask = sprintf("apply(%s,%s,%s)",mask,dims[2],dims[1])
    mask = sprintf("redimension(%s,%s)",mask, bschema)
    mask = sprintf("attribute_rename(%s,%s,%s)",mask,atr,make.unique_(atr,"v"))
    GEMM.BUG = ifelse(is.logical(options("scidb.gemm_bug")[[1]]),options("scidb.gemm_bug")[[1]],FALSE)
    if(GEMM.BUG) query = sprintf("project(join(sg(subarray(%s,null,null,null,null),1,-1),%s),%s)",x@name,mask,atr)
    else query = sprintf("project(join(subarray(%s,null,null,null,null),%s),%s)",x@name,mask,atr)
    query = sprintf("unpack(%s,%s)",query, make.unique_(dims, "i"))
    query = sprintf("project(%s,%s)",query,atr)
    if(GEMM.BUG) query = sprintf("sg(subarray(%s,0,%s),1,-1)",query,noE(min(len)-1))
    else query = sprintf("subarray(%s,0,%s)",query,noE(min(len)-1))
    return(.scidbeval(query, eval=FALSE, gc=TRUE, depend=list(x)))
  }
# Case 2: Given a vector return  diagonal matrix.
  dim2 = make.unique_(dims, "j")
  query = sprintf("build(<%s:%s>[%s=0:%s,%s,%s],nan)",
           make.unique_(x@attributes,"v"), xtypes,
           dims[1], noE(len[1] - 1), chunk[1], overlap[1])
  query = sprintf("apply(%s,%s,%s)",query,dim2,dims[1])
  query = sprintf("redimension(%s,<%s:%s>[%s=0:%s,%s,%s,%s=0:%s,%s,%s])",
           query, make.unique_(x@attributes,"v"), scidb_types(x),
           dims[1], noE(len[1] - 1) , chunk[1], overlap[1],
           dim2, noE(len[1] - 1), chunk[1], overlap[1])
  query = sprintf("cross_join(%s as __X,%s as __Y,__X.%s,__Y.%s)",query,x@name,dims[1],dims[1])
  query = sprintf("project(%s,%s)",query,atr)
  ans = .scidbeval(query, eval=FALSE, gc=TRUE, depend=list(x))
  attr(ans, "sparse") = TRUE
  return(ans)
})

setMethod("head", signature(x="scidb"),
function(x, n=6L, ...)
{
  class(x)="scidbdf"  # to supress warning
  iqdf(x,n)
})

setMethod("tail", signature(x="scidb"),
function(x, n=6L, ...)
{
  bounds = scidb_coordinate_bounds(x)
  xstart = as.numeric(bounds$start)
  xlen   = as.numeric(bounds$length)
  p = xstart + xlen - 1
  m = xstart + xlen - n
  m = unlist(lapply(1:length(m),function(j) max(m[j],xstart[j])))
  limits = lapply(1:length(m), function(j) seq(m[j],p[j]))
  do.call(dimfilter,args=list(x=x,i=limits,eval=FALSE,drop=TRUE))[]
})

setMethod("is.scidb", signature(x="ANY"),
  function(x) 
  {
    if(inherits(x, "scidb")) return(TRUE)
    FALSE
  }
)

setMethod("print", signature(x="scidb"),
  function(x) {
    show(x)
  })

setMethod("show", "scidb",
  function(object) {
    if(is.null(dim(object)) || length(dim(object))==1)
      cat("Reference to a SciDB vector of length",
           scidb_coordinate_bounds(object)$length,"\n")
    else
      cat("A reference to a ",
          paste(dim(object),collapse="x"),
          "SciDB array\n")
  })

setMethod("image", signature(x="scidb"),
function(x, grid=c(500,500), op=sprintf("sum(%s)", .get_attribute(x)), na=0, ...)
{
  if(length(dim(x))!=2) stop("Sorry, array must be two-dimensional")
  if(length(grid)!=2) stop("The grid parameter must contain two values")
  if(any(is.infinite(dim(x)))) x = bound(x)
  el = list(...)
  if("plot" %in% names(el))
  {
    plot = as.logical(el$plot)
  }
  else plot=TRUE
  blocks = as.numeric(scidb_coordinate_bounds(x)$length)
  blocks = blocks/grid
  if(any(blocks<1)) blocks[which(blocks<1)] = 1
  query = sprintf("regrid(project(%s,%s),%.0f,%.0f,%s)",x@name,.get_attribute(x),blocks[1],blocks[2],op)
  A = iquery(query,return=TRUE,n=Inf)
  A[is.na(A[,3]),3] = na
  m = max(A[,1]) + 1
  n = max(A[,2]) + 1
  B = matrix(0,m,n)
  B[A[,1] + A[,2]*m + 1] = A[,3]
  if(!plot) return (B)
  xlbl=(1:ncol(B))*blocks[2]
  xat=seq(from=0,to=1,length.out=ncol(B))
  ylbl=(nrow(B):1)*blocks[1]
  yat=seq(from=0,to=1,length.out=nrow(B))
  image(B,xaxt='n',yaxt='n',...)
  axis(side=1,at=xat,labels=xlbl)
  axis(side=2,at=yat,labels=ylbl)
  B
})

setMethod("aggregate", signature(x="scidb"), aggregate_scidb)
setMethod("sweep", signature(x="scidb"), sweep_scidb)
setMethod("apply", signature(X = "scidb"), apply_scidb)
setMethod("unpack",signature(x="scidb"),unpack_scidb)
setMethod("reshape", signature(data="scidb"), reshape_scidb)
setMethod("dist", signature(x="scidb"), dist_scidb)
setMethod("kmeans", signature(x="scidb"), kmeans_scidb)
setMethod("svd", signature(x="scidb"), svd_scidb)
setMethod("glm.fit", signature(x="scidb",y="ANY",weights="MNSN"), glm.fit_scidb)
setMethod("na.locf",signature(object="scidb"), na.locf_scidb)
setMethod("hist",signature(x="scidb"), hist_scidb)
setMethod("rank",signature(x="scidb"), rank_scidb)
setGeneric("order")
setMethod("order",signature("ANY"), order_scidb)


# Transpose a matrix or vector
setMethod("t", signature(x="scidb"), 
  function(x)
  {
    query = sprintf("transpose(%s)",x@name)
    .scidbeval(query, eval=FALSE, gc=TRUE, depend=list(x))
  }
)

# Lead or lag a time series
setMethod("lag",signature(x="scidb"),
  function(x,k=1,dim=1,eval=FALSE)
  {
    names = dimensions(x)
    start = scidb_coordinate_start(x)
    xstart = start
    n = make.unique_(c(scidb_attributes(x), names), "n")
    expr = sprintf("%s - %s", names[dim], k)
    y = bind(x,n,expr)
    start[dim] = noE(as.numeric(start[dim]) - k)
    names[dim] = n
    schema = sprintf("%s%s", build_attr_schema(x),
      build_dim_schema(x,newstart=start,newnames=names))
    y = redimension(y,schema)
    cschema = sprintf("%s%s",build_attr_schema(x),
                build_dim_schema(y,newnames=dimensions(x)))
    y = cast(y,cschema)
    query = sprintf("between(%s,%s)",y@name,between_coordinate_bounds(x))
    query = sprintf("redimension(%s,%s%s)",query, build_attr_schema(x),build_dim_schema(x))
    .scidbeval(query,eval=FALSE,gc=TRUE,depend=list(x,y))
  })

# SciDB's regrid and xgrid operators (simple wrappers)
setMethod("regrid", signature(x="scidb"),
  function(x, grid, expr)
  {
    if(missing(expr)) expr = paste(sprintf("avg(%s)",x@attributes),collapse=",")
    query = sprintf("regrid(%s, %s, %s)",
               x@name, paste(noE(grid),collapse=","), expr)
    .scidbeval(query, eval=FALSE, gc=TRUE, depend=list(x))
  })
setMethod("xgrid", signature(x="scidb"),
  function(x, grid)
  {
    query = sprintf("xgrid(%s, %s)", x@name, paste(noE(grid),collapse=","))
    .scidbeval(query, eval=FALSE, gc=TRUE, depend=list(x))
  })

setMethod("Filter",signature(f="character",x="scidb"),
  function(f, x)
  {
    filter_scidb(x,f)
  })

setMethod("sin",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "sin")
  })
setMethod("cos",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "cos")
  })
setMethod("tan",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "tan")
  })
setMethod("asin",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "asin")
  })
setMethod("acos",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "acos")
  })
setMethod("atan",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "atan")
  })
setMethod("abs",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "abs")
  })
setMethod("sqrt",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "sqrt")
  })
setMethod("exp",signature(x="scidb"),
  function(x)
  {
    fn_scidb(x, "exp")
  })
