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

setGeneric("c")
setGeneric('is.scidbdf', function(x) standardGeneric('is.scidbdf'))

# XXX This is a prototype. Right now the attribute list must match.
# XXX rbind on dataframes should call this
setMethod(c,signature(x="scidbdf"),
function(x,y,`eval`=FALSE)
{
  if(is.scidb(y)) y = scidb(y,`data.frame`=TRUE)
  if(!is.scidbdf(y)) y = as.scidb(y)
  if(as.numeric(scidb_coordinate_bounds(x)$length) < as.numeric(.scidb_DIM_MAX))
  {
    s = sprintf("%s%s",build_attr_schema(x),build_dim_schema(x,newend=.scidb_DIM_MAX))
    x = redimension(x,s)
  }
  i = count(x) + as.numeric(scidb_coordinate_start(x)) - as.numeric(scidb_coordinate_start(y))
  j = make.unique_(y@attributes, "j")
  fun = sprintf("%s + %.0f", dimensions(y), i)
  s = sprintf("apply(%s, %s, %s)",y@name, j, fun)
  scma = sprintf("%s%s",build_attr_schema(y), build_dim_schema(x,newname=j))
  s = sprintf("redimension(%s, %s)",s, scma)
  s = sprintf("cast(%s, %s%s)", s,build_attr_schema(y), build_dim_schema(x))
  s = sprintf("merge(%s, %s)", x@name, s)
  .scidbeval(s, `data.frame`=TRUE, gc=TRUE, `eval`=eval, depend=list(x,y))
})

# Head and tail are not very efficient, but they're nifty!  XXX
setMethod("head", signature(x="scidbdf"),
function(x, n=6L, ...)
{
  iqdf(x, n)[,-c(1,2)]
})

setMethod("tail", signature(x="scidbdf"),
function(x, n=6L, ...)
{
  ans = x[x[,1],]
  end = as.numeric(scidb_coordinate_bounds(ans)$end)
  start = max(0, end-n+1)
  ans[start:end,][]
})

setMethod("Filter",signature(f="character",x="scidbdf"),
  function(f, x)
  {
    filter_scidb(x,f)
  })

setMethod('is.scidbdf', signature(x='scidbdf'),
  function(x) return(TRUE))
setMethod('is.scidbdf', definition=function(x) return(FALSE))

setMethod('print', signature(x='scidbdf'),
  function(x) {
    show(x)
  })

setMethod("na.locf",signature(object="scidbdf"), na.locf_scidb)
setMethod("hist",signature(x="scidbdf"), hist_scidb)

setMethod('show', 'scidbdf',
  function(object) {
    v = ifelse(length(object@attributes)<2, "variable", "variables")
    l = scidb_coordinate_bounds(object)$length
    if(as.numeric(l) > 4e18) l = "*"
    cat(sprintf("SciDB 1-D array: %s obs. of %d %s.\n", l,
        length(object@attributes),v))
  })

setMethod("regrid", signature(x="scidbdf"),
  function(x, grid, expr)
  {
    if(missing(expr)) expr = paste(sprintf("max(%s)",x@attributes),collapse=",")
    query = sprintf("regrid(%s, %s, %s)",
               x@name, paste(noE(grid),collapse=","), expr)
    .scidbeval(query, eval=FALSE, gc=TRUE, depend=list(x))
  })
setMethod("xgrid", signature(x="scidbdf"),
  function(x, grid)
  {
    query = sprintf("xgrid(%s, %s)", x@name, paste(noE(grid),collapse=","))
    .scidbeval(query, eval=FALSE, gc=TRUE, depend=list(x))
  })
setMethod("unpack",signature(x="scidbdf"),unpack_scidb)

setMethod("aggregate", signature(x="scidbdf"), aggregate_scidb)
setMethod("reshape", signature(data="scidbdf"), reshape_scidb)

scidbdf_grand = function(x, op)
{
  query = sprintf("aggregate(%s, %s(%s) as %s)", x@name, op, x@attributes[1], x@attributes[1])
  iquery(query, `return`=TRUE)[,2]
}

# The following methods return data to R
setMethod("sum", signature(x="scidbdf"),
function(x)
{
  scidbdf_grand(x, "sum")
})

setMethod("median", signature(x="scidbdf"),
function(x)
{
  scidbdf_grand(x, "median")
})

setMethod("mean", signature(x="scidbdf"),
function(x)
{
  scidbdf_grand(x, "avg")
})

setMethod("min", signature(x="scidbdf"),
function(x)
{
  scidbdf_grand(x, "min")
})

setMethod("max", signature(x="scidbdf"),
function(x)
{
  scidbdf_grand(x, "max")
})

setMethod("sd", signature(x="scidbdf"),
function(x)
{
  scidbdf_grand(x, "stdev")
})

setMethod("var", signature(x="scidbdf"),
function(x)
{
  scidbdf_grand(x, "var")
})

log.scidbdf = function(x, base=exp(1))
{
  log_scidb(x,base) 
}

setMethod("sin", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "sin")
  })
setMethod("cos", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "cos")
  })
setMethod("tan", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "tan")
  })
setMethod("asin", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "asin")
  })
setMethod("acos", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "acos")
  })
setMethod("atan", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "atan")
  })
setMethod("abs", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "abs")
  })
setMethod("sqrt", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "sqrt")
  })
setMethod("exp", signature(x="scidbdf"),
  function(x)
  {
    fn_scidb(x, "exp")
  })
# Non-traditional masking binary comparison operators
setMethod("%<%",signature(x="scidbdf", y="ANY"),
  function(x,y)
  {
    .compare(x,y,"<",traditional=FALSE)
  },
  valueClass="scidbdf"
)
setMethod("%>%",signature(x="scidbdf", y="ANY"),
  function(x,y)
  {
    .compare(x,y,">",traditional=FALSE)
  },
  valueClass="scidbdf"
)
setMethod("%<=%",signature(x="scidbdf", y="ANY"),
  function(x,y)
  {
    .compare(x,y,"<=",traditional=FALSE)
  },
  valueClass="scidbdf"
)
setMethod("%>=%",signature(x="scidbdf", y="ANY"),
  function(x,y)
  {
    .compare(x,y,">=",traditional=FALSE)
  },
  valueClass="scidbdf"
)
setMethod("%==%",signature(x="scidbdf", y="ANY"),
  function(x,y)
  {
    .compare(x,y,"==",traditional=FALSE)
  },
  valueClass="scidbdf"
)
setMethod("%!=%",signature(x="scidbdf", y="ANY"),
  function(x,y)
  {
    .compare(x,y,"!=",traditional=FALSE)
  },
  valueClass="scidbdf"
)
setGeneric("head")
