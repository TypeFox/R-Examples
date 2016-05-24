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
Ops.scidb = function(e1,e2) {
  switch(.Generic,
    '^' = .binop(e1,e2,"^"),
    '+' = .binop(e1,e2,"+"),
    '-' = .binop(e1,e2,"-"),
    '*' = .binop(e1,e2,"*"),
    '/' = .binop(e1,e2,"/"),
    '<' = .compare(e1,e2,"<"),
    '<=' =.compare(e1,e2,"<="),
    '>' = .compare(e1,e2,">"),
    '>=' = .compare(e1,e2,">="),
    '==' = .compare(e1,e2,"="),
    '!=' = .compare(e1,e2,"<>"),
    default = stop("Unsupported binary operation.")
  )
}

# e1 and e2 must each already be SciDB arrays.
scidbmultiply = function(e1,e2)
{
# As of SciDB version 13.12, SciDB exhibits nasty bugs when gemm is nested
# within other SciDB operators, in particular subarray. We use sg to avoid
# this problem.
  GEMM.BUG = ifelse(is.logical(options("scidb.gemm_bug")[[1]]),options("scidb.gemm_bug")[[1]],FALSE)
  `eval` = FALSE
# Check for availability of spgemm
  SPGEMM = length(grep("spgemm",.scidbenv$ops[,2]))>0
  a1 = .get_attribute(e1)
  a2 = .get_attribute(e2)
  e1.sparse = is.sparse(e1)
  e2.sparse = is.sparse(e2)
  SPARSE = e1.sparse || e2.sparse

# Up to at least SciDB 13.12, gemm does not accept nullable attributes.
  if(any(scidb_nullable(e1)))
  {
    warnonce("missing")
    e1 = replaceNA(e1)
  }
  if(any(scidb_nullable(e2)))
  {
    warnonce("missing")
    e2 = replaceNA(e2)
  }

# Promote vectors to row- or column-vectors as required.
  if(length(dim(e1))<2)
  {
    e2len = as.numeric(scidb_coordinate_bounds(e2)$length)
    e1chunk = scidb_coordinate_chunksize(e1)
    e2chunk = scidb_coordinate_chunksize(e2)
    L = dim(e1)
    M = L/e2len[1]
    if(M != floor(M)) stop("Non-conformable dimensions")

    as = build_attr_schema(e1)
    osc = sprintf("%s[i=0:%s,%s,0,j=0:%s,%s,0]",
            as, noE(M-1), noE(e1chunk[1]),
              noE(e2len[1]-1), noE(e2chunk[1]))
    e1 = reshape(e1, schema=osc)
  }
  if(length(dim(e2)) < 2)
  {
    e1len = as.numeric(scidb_coordinate_bounds(e1)$length)
    e1chunk = scidb_coordinate_chunksize(e1)
    e2chunk = scidb_coordinate_chunksize(e2)
    L = dim(e2)
    N = L/e1len[2]
    if(N != floor(N)) stop("Non-conformable dimensions")

    as = build_attr_schema(e2)
    osc = sprintf("%s[i=0:%s,%s,0,j=0:%s,%s,0]",
            as, noE(e1len[2]-1), noE(e1chunk[2]),
                     noE(N-1), noE(e2chunk[1]))
    e2 = reshape(e2, schema=osc)
  }

# We use subarray to handle starting index mismatches (subarray
# returns an array with dimension indices starting at zero).
  l1 = length(dim(e1))
  lb = paste(rep("null",l1),collapse=",")
  ub = paste(rep("null",l1),collapse=",")
  if(GEMM.BUG) op1 = sprintf("sg(subarray(%s,%s,%s),1,-1)",e1@name,lb,ub)
  else op1 = sprintf("subarray(%s,%s,%s)",e1@name,lb,ub)
  l2 = length(dim(e2))
  lb = paste(rep("null",l2),collapse=",")
  ub = paste(rep("null",l2),collapse=",")
  if(GEMM.BUG) op2 = sprintf("sg(subarray(%s,%s,%s),1,-1)",e2@name,lb,ub)
  else op2 = sprintf("subarray(%s,%s,%s)",e2@name,lb,ub)

  if(!SPARSE)
  {
# Adjust the arrays to conform to GEMM requirements
    dnames = make.names_(c(dimensions(e1)[[1]],dimensions(e2)[[2]]))
    CHUNK_SIZE = options("scidb.gemm_chunk_size")[[1]]
    op1 = sprintf("repart(%s,<%s:%s>[%s=0:%s,%s,0,%s=0:%s,%s,0])",op1,a1,scidb_types(e1)[1],
            dimensions(e1)[1],noE(as.numeric(scidb_coordinate_bounds(e1)$length[1])-1),noE(CHUNK_SIZE),
            dimensions(e1)[2],noE(as.numeric(scidb_coordinate_bounds(e1)$length[2])-1),noE(CHUNK_SIZE))
    op2 = sprintf("repart(%s,<%s:%s>[%s=0:%s,%s,0,%s=0:%s,%s,0])",op2,a2,scidb_types(e2)[1],
            dimensions(e2)[1],noE(as.numeric(scidb_coordinate_bounds(e2)$length[1])-1),noE(CHUNK_SIZE),
            dimensions(e2)[2],noE(as.numeric(scidb_coordinate_bounds(e2)$length[2])-1),noE(CHUNK_SIZE))
    osc = sprintf("<%s:%s>[%s=0:%s,%s,0,%s=0:%s,%s,0]",a1,scidb_types(e1)[1],
              dnames[[1]],noE(as.numeric(scidb_coordinate_bounds(e1)$length[1])-1),noE(CHUNK_SIZE),
              dnames[[2]],noE(as.numeric(scidb_coordinate_bounds(e2)$length[2])-1),noE(CHUNK_SIZE))
    op3 = sprintf("build(%s,0)",osc)
  } else
  {
# Adjust array partitions as required by spgemm
    op2 = sprintf("repart(%s, <%s:%s>[%s=0:%s,%s,0,%s=0:%s,%s,0])",
            op2, a2, scidb_types(e2)[1],
            dimensions(e2)[1],noE(as.numeric(scidb_coordinate_bounds(e2)$length[1])-1),
                              noE(scidb_coordinate_chunksize(e1)[2]),
            dimensions(e2)[2], noE(as.numeric(scidb_coordinate_bounds(e2)$length[2])-1),
                              noE(scidb_coordinate_chunksize(e2)[2]))
  }

# Decide which multiplication algorithm to use
  if(SPARSE && !SPGEMM)
  {
    stop("Sparse matrix multiplication not supported")
  }
  else if (SPARSE && SPGEMM)
  {
    query = sprintf("spgemm(%s, %s)", op1, op2)
  }
  else
  {
    query = sprintf("gemm(%s, %s, %s)",op1,op2,op3)
    if(GEMM.BUG) query = sprintf("sg(gemm(%s, %s, %s),1,-1)",op1,op2,op3)
  }
  ans = .scidbeval(query,gc=TRUE,eval=eval,depend=list(e1,e2))
# Some SciDB array operators produce invalid-named output. Here is a fix:
  if(length(unique(dimensions(ans))) != length(dimensions(ans)))
  {
    new_dimension_names = make.unique_(c(dimensions(ans),ans@attributes),dimensions(ans))
    ans = dimension_rename(ans, new=new_dimension_names)
  }
  ans
}

# Element-wise binary operations
# These are terribly inefficient for various reasons.
.binop = function(e1,e2,op)
{
  e1a = "scalar"
  e2a = "scalar"
  dnames = c()
  depend = c()
# Handle unary minus
  if(missing(e2) && op=="-")
  {
    e2 = -1
    op = "*"
  }
# Check for non-scalar, non-scidb object arguments and convert to scidb
  if(!inherits(e1,"scidb") && length(e1)>1)
  {
    e1 = as.scidb(e1)
  }
  if(!inherits(e2,"scidb") && length(e2)>1)
  {
    e2 = as.scidb(e2)
  }
  if(inherits(e1,"scidb"))
  {
    e1 = make_nullable(e1)
    e1a = .get_attribute(e1)
    depend = c(depend, e1)
    dnames = c(dnames, dimensions(e1))
  }
  if(inherits(e2,"scidb"))
  {
    e2 = make_nullable(e2)
    e2a = .get_attribute(e2)
    depend = c(depend, e2)
    dnames = c(dnames, dimensions(e2))
  }
# Determines if we need to fill in sparse entries or not:
  fill = op %in% c("+","-","/")
  fill_div = op %in% "/"
  fill = fill && (is.sparse(e1,count=FALSE) || is.sparse(e2,count=FALSE))

# v holds new attribute name
  v = make.unique_(c(e1a,e2a,dnames),"v")
# Handle special scalar multiplication cases, including sparse cases:
  if(is.null(dim(e1)))
  {
    e2 = project(e2,e2a)
    if(op=="^")
    {  
      op = sprintf("pow(%s)",paste(c(sprintf("%.15f",e1), e2a),collapse=","))
      if(e1<0) e2 = merge(e2,build(0,e2),merge=TRUE)
    }
    else
    {
      op = paste(c(sprintf("%.15f",e1),e2a),collapse=op)
    }
    if(fill) return(project(bind(merge(e2,build(0,e2),merge=TRUE),v,op),v, eval=TRUE))
    return(project(bind(e2,v,op),v))
  }
  if(is.null(dim(e2)))
  {
    e1 = project(e1,e1a)
    if(op=="^")
    {
      op = sprintf("pow(%s)",paste(c(e1a,sprintf("%.15f",e2)),collapse=","))
      if(e2<0) e1 = merge(e1,build(0,e1),merge=TRUE)
    }
    else
      op = paste(c(e1a, sprintf("%.15f",e2)),collapse=op)
    if(fill) return(project(bind(merge(e1,build(0,e1),merge=TRUE),v,op),v, eval=TRUE))
    else return(project(bind(e1,v,op),v, eval=TRUE))
  }

# OK, we've got two scidb arrays, op them.
  l1 = length(dim(e1))
  l2 = length(dim(e2))

# Re-write op in a form compatible with SciDB apply
  rewrite_op = function(A, op)
  {
    if(op=="^")
    {
      op=sprintf("pow(%s)",paste(A@attributes,collapse=","))
    } else
    {
       op = paste(A@attributes, collapse=op)
    }
  }
# We *can't* avoid fully dense fills with "/"
  if(fill_div)
  {
    e1 = merge(e1,build(0,e1),merge=TRUE)
    e2 = merge(e2,build(0,e2),merge=TRUE)
  }
# Check special case row/column vector and orinary vector dim mismatch
  if(l1>l2 && dim(e1)[1] == dim(e2)[1] && dim(e1)[2] == 1)
  {
    e1 = reshape(e1,e2)
    l1 = length(dim(e1))
  }
  if(l2>l1 && dim(e2)[1] == dim(e1)[1] && dim(e2)[2] == 1)
  {
    e2 = reshape(e2,e1)
    l2 = length(dim(e2))
  }
# Handle conformable-dimension case
  if(l1 == l2)
  {
# Handle mis-matched coordinates :/ arrggh
    if(!compare_schema(e1,e2, ignore_dimnames=TRUE,
                              ignore_attributes=TRUE,
                              ignore_nullable=TRUE,
                              ignore_types=TRUE))
    {
      schema = sprintf("%s%s",build_attr_schema(e2),build_dim_schema(e1))
      e2 = reshape(e2,schema)
    }
# Note that we use outer join here with the special fillin option.
    M = merge(e1,e2,by.x=dimensions(e1),by.y=dimensions(e2),all=fill,fillin=0)
    v = make.unique_(c(M@attributes),"v")
    op = rewrite_op(M, op)
    return(project(bind(M,v,op), v)) # see note at end
  }

# Left now with very special vector-recycling cases.
  if(l1==1)
  {
# Handle special case similar to, but a bit different than vector recylcing.
# This case requires a dimensional match along the 1st dimensions, and it's
# useful for matrix row scaling. This is limited to vectors that match the
# number of rows of the array.
# 
# The default case:
# Handle special case similar to, but a bit different than vector recylcing.
# This case requires a dimensional match along the 2nd dimension, and it's
# useful for matrix column scaling. This is not R standard but very useful.
    x  = e1
    e1 = e2
    e2 = x
    newschema = build_dim_schema(e1,I=1,newnames=dimensions(e2)[1])
    newschema = sprintf("%s%s",build_attr_schema(e2),newschema)
    e2 = reshape(e2, newschema)
    along = dimensions(e1)[1]
  } else
  {
    newschema = build_dim_schema(e1,I=2,newnames=dimensions(e2)[1])
    newschema = sprintf("%s%s",build_attr_schema(e2),newschema)
    e2 = reshape(e2, newschema)
    along = dimensions(e1)[2]
  }
  if(fill) e1 = merge(e1, build(0,e1), merge=TRUE)
  M = merge(e1,e2,by.x=along,by.y=dimensions(e2))
  v = make.unique_(c(M@attributes),"v")
  op = rewrite_op(M,op)
# XXX We eval here because this is so inefficient...
  project(bind(M, v, op), v, eval=TRUE)
}

# Very basic comparisons. See also filter.
# e1: A scidb array
# e2: A scalar or a scidb array. If a scidb array, the return
# .joincompare(e1,e2,op) (q.v.)
# op: A comparison infix operator character
#
# Return a scidb object
# Can throw a query error.
.compare = function(e1,e2,op,traditional)
{
  if(missing(traditional)) traditional=TRUE
  if(!(inherits(e1,"scidb") || inherits(e1,"scidbdf"))) stop("Sorry, not yet implemented.")
  if(inherits(e2,"scidb")) return(.joincompare(e1,e2,op))
  op = gsub("==","=",op,perl=TRUE)
  op = gsub("!=","<>",op,perl=TRUE)
# Automatically quote characters
  if(is.character(e2)) e2 = sprintf("'%s'",e2)
  q1 = paste(paste(e1@attributes,op,e2),collapse=" and ")
# Traditional R comparisons return an array of the same shape with a true/false
# value.
  if(traditional)
  {
    newattr = make.unique_(e1@attributes, "condition")
    query = sprintf("project(apply(%s, %s, %s), %s)", e1@name, newattr, q1, newattr)
    return(.scidbeval(query, eval=FALSE, gc=TRUE, depend=list(e1)))
  }
# Alternate comparisons return a sparse mask
  query = sprintf("filter(%s, %s)", e1@name, q1)
  .scidbeval(query, eval=FALSE, gc=TRUE, depend=list(e1))
}

.joincompare = function(e1,e2,op)
{
  stop("Yikes! Not implemented yet...")
}

tsvd = function(x,nu,tol=0.1,maxit=20,tx,v,pca=FALSE)
{
  m = ceiling(1e6/nrow(x))
  n = ceiling(1e6/ncol(x))
  if(missing(v))
  {
    v = build(1,ncol(x),type="double",chunksize=ncol(x))@name
  } else v=v@name
  if(pca)
  {
    w = build(1,nrow(x),type="double",chunksize=m)@name
    v = paste(c(v,sprintf("%s, substitute(project(apply(aggregate(%s,sum(%s) as colsum,count(%s) as colcount,%s),colmean,colsum/colcount),colmean),build(<v:double>[i=0:0,1,0],0))",w,x@name, x@attributes[1], x@attributes[1], dimensions(x)[2])),collapse=",")
  }
  if(!missing(tx))
  {
    query  = sprintf("tsvd(%s, %s, %.0f, %f, %.0f, %s)", x@name, tx@name, nu,tol,maxit, v)
  } else
  {
    schema = sprintf("[%s=0:%s,%s,0,%s=0:%s,%s,0]",
                       dimensions(x)[1], noE(nrow(x)-1), noE(m),
                       dimensions(x)[2], noE(ncol(x)-1), noE(ncol(x)))
    tschema = sprintf("[%s=0:%s,%s,0,%s=0:%s,%s,0]",
                       dimensions(x)[2], noE(ncol(x)-1), noE(n),
                       dimensions(x)[1], noE(nrow(x)-1), noE(nrow(x)))
    schema = sprintf("%s%s",build_attr_schema(x), schema)
    tschema = sprintf("%s%s",build_attr_schema(x), tschema)
    query  = sprintf("tsvd(redimension(unpack(%s,row),%s), redimension(unpack(transpose(%s),row),%s), %.0f, %f, %.0f, %s)", x@name, schema, x@name, tschema, nu,tol,maxit, v)
  }
  narray = .scidbeval(query, eval=TRUE, gc=TRUE)
  ans = list(u=slice(narray, "matrix", 0,eval=FALSE)[,between(0,nu-1)],
             d=slice(narray, "matrix", 1,eval=FALSE)[between(0,nu-1),between(0,nu-1)],
             v=slice(narray, "matrix", 2,eval=FALSE)[between(0,nu-1),],
             narray=narray)
  attr(ans$u,"sparse") = TRUE
  attr(ans$d,"sparse") = TRUE
  attr(ans$v,"sparse") = TRUE
  ans
}

svd_scidb = function(x, nu=min(dim(x)), nv=nu)
{
  got_tsvd = length(grep("tsvd",.scidbenv$ops[,2]))>0
  if(missing(nu)) nu = min(dim(x))
  if(!is.sparse(x) && (nu > (min(dim(x))/3)) || !got_tsvd)
  {
# Compute the full SVD
    u = tmpnam()
    d = tmpnam()
    v = tmpnam()
    xend = as.numeric(scidb_coordinate_bounds(x)$length) - 1
    schema = sprintf("[%s=0:%s,1000,0,%s=0:%s,1000,0]",
                     dimensions(x)[1],noE(xend[1]),
                     dimensions(x)[2],noE(xend[2]))
    schema = sprintf("%s%s",build_attr_schema(x),schema)
    iquery(sprintf("store(gesvd(reshape(%s,%s),'left'),%s)",x@name,schema,u))
    iquery(sprintf("store(gesvd(reshape(%s,%s),'values'),%s)",x@name,schema,d))
    iquery(sprintf("store(transpose(gesvd(reshape(%s,%s),'right')),%s)",x@name,schema,v))
    ans = list(u=scidb(u,gc=TRUE),d=scidb(d,gc=TRUE),v=scidb(v,gc=TRUE))
    attr(ans$u,"sparse") = FALSE
    attr(ans$d,"sparse") = FALSE
    attr(ans$v,"sparse") = FALSE
    return(ans)
  }
  warning("Using the IRLBA truncated SVD algorithm")
  return(tsvd(x,nu))
}


# Miscellaneous functions
log_scidb = function(x, base=exp(1))
{
  w = scidb_types(x) == "double"
  if(!any(w)) stop("requires one double-precision valued attribute")
  if(class(x) %in% "scidb") attr = .get_attribute(x)
  else attr = x@attributes[which(w)[[1]]]
  new_attribute = sprintf("%s_log",attr)
  if(base==exp(1))
  {
    query = sprintf("apply(%s, %s, log(%s))",x@name, new_attribute, attr)
  } else if(base==10)
  {
    query = sprintf("apply(%s, %s, log10(%s))",x@name, new_attribute, attr)
  }
  else
  {
    query = sprintf("apply(%s, %s, log(%s)/log(%.15f))",x@name, new_attribute, attr, base)
  }
  .scidbeval(query,`eval`=FALSE)
}

# S4 method conforming to standard generic trig functions. See help for
# details about attribute selection and naming.
fn_scidb = function(x,fun,attr)
{
  if(missing(attr))
  {
    attr = x@attributes
  }
  new_attributes = make.unique_(x@attributes,attr)
  expr = paste(sprintf("%s, %s(%s)", new_attributes, fun, attr),collapse=",")
  query = sprintf("apply(%s, %s)",x@name, expr)
  query = sprintf("project(%s, %s)",query, paste(new_attributes,collapse=","))
  ren = paste(sprintf("%s,%s",new_attributes,attr),collapse=",")
  query = sprintf("attribute_rename(%s,%s)",query,ren)
  .scidbeval(query,`eval`=FALSE,gc=TRUE,depend=list(x),`data.frame`=(is.scidbdf(x)))
}

# S3 Method conforming to usual diff implementation. The `differences`
# argument is not supported here.
diff.scidb = function(x, lag=1, ...)
{
  y = lag(x,lag)
  n = make.unique_(c(x@attributes,y@attributes),"diff")
  z = merge(y,x,by=dimensions(x),all=FALSE)
  expr = paste(z@attributes,collapse=" - ")
  project(bind(z, n, expr), n)
}
