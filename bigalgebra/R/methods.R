setMethod("%*%",signature(x="big.matrix", y="big.matrix"),
  function(x,y) 
  {
    dgemm(A=x, B=y)
  },
  valueClass="big.matrix"
)

setMethod("%*%",signature(x="matrix", y="big.matrix"),
  function(x,y)
  {
    if(dim(x)[2] != dim(y)[1]) stop("non-conformant matrices")
    R = options("bigalgebra.mixed_airthmetic_returns_R_matrix")[[1]]
    if(!is.null(R) && R) return(dgemm(A=x, B=y, C=0))
    dgemm(A=x, B=y)
  },
  valueClass="matrix"
)

setMethod("%*%",signature(x="big.matrix", y="matrix"),
  function(x,y) 
  {
    if(dim(x)[2] != dim(y)[1]) stop("non-conformant matrices")
    R = options("bigalgebra.mixed_airthmetic_returns_R_matrix")[[1]]
    if(!is.null(R) && R) return(dgemm(A=x, B=y, C=0))
    dgemm(A=x, B=y)
  },
  valueClass="matrix"
)

setMethod("Arith",c(e1="big.matrix", e2="big.matrix"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    switch(op,
      `+` = daxpy(1.0,e1,e2),
      `-` = daxpy(-1.0,e2,e1),
      stop("Undefined operation")
    )
  }
)
setMethod("Arith",c(e1="big.matrix", e2="matrix"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    switch(op,
      `+` = daxpy(1.0,e1,e2),
      `-` = daxpy(-1.0,e2,e1),
      stop("Undefined operation")
    )
  }
)
setMethod("Arith",c(e1="matrix", e2="big.matrix"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    switch(op,
      `+` = daxpy(1.0,e1,e2),
      `-` = daxpy(-1.0,e2,e1),
      stop("Undefined operation")
    )
  }
)
setMethod("Arith",c(e1="numeric", e2="big.matrix"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    if(length(e1)==1 && op=="*") return(daxpy(e1,e2))
    switch(op,
      `+` = daxpy(1.0,e1,e2),
      `-` = daxpy(-1.0,e2,e1),
      stop("Undefined operation")
    )
  }
)
setMethod("Arith",c(e1="big.matrix", e2="numeric"),
  function(e1,e2) 
  {
    op = .Generic[[1]]
    if(length(e2)==1 && op=="*") return(daxpy(e2,e1))
    switch(op,
      `+` = daxpy(1.0,e1,e2),
      `-` = daxpy(-1.0,e2,e1),
      stop("Undefined operation")
    )
  }
)
