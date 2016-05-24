## Helper functions for Mirror operations with respect to a given Cartesian plan
## or a plan defined by two lattice vectors.

# With respect to the xy-plan
Mxy <- function(...)
  UseMethod("Mxy")

Mxy.coords <- function(x, mask = TRUE, cryst1 = NULL, ...)
  mirror(x, c(0,0,0), c(1,0,0), c(0,1,0), mask=mask, cryst1=cryst1)

Mxy.pdb <- function(x, mask = TRUE, cryst1 = x$cryst1, ...)
  mirror(x, c(0,0,0), c(1,0,0), c(0,1,0), mask=mask, cryst1=cryst1)

# With respect to the yz-plan
Myz <- function(...)
  UseMethod("Myz")

Myz.coords <- function(x, mask = TRUE, cryst1 = NULL, ...)
  mirror(x, c(0,0,0), c(0,1,0), c(0,0,1), mask=mask, cryst1=cryst1)

Myz.pdb <- function(x, mask = TRUE, cryst1 = x$cryst1, ...)
  mirror(x, c(0,0,0), c(0,1,0), c(0,0,1), mask=mask, cryst1=cryst1)

# With respect to the zx-plan
Mzx <- function(...)
  UseMethod("Mzx")

Mzx.coords <- function(x, mask = TRUE, cryst1 = NULL, ...)
  mirror(x, c(0,0,0), c(0,0,1), c(1,0,0), mask=mask, cryst1=cryst1)

Mzx.pdb <- function(x, mask = TRUE, cryst1 = x$cryst1, ...)
  mirror(x, c(0,0,0), c(0,0,1), c(1,0,0), mask=mask, cryst1=cryst1)


# With respect to the ab-plan
Mab <- function(...)
  UseMethod("Mab")

Mab.coords <- function(x, cryst1, mask = TRUE, ...){
  if(missing(cryst1))
    stop("'cryst1' is required to defined the mirror plan")
  cell <- cell.coords(x, cryst1=cryst1)
  mirror(x, c(0,0,0), cell[,"a"], cell[,"b"], mask=mask, cryst1=cryst1)
}
  
Mab.pdb <- function(x, cryst1 = x$cryst1, mask = TRUE, ...){
  if(is.null(cryst1))
    stop("'cryst1' is required to defined the mirror plan")
  cell <- cell.coords(x, cryst1=cryst1)
  mirror(x, c(0,0,0), cell[,"a"], cell[,"b"], mask=mask, cryst1=cryst1)
}

# With respect to the bc-plan
Mbc <- function(...)
  UseMethod("Mbc")

Mbc.coords <- function(x, cryst1, mask = TRUE, ...){
  if(missing(cryst1))
    stop("'cryst1' is required to defined the mirror plan")
  cell <- cell.coords(x, cryst1=cryst1)
  mirror(x, c(0,0,0), cell[,"b"], cell[,"c"], mask=mask, cryst1=cryst1)
}

Mbc.pdb <- function(x, cryst1 = x$cryst1, mask = TRUE, ...){
  if(is.null(cryst1))
    stop("'cryst1' is required to defined the mirror plan")
  cell <- cell.coords(x, cryst1=cryst1)
  mirror(x, c(0,0,0), cell[,"b"], cell[,"c"], mask=mask, cryst1=cryst1)
}

# With respect to the ca-plan
Mca <- function(...)
  UseMethod("Mca")

Mca.coords <- function(x, cryst1, mask = TRUE, ...){
  if(missing(cryst1))
    stop("'cryst1' is required to defined the mirror plan")
  cell <- cell.coords(x, cryst1=cryst1)
  mirror(x, c(0,0,0), cell[,"c"], cell[,"a"], mask=mask, cryst1=cryst1)
}

Mca.pdb <- function(x, cryst1 = x$cryst1, mask = TRUE, ...){
  if(is.null(cryst1))
    stop("'cryst1' is required to defined the mirror plan")
  cell <- cell.coords(x, cryst1=cryst1)
  mirror(x, c(0,0,0), cell[,"c"], cell[,"a"], mask=mask, cryst1=cryst1)
}
