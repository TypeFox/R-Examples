# 
# Copyright (c) 2010, 2014, IBM Corp. All rights reserved. 
# 		
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>. 
# 

################ Stats ############################

setMethod("mean", signature(x="ida.data.frame"),
    function(x) {
      query <- paste(paste("AVG(\"", x@cols, "\") AS ", x@cols, collapse=",",sep=''),
          " FROM ", idadf.from(x))
      if (nchar(x@where)) {
        query <- paste(query, " WHERE ", x@where)
      }
      return(idaQuery ("SELECT ", query))
    }
)


setMethod("min", signature(x="ida.data.frame"),
    function(x) {
      query <- paste(paste("MIN(\"", x@cols, "\") AS ", x@cols, collapse=",",sep=''),
          " FROM ", idadf.from(x))
      if (nchar(x@where)) {
        query <- paste(query, " WHERE ", x@where)
      }
      return(idaQuery ("SELECT ", query))
    }
)


setMethod("max", signature(x="ida.data.frame"),
    function(x) {
      query <- paste(paste("MAX(\"", x@cols, "\") AS ", x@cols, collapse=",",sep=''),
          " FROM ", idadf.from(x))
      if (nchar(x@where)) {
        query <- paste(query, " WHERE ", x@where)
      }
      return(idaQuery ("SELECT ",query))
    }
)


################ hist ############################

hist.ida.data.frame <- function (x, breaks="Sturges",...) {
  if (!is.ida.data.frame(x))
    stop("x is not a proper ida.data.frame object")
  if (length(x@cols) > 1)
    stop("there can be only one column selected")
  
  # we can do smth with breaks if we want it run on SPUs
  # pass the rest of the arguments
  x <- as.data.frame(x)[[1]]
  invisible(hist(x, breaks=breaks,...))
}

setMethod("hist", signature(x="ida.data.frame"), hist.ida.data.frame)

################ Basic ############################

setMethod("dim", signature(x="ida.data.frame"),
    function(x) {
      rowsno <- idaScalarQuery("SELECT CAST(COUNT(*) AS INTEGER) FROM ",
          idadf.from(x), ifelse(nchar(x@where), paste(" WHERE ", x@where), ""))
      
      # Integers are 32bit in R, there are tables with more rows than that
      # In this case, we need to return double instead
      if(as.integer(rowsno) == as.numeric(rowsno)){	
        return(c(as.integer(rowsno), length(x@cols)))
      } else {
        return(c(as.double(rowsno), length(x@cols)))
      }	
    }
)

# ---------------------------------------------------------------------

setMethod("length", signature(x="ida.data.frame"),
    function(x) { return(length(x@cols)) }
)

# ---------------------------------------------------------------------

setMethod("NROW", signature(x="ida.data.frame"),
    function(x) { return(nrow(x)) }
)

# ---------------------------------------------------------------------

setMethod("NCOL", signature(x="ida.data.frame"),
    function(x) { return(ncol(x)) }
)

# ---------------------------------------------------------------------

setMethod("colnames", signature(x="ida.data.frame"),
    function(x) { x@cols }
)

# ---------------------------------------------------------------------
setMethod("head", signature(x="ida.data.frame"),
    function(x, n = 6, ...) {
      if (n >= 0) {
        ans <- idaQuery(idadf.query(x), " FETCH FIRST ", format(n, scientific = FALSE), " ROWS ONLY");
        
        if (nrow(ans) > 0) rownames(ans) <- 1:nrow(ans);
        return(ans)
      } else {
        nr <- nrow(x)
        n <- abs(n)
        ans <- idaQuery(idadf.query(x), " FETCH FIRST ", format(nr - n, scientific = FALSE), " ROWS ONLY")
        
        if ((nr-n) != 0) rownames(ans) <- 1:(nr-n);
        return(ans)
      }
    }
)


# ---------------------------------------------------------------------

setMethod("names", signature(x="ida.data.frame"),
    function(x) { return(x@cols) }
)
