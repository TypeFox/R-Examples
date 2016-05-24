fac.match <- function(x, table, col.names, nomatch = NA_integer_, multiples.allow = FALSE)
#Function to find, for each combination of col.names in x, the row that has the same combination in table
#It can be viewed as a generalization of the match function from a single vector to multiple vectors
{ if (class(col.names) != "character") 
    stop("Must supply a character vector of column names")
  ncols <- length(col.names)
  if (any(!(col.names %in% names(x))))
    stop("All column names must be in x")
  if (any(!(col.names %in% names(table))))
    stop("All column names must be in x")
  if (length(dim(x)) != 2 | length(dim(table)) != 2)
    stop("Both x and table must have two dimensions")
  col.list <- as.list(x[col.names])
  index <- as.data.frame.table(with(x, 
                    by(x, 
                       col.list, 
                       function(x, table, col.names)
                       {  k <- rep(TRUE, nrow(table))
                          i <- 1
                          for (i in 1:ncols)
                          { k <- k & (table[[col.names[i]]] == x[[col.names[i]]])
                          }
                          if (anyNA(k))
                            k[is.na(k)] <- FALSE
                          if (sum(k) > 1)
                          { if (!multiples.allow)
                            { print(table[k, col.names])
                              stop("The above combination in x has ",
                                   "more than one combination in table")
                            }
                            else
                            { warning("The following combination in x has ",
                                      "more than one combination in table\n",
                                                "- only the first matched")
                              print(table[k, col.names])
                              nk1 <- sum(k)-1
                              k[k] <- c(TRUE, rep(FALSE, nk1))
                            }
                          }
                          if (all(!k))
                            k <- NA
                          else
                            k <- c(1:nrow(table))[k]
                          return(k)
                         },
                        table = table, col.names = col.names)))
  index <- na.omit(index)
  index <- merge(x[col.names], index, all.x=TRUE, sort=FALSE)$Freq
  if (anyNA(index))
    index[is.na(index)] <- nomatch
  index <- as.integer(index)
  return(index)
}
