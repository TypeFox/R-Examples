newSobject <-
function(meta.data) {
  special.indexes.1 <- which(colnames(meta.data) %in% c('$C','$I'))
  special.indexes.2 <- which(colnames(meta.data) %in% c('$H','$S'))
  special.indexes.all <- sort(c(special.indexes.1, special.indexes.2))
  index.colnames <- colnames(meta.data)[special.indexes.all]
  
  sym.var.length <- integer()
  for (i in 1:length(index.colnames)) {
    switch (index.colnames[[i]],
            '$C' = {
              sym.var.length <- c(sym.var.length, 1)
            },
            '$I' = {
              sym.var.length <- c(sym.var.length, 2)
            },
            '$H' = {
              sym.var.length <- c(sym.var.length, meta.data[[1, special.indexes.all[[i]] + 1]])
            },
            '$S' = {
              sym.var.length <- c(sym.var.length, meta.data[[1, special.indexes.all[[i]] + 1]])
            },
            stop("Invalid argument!"))
  }
  
  return (list(N = nrow(meta.data), M = length(special.indexes.all),
               sym.obj.names = row.names(meta.data), 
               sym.var.names = colnames(meta.data)[special.indexes.all + 1],
               sym.var.types = colnames(meta.data)[special.indexes.all],
               sym.var.length = sym.var.length,
               sym.var.starts = sort(c(special.indexes.1 + 1, special.indexes.2 + 2)),
               meta = meta.data, data = meta.data[, -c(special.indexes.all, special.indexes.2 + 1)]))
}
