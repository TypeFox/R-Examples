writeData <- function(data, opt) {
  for(i in 1:length(data)) {
    write.table(data[[i]]@psi.df, file=paste(opt@makeps, "_dataset_",
                                    i, ".txt", sep=""),
                quote = FALSE,
                row.names = data[[i]]@x, col.names = data[[i]]@x2)
  }

}
