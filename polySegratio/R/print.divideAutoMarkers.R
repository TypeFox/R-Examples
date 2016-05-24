`print.divideAutoMarkers` <-
function(x, ..., row.index=c(1:10), col.index=c(1:10),
         tabulate.extras=FALSE )
{
  cat(x$description,"\n\n")

  index <- c(1:length(names(x)))[names(x) %in% c("p10","p01","p11")]

  for (I in index){

    cat("**** data set:",x[[I]]$description,"\n")
    dim.m <- dim(x[[I]]$markers)
    cat("Dimension of marker data:", dim.m ,"\nData:\n")
    rows <- intersect(row.index,1:dim.m[1])
    cols <- intersect(col.index,1:dim.m[2])
    if (length(rows)> 0 & length(cols)> 0) {
      print(cbind(x[[I]]$markers[rows, cols],
                  cbind(r=x[[I]]$seg.ratios$r[rows],
                        n=x[[I]]$seg.ratios$n[rows],
                        ratio=x[[I]]$seg.ratios$seg.ratio[rows])),...)
    }
    cat("No. markers")
    print(table(x[[I]]$markers))
    if (length(x[[I]]$extras)>0) {
      if (tabulate.extras)
        print(table(x[[I]]$extras), ...)
      else
        if (length(rows)>0) {
          print(x[[I]]$extras[rows], ...)
        }
    }
  }
  
  cat("Call:\n")
  print(x$call)
  
}

