merge_order <-
function(x, y, ..., sort = T, keep_order)
  {
    add.id.column.to.data <- function(DATA)
      {
        data.frame(DATA, id... = seq_len(nrow(DATA)))
      }
    order.by.id...and.remove.it <- function(DATA)
      {
        if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
        ss_r <- order(DATA$id...)
        ss_c <- colnames(DATA) != "id..."
        DATA[ss_r, ss_c]
      }
    if(!missing(keep_order))
      {
        if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
        if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
        warning("The function .merge_order only accepts NULL/1/2 values for the keep_order variable")
      }
    else {return(merge(x=x, y=y, ..., sort = sort))}
  }
