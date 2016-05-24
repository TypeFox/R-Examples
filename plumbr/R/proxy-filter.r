#' Proxy: filtering
#' @noRd
filter_proxy <- function(mf, i = NULL, j = NULL, rn) {
  if (is.null(i)) {
    varlist <- proxy_bindings(mf, j)
  } else {
    varlist <- filter_bindings(mf, j, i)
  }

  child <- .mutaframe(varlist, rn)

  if (is.null(i)) {
    add_listener(mf, function(event_i, event_j) {
      if (shape_changed(event_i, event_j)) {
        notify_listeners(child, NULL, NULL)
      } else {
        incl <- event_j %in% j
        notify_listeners(child, event_i[incl], event_j[incl])
      }
      
    })
  } else {    
    add_listener(mf, function(event_i, event_j) {
      if (shape_changed(event_i, event_j)) {
        notify_listeners(child, NULL, NULL)
      } else {
        new_i <- match(event_i, i)      
        incl <- !is.na(new_i)
      
        notify_listeners(child, new_i[incl], event_j[incl])        
      }
    })
  }
  
  child
}
