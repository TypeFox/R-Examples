plot.quint <-
function(x, digits=2,...){
  #x is an object of class "quint"
   x$si[,4] <- round(x$si[,4],digits=digits)
   party.quint <- as.party(x)
   plot(party.quint, inner_panel= node_quint, terminal_panel=terminal_quint, ...)
 return(invisible(party.quint))
}
