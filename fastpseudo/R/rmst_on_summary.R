rmst_on_summary <- function(df, tmax) {
  not_dead <- mapply(
    function(eligible, event) {
      if (eligible == 0) {
        return(1)
      } else {
        return((eligible - event) / eligible)
      }
    },
    df$eligible[df$time < tmax],
    df$event[df$time < tmax]
  )
  
  return(sum(cumprod(not_dead)))
}