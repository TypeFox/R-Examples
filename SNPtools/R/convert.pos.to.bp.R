convert.pos.to.bp <-
function(pos) {

  if(max(pos) <= 200) {
    pos = pos * 1e6
  } # if(max(pos] <= 200))

  return(pos)
}
