xyz2atom <-
function(xyz.ind) {
  return( unique( ceiling(xyz.ind/3) ) )
}
