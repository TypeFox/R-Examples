"lines.gp" <-
function(x,...){
  # add line plot for one-dimensional gp
  if(x$d!=1){
    stop(" lines function only works for one-dimensional processes")
  }
  lines(getgrid(x),predict(x),...)
}
