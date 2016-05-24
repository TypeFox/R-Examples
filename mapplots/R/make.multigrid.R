make.multigrid <-
function(x,y,z,group,...){
  grd <- NULL
  for(g in unique(as.character(group))){
    i <- which(group==g)
    grd[[g]] <- make.grid(x[i],y[i],z[i],...)
  }
  return(grd)
}

