`meshgrid` <-
function(a,b) {
  return(list( x=outer(b*0,a,FUN="+"), y=outer(b,a*0,FUN="+") ))
}

