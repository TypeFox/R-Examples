`circle` <-
function(n=1)
  {
    if(missing(n)) n<-1
    ##%   create a circle for plotting
    
    i=pi*seq(from=0, to=360, by=n)/180
    
    cx = cos(i);
    cy = sin(i);
    
    return(list(x=cx, y=cy))
    
}

