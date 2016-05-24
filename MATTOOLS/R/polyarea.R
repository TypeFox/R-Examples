 
polyarea <-  function(inpolygon)
  {
    a=inpolygon
    tlen=length(a[,1])
    xi=a[1:(tlen-1),1]
    xiPLUS1=a[2:tlen,1]
    yi=a[1:(tlen-1),2]
    yiPLUS1=a[2:tlen,2]
    parea=sum( (xi*yiPLUS1)-(yi*xiPLUS1)  )/2
    parea
  }