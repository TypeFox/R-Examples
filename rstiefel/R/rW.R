rW <-
function(kap,m)
{
  #simulate W as described in Wood(1994)
  .C("rW",kap=as.double(kap),m=as.integer(m),w=double(1))$w
}
