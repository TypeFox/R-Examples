ry_bmf <-
function(y,l,d)
{
  .C("ry_bmf",y=as.double(y),l=as.double(l),d=as.double(d),
               n=as.integer(length(y)))$y
}
