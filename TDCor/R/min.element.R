min.element <-
function(u,v)

{

k=ceiling(sign(u-v)/2+0.1)

return(u*(1-k)+v*k)

}
