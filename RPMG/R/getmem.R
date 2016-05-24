`getmem` <-
function(v, mem=1)
{
  ########    get the mem element of each vector in a  list of vectors
if(missing(mem)) { mem=1 }
return(v[mem])
}

