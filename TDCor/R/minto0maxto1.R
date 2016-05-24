minto0maxto1 <-
function(x)

{

x[sign(x)==-1]=0

x[sign(x-1)==1]=1

return(x)

}
