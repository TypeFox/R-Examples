`screens` <-
function(n)
{
 ##X##   set up n screens for plotting in R
  ##X##  
 devl = dev.list()
 j = n-length(devl)
 if(j>0)
   {
for(i in 1:j)
  {
  ##  get(getOption("device"))()
   dev.new()
  }
}
}

