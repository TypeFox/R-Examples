## Function written to match MATLAB function
## Author: Andrew Hooker

write_matrix <- function (f,x) 
{
  for(i in 1:size(x,1)){
    fprintf(f,"%6e",x[i,])
    fprintf(f,"\n")
  }
  return()
}
