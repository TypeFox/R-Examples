## Function written to match MATLAB function
## Author: Andrew Hooker

fprintf <- function(file="",...){
  if(any(class(file)=="file")){
    cat(sprintf(...),file=file)
  } else {
    if(file==""){cat(sprintf(...))} else {
      cat(sprintf(file,...))
    }
  }
}
