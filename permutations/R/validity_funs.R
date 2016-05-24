singleword_valid <- function(w){  # takes an integer vector
    if(length(w)==0){return(TRUE)}
    if(identical(seq_len(max(w)),sort(w))){
        return(TRUE)
    } else {
          stop("invalid word: rows should be permutations of seq_len(ncol(w))")
      }
}

cyclist_valid <- function(x){   # takes a cyclist and checks it for validity
    if(length(unlist(x))==0){return(TRUE)}
    x <- c(x,recursive=TRUE)
    if(any(x<0)){
        stop("negative elements")
    } else if (any(x==0)){
        stop("zero element")
    } else if (any(abs(x-round(x))>0.1)) {
        stop('non-integer entries')
    } else if (any(table(x)>1)){
        stop("repeated value")
    } else {
        return(TRUE)
    }
}
