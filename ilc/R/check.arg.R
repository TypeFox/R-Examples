check.arg <-
function(arg, choices, what = "", n = length(arg), value = T){
    index <- T
    if(length(arg) > n) {
        warning(paste("Only the first", n, "values tested in",
                      deparse(substitute(arg)), "!\n"))
        arg <- arg[1:n]
    }
    while(index) {
        index <- pmatch(arg, choices, F, duplicates.ok = T)
        if(!all(index)) {
            cat("Sorry, it seems you have wrong argument(s)! in ",
                deparse(substitute(arg)), " ...\n")
            arg <- enter(deparse(substitute(arg)), choices, what = what, n = n)
            index <- 1
        }
        else break
    }
    if(value) choices[index]
    else index
}
