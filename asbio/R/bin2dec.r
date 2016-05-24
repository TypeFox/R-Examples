bin2dec <- function(digits, round = 4){
x <- as.character(as.numeric(digits))
if(as.numeric(x)!= 0 & as.numeric(x) < 1){
h <- unlist(strsplit(x, ""))
b <- c(0, as.numeric(unlist(h[-c(1, 2)])))
if(any(b != "0" & b != "1")) stop("Digits must be binary")
pow <- -1*(1:length(b))[-(length(b))]
res <- round(sum(c(0,2^pow)*b),round)
}
if(as.numeric(x) == 0 | as.numeric(x) >= 1){
b <- as.numeric(unlist(strsplit(x, "")))
if(any(b != "0" & b != "1")) stop("Digits must be binary")
pow <- (length(b) - 1):0
res <- round(sum((2^pow)*b),round)
}
res
}

