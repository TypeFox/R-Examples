aspectAbs <- function(r,pow = 1) {
list(f=sum(abs(r)^pow)/2,g=pow*sign(r)*(abs(r)^(pow-1)))
}

