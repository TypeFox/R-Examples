qlogunif <- function(p,a=0,b=1,base=10) {
        -log(1-qunif(p,a,b),base)
}
p <- qqmath(~-log10(x) | group, data=someData, distribution=qlogunif)
print(p)
