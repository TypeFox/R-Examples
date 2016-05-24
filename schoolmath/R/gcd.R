gcd <-
function(x,y){
gcd <- 1
liste <- NULL
a <- prime.factor(x)
b <- prime.factor(y)
a1 <- unique(a)
b1 <- unique(b)
c <- length(a1)
d <- length(b1)
if(c<d){
gemeinsam <- is.element(a1,b1)
switch <- a1
}else{
gemeinsam <- is.element(b1,a1)
switch <- b1
}
gem.l <- length(gemeinsam)
start <- 1
end <- gem.l+1
while(start<end){
if(gemeinsam[start]==TRUE){
liste <- c(liste, switch[start])
}
start <- start+1
}
liste.l <- length(liste)
start <-1
end <- liste.l+1
while(start<end){
dieser <- liste[start]
test1 <- sum(as.integer(is.element(a, dieser)))
test2 <- sum(as.integer(is.element(b, dieser)))
if(test1<test2){
zwischen <- liste[start]^test1
}else{
zwischen <- liste[start]^test2
}
gcd <- gcd*zwischen
start <- start+1
}
return(gcd)
}

