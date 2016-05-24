scm <-
function(x,y){
scm <- 1
liste <- NULL
a <- prime.factor(x)
b <- prime.factor(y)
aa <-a
a1 <- unique(a)
bb <-b
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
aa <- aa[aa!=dieser]
bb <- bb[bb!=dieser]
if(test1>test2){
zwischen <- liste[start]^test1
}else{
zwischen <- liste[start]^test2
}
scm <- scm*zwischen
start <- start+1
}
rest <- c(aa,bb)
rest.ende <- length(rest)+1
start<-1
while(start<rest.ende){
scm <- scm*rest[start]
start <- start+1
}
return(scm)
}

