primes <-
function(start=12, end=9999){
dummy  <- 0
dummyr <- 0
first <- 0
nixda <- 0
test <- is.negative(start)
if(test==TRUE) {start <- start*(-1)}
test <- is.negative(end)
if(test==TRUE){end <- end*(-1)}
if(end<start) {
xyz <- start
start <- end
end <- xyz
}
if(start<18){
primzahlen <- c(1,2,3,5,7,11,13,17)
first <- 1
start <- 18
}
# befinde ich mich in der 3er-Reihe?
test <- is.even(start)
if(test==TRUE){start <- start+1}
test <- start/3
test2 <- floor(test)
if(test==test2){
dummy3 <-1
}else{
dummy3 <- 0
dummy <- 1
}
## Test auf Primzahl
while(start<(end+1)){
    test1 <- 5
    test2 <- ceiling(sqrt(start))
    while(test1 < (test2+1)){
    
    primtest1 <- start/test1
    primtest2 <- floor(primtest1)
    if(primtest1==primtest2){
    nixda <- nixda + 1
    }
 
    if(dummyr==0){
    test1 <- test1 + 2
    dummyr <- 1
    }else{
    test1 <- test1 + 4
    dummyr <- 0
    }
    }
    
        if(nixda==0){
        #cat("Primzahl gefunden: ",start,"\r")
        if(first==0){
    primzahlen <- start
    first <- 1
    }else{
    primzahlen <- c(primzahlen, start)
    }

        }
nixda <- 0
## 2er und 3er auslassen
#cat (c(start, "\r"))
if(dummy==1){
start <- start+4
dummy<-0
}else{
start <- start+2
dummy <- 1
}
if(dummy3==1){
dummy <- 0
dummy3 <- 0
}
}

return(primzahlen)
}

