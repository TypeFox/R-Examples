is.prim <-
function(y){
starten <- 1
enden <- length(y)+1
while(starten<enden){
x <- y[starten]

error <- 0

test <- is.negative(x)
if(test==TRUE){
x <- x*(-1)
}

#test ob ganze Zahl
test <- is.whole(x)
if(test==FALSE) {
#cat("keine ganze Zahl\r")
error <- error+1
}

# 0 ist keine Primzahl
if(x==0){
##cat("ist 0 \r")
error <- error+1
}

#test ob gerade Zahle
test <- is.even(x)
if(test==TRUE & x!=2){
#cat("ist gerade\r")
error <- error+1
}

# test ob Quadratwurzel
test <- is.whole(sqrt(x))
if(test==TRUE & x!=1){
#cat("ist wurzel\r")
error <- error+1
}

#### los gehtz
if(error>0){
#cat("setze FALSE\r")
if (starten==1){
result=FALSE
}else{
result <- c(result, FALSE)
}
}else{
##### Teste, ob Primzahl
if(x==1 |x==2|x==3|x==5|x==7){
#cat("ist 2, setze TRUE\r")
if (starten==1){
result=TRUE
}else{
result <- c(result, TRUE)
}
}else{
anfang <- 3
ende <- ceiling(sqrt(x))
#cat("potentiell\r")
while(anfang<ende){

test1 <- x/anfang
test2 <- floor(test1)

if(test1==test2){
error <- error+1
}

anfang <- anfang+2
}
    if(error==0){
    if (starten==1){
result=TRUE
}else{
result <- c(result, TRUE)
}
    } else{
    if (starten==1){
result=FALSE
}else{
result <- c(result, FALSE)
}
    }
   }

}

starten <- starten+1
}
return(result)
}

