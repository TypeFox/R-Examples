decimal2fraction <-
function(decimal, period=0){
dummy<-0
ganze <- floor(decimal)
anzahl.ganze <- nchar(ganze)
nuller <- nchar(decimal) - (anzahl.ganze+1) # wird abgezogen
nenner1 <- 10^nuller
zaehler1 <- floor(decimal*nenner1)
if(period<0) period <- period*(-1)
if(period>0){
test <- is.whole(period)
if(test==FALSE){
msg <- "period must be a whole number\r"
return(msg)
}
dummy <-1
zaehler2 <- period
nenner2 <- nenner1*9

zaehler1 <- zaehler1*9
nenner1 <- nenner1*9

zaehler3 <- zaehler1 + zaehler2
nenner3 <- nenner2
}

if(dummy==0) {
zaehler3 <- zaehler1
nenner3 <- nenner1
}
test1 <- is.prim(zaehler3)
if(test1==FALSE){
ggT <- gcd(zaehler3,nenner3)
zaehler3 <- zaehler3/ggT
nenner3 <- nenner3/ggT

}

cat("\r--------------\r")
cat(zaehler3,"/",nenner3)

}

