"getAlloc" <-
function(len,lmin=2){

#### this makes all alloc matrices from lmin up to len
  
write("",file="tmp-alloc.R")

## lines to declare alloc matrices
for(i in lmin:len){
 text <- paste("alloc",i," <- array(0,c(1,",i,"))",sep="")
 write(text,file="tmp-alloc.R",append=T)
}

## lines for i=1
text1 <- "for(n1 in 0:1){"
text2 <- "ndum1 <- n1"
text3 <- "nlast2 <- 2 - ndum1"
if(lmin==2) text4 <- "alloc2 <- rbind(alloc2,c(n1,nlast2))"
write(text1,file="tmp-alloc.R",append=T)
write(text2,file="tmp-alloc.R",append=T)
write(text3,file="tmp-alloc.R",append=T)
if(lmin==2) write(text4,file="tmp-alloc.R",append=T)

if(len>2){
## loop to write lines for i=2,len-1
for(i in 2:(len-1)){
 text1 <- paste("for(n",i," in 0:(",i,"-ndum",i-1,")){",sep="")
 text2 <- paste("ndum",i," <- ndum",i-1," + n",i,sep="")
 text3 <- paste("nlast",i+1," <- ",i+1," - ndum",i,sep="")
 write(text1,file="tmp-alloc.R",append=T)
 write(text2,file="tmp-alloc.R",append=T)
 write(text3,file="tmp-alloc.R",append=T)

 if(i>=lmin-1){
   text <- paste("n",1:i,sep="",collapse=",")
   text4 <- paste("alloc",i+1," <- rbind(alloc",
                  i+1,",c(",text,",nlast",i+1,"))",sep="")
   write(text4,file="tmp-alloc.R",append=T)
 }
}
}

## writes close of loops
text <- paste(rep("}",len-1),collapse="")
write(text,file="tmp-alloc.R",append=T)

## lines to remove first row of alloc matrices
for(i in lmin:len){
 text <- paste("alloc",i," <- alloc",i,"[-1,]",sep="")
 write(text,file="tmp-alloc.R",append=T)
}

source("tmp-alloc.R")

}

