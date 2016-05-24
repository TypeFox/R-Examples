a.glc <-
function(side=1,  line=0){

#########################################
# with vector
if(length(side)>1 | length(line)>1){
if(length(side)==1)  side<-rep(side, length(line))
if(length(line)==1)  line<-rep(line, length(side))

at<-side
for(i in 1:length(side)){
at[i] <- a.glc(side[i], line[i])
}
return(at)
}
#########################################

if(length(side)==1 & length(line)==1){
if(side==0){
at <- par('usr')[1]  +  (par('usr')[2]-par('usr')[1])/2
}
if(side==1){
one.line<- yinch(y = 0.2, warn.log = TRUE)
at <- par('usr')[3]  -  one.line * line
}
if(side==2){
one.line<- xinch(x = 0.2, warn.log = TRUE)
at <- par('usr')[1]  -  one.line * line
}
if(side==3){
one.line<- yinch(y = 0.2, warn.log = TRUE)
at <- par('usr')[4]  +  one.line * line
}
if(side==4){
one.line<- xinch(x = 0.2, warn.log = TRUE)
at <- par('usr')[2]  +  one.line * line
}
if(side==5){
at <- par('usr')[3]  +  (par('usr')[4]-par('usr')[3])/2
}

return(at)
}
}
