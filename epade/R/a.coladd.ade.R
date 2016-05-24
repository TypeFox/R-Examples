a.coladd.ade <-
function(color, add=0){

#########################################
# make single values
if(length(color)>1){
y<-color
for(i in 1:length(color)){
y[i] <- a.coladd.ade(color[i], add=add)
}
return(y)
}
#########################################


if(length(color)==1){
rgba<-col2rgb(color, alpha = T)
rgba[1]=rgba[1]+add
rgba[2]=rgba[2]+add
rgba[3]=rgba[3]+add
if(rgba[1]<0)   rgba[1]<-0
if(rgba[2]<0)   rgba[2]<-0
if(rgba[3]<0)   rgba[3]<-0
if(rgba[1]>255)   rgba[1]<-255
if(rgba[2]>255)   rgba[2]<-255
if(rgba[3]>255)   rgba[3]<-255

color<-rgb(rgba[1],rgba[2],rgba[3],names = NULL,rgba[4],maxColorValue = 255)

return(color)
}
}
