a.alpha.ade <-
function(color, alpha=1){

#########################################
# make single values
if(length(color)>1 & length(alpha)==1){
y<-color
for(i in 1:length(color)){
y[i] <- a.alpha.ade(color[i], alpha=alpha)
}
return(y)
}
#########################################

#########################################
# make single values
if(length(color)>1 & length(alpha)==length(color)){
y<-color
for(i in 1:length(color)){
y[i] <- a.alpha.ade(color[i], alpha=alpha[i])
}
return(y)
}
#########################################

if(length(color)==1){
rgba<-grDevices::col2rgb(color, alpha = TRUE)
color<-rgb(rgba[1],rgba[2],rgba[3],names = NULL, alpha*255 , maxColorValue = 255)
return(color)
}
}
