dms <-
function (dom,paired=FALSE) { #produces character matrix from dominance matrix with signs
dx <- matrix(nrow=nrow(dom),ncol=ncol(dom))
for (i in 1:nrow(dom)) {
  for (j in 1:ncol(dom))  {
     if(sign(dom[i,j])==-1) {dx[i,j]<-"-"}
     if(sign(dom[i,j])==1) {dx[i,j]<-"+"}
     if(sign(dom[i,j])==0) {dx[i,j]<-"O"} }}
 if ((paired==TRUE)&&(nrow(dom)==ncol(dom))){
   for (i in 1:nrow(dom)) {
      if(sign(dom[i,i])==-1) {dx[i,i]<-"<"}
      if(sign(dom[i,i])==1) {dx[i,i]<-">"}
      if(sign(dom[i,i])==0) {dx[i,i]<-"="} }} 
return(dx)}

