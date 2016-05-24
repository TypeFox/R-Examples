`order.stat` <-
function(treatment,means,minimum,console) {
n<-length(means)
z<-data.frame(treatment,means)
letras<-c(letters[1:26],LETTERS[1:26])
# ordena tratamientos
w<-z[order(z[,2],decreasing=TRUE),]
M<-rep("",n)
k<-1
k1<-0
j<-1
i<-1
cambio<-n
cambio1<-0
chequeo=0
M[1]<-letras[k]
while(j<n) {
chequeo<-chequeo+1
if (chequeo > n) break
for(i in j:n) {
s<-abs(w[i,2]-w[j,2])<=minimum
if(s) {
if(lastC(M[i]) != letras[k])M[i]<-paste(M[i],letras[k],sep="")
}
else {
k<-k+1
cambio<-i
cambio1<-0
ja<-j
for(jj in cambio:n) M[jj]<-paste(M[jj],"",sep="") # El espacio
M[cambio]<-paste(M[cambio],letras[k],sep="")
for( v in ja:cambio) {
if(abs(w[v,2]-w[cambio,2])>minimum) {j<-j+1
cambio1<-1
}
else break
}
break
}
}
if (cambio1 ==0 )j<-j+1
}
#-----------
w<-data.frame(w,stat=M)
trt<-as.character(w$treatment)
means<-as.numeric(w$means)
for(i in 1:n){
if(console)cat(M[i],"\t",trt[i],"\t",means[i],"\n")
}
output<-data.frame(trt,means,M)
return(output)
}

