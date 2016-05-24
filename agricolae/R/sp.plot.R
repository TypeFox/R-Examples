`sp.plot` <-
function(block, pplot, splot, Y)
{
name.y <- paste(deparse(substitute(Y)))
name.r <- paste(deparse(substitute(block)))
name.p <- paste(deparse(substitute(pplot)))
name.sp <- paste(deparse(substitute(splot)))
block<-as.factor(block)
pplot<-as.factor(pplot)
splot<-as.factor(splot)
cat("\nANALYSIS SPLIT PLOT: ", name.y, "\nClass level information\n\n")
nrep<- length(unique(block))
np  <- length(unique(pplot))
nsp <- length(unique(splot))
cat(name.p,  "\t: ",unique(as.character(pplot)),"\n")
cat(name.sp, "\t: ",unique(as.character(splot)),"\n")
cat(name.r,  "\t: ",unique(as.character(block)),"\n")
cat("\nNumber of observations: ", length(Y), "\n\n")
model<- aov(Y ~ block*pplot*splot)
B<-suppressWarnings(anova(model))
W<-NULL
W<-B[c(1,2,7,3,6,7),]
for (j in 1:2){
W[3,j]<-B[4,j]
W[6,j]<-B[5,j]+B[7,j]
}
W[,3]<-W[,2]/W[,1]
W[1:2,4]<-W[1:2,3]/W[3,3]
W[4:5,4]<-W[4:5,3]/W[6,3]
# Pvalue
W[1:2,5]<-1-pf(W[1:2,4],W[1:2,1],W[3,1])
W[4:5,5]<-1-pf(W[4:5,4],W[4:5,1],W[6,1])
N<-NULL
N[1]<- name.r
N[2]<- name.p
N[3]<- "Ea"
N[4]<- name.sp
N[5]<- paste(name.p,":",name.sp,sep="")
N[6]<- "Eb"
rownames(W)<-N
attributes(W)$heading[2]<-paste("Response:",name.y)
print(W)
medy <- mean(Y,na.rm=TRUE)
gl.a<-W[3,1]; Ea<-W[3,3]
gl.b<-W[6,1]; Eb<-W[6,3]
cat("\ncv(a) =",round(sqrt(Ea)*100/medy,1),"%,", 
"cv(b) =",round(sqrt(Eb)*100/medy,1),"%,",
"Mean =", medy,"\n\n")
output<-list(anva=W, gl.a=gl.a,gl.b=gl.b,Ea=Ea,Eb=Eb)
invisible(output)
}


