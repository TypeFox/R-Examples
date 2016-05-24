`strip.plot` <-
function (BLOCK, COL, ROW, Y)
{
name.y <- paste(deparse(substitute(Y)))
name.col <- paste(deparse(substitute(COL)))
name.row <- paste(deparse(substitute(ROW)))
name.r <- paste(deparse(substitute(BLOCK)))
cat("\nANALYSIS STRIP PLOT: ", name.y, "\nClass level information\n\n")
COL <- as.factor(COL)
ROW <- as.factor(ROW)
BLOCK <- as.factor(BLOCK)
nrep <- length(unique(BLOCK))
nCOL <- length(unique(COL))
nROW <- length(unique(ROW))
data <-data.frame(BLOCK, COL, ROW, Y)
names(data)<- c(name.r ,name.col,name.row, name.y)
cat(name.col,  "\t: ",unique(as.character(COL)),"\n")
cat(name.row,  "\t: ",unique(as.character(ROW)),"\n")
cat(name.r,    "\t: ",unique(as.character(BLOCK)),"\n")
cat("\nNumber of observations: ", length(Y), "\n\n")
model <-aov(Y~BLOCK+COL+  BLOCK:COL + ROW + BLOCK:ROW+ COL:ROW)
mm <- anova(model)
nn <- mm[3, ]
nn1<- row.names(mm)[3]
nn2<- row.names(mm)[4]
row.names(mm)[4] <- " "
mm[3, ] <- mm[4, ]
mm[4, ] <- nn
row.names(mm)[3] <- nn2
row.names(mm)[4] <- nn1
mm[2, 4] <- mm[2, 3]/mm[3, 3]
mm[2, 5] <- 1 - pf(mm[2, 4], mm[2, 1], mm[3, 1])
mm[4, 4] <- mm[4, 3]/mm[5, 3]
mm[4, 5] <- 1 - pf(mm[4, 4], mm[4, 1], mm[5, 1])
N<-NULL
N[1]<- name.r
N[2]<- name.col
N[3]<- "Ea"
N[4]<- name.row
N[5]<- "Eb"
N[6]<- paste(name.row,":",name.col,sep="")
N[7]<- "Ec"
NN<-paste(name.y,"~",N[1],"+",N[2], "+ Ea +",N[4],"+ Eb +",N[6],"+ Ec")
cat("model Y:",NN,"\n\n")
rownames(mm)<-N
attributes(mm)$heading[2]<-paste("Response:",name.y)
print(mm)
DFE <- df.residual(model)
MSE <- deviance(model)/DFE
medy <- mean(Y,na.rm=TRUE)
gl.a<-mm[3,1]; Ea<-mm[3,3]
gl.b<-mm[5,1]; Eb<-mm[5,3]
gl.c<-mm[7,1]; Ec<-mm[7,3]
cat("\ncv(a) =",round(sqrt(Ea)*100/medy,1),"%,", 
"cv(b) =",round(sqrt(Eb)*100/medy,1),"%,",
"cv(c) =",round(sqrt(Ec)*100/medy,1),"%,","Mean =", medy,"\n\n")
output<-list(data=data,anva=mm, gl.a=gl.a,gl.b=gl.b,gl.c=gl.c,Ea=Ea,Eb=Eb,Ec=Ec)
invisible(output)
}

