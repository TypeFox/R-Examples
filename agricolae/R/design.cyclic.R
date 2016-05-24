`design.cyclic` <-
function (trt, k, r, serie = 2, rowcol=FALSE, seed = 0, kinds = "Super-Duper",randomization=TRUE)
{
number<-10
if(serie>0) number<-10^serie
name.trt <- deparse(substitute(trt))
ntr <- length(trt)
# Control
if ((k*(r%/%k) != r )| (ntr >30) | (ntr<6) | (k > 10) | (r > 10) | (k==1) |(r==1))
{
cat("\nsee help(design.cyclic\n")
return()
}
if (seed == 0) {
genera<-runif(1)
seed <-.Random.seed[3]
}
set.seed(seed,kinds)
parameters<-list(design="cyclic",trt=trt,k=k,r=r,serie=serie,rowcol=rowcol,seed=seed,kinds=kinds)

ultimo<-rbind(rep(2,10),c(3,rep(4,6),6,5,5), c(4,3,3,3,3,6,6,3,7,3),
c(6,rep(5,4),3,3,5,4,8), c(5,6,6,6,6,5,5,4,6,6), c(4,4,4,4,rep(5,6)),
c(2,4,8,7,8,8,6,8,8,9), c(4,4,5,6,4,4,7,5,7,6), c(3,7,8,8,7,8,8,10,8,8), 
c(3,7,8,9,3,rep(7,5)), c(6,6,rep(8,6),10,11), c(6,6,8,9,10,9,rep(0,4)), 
c(rep(0,6),8,12,13,11),c(0,6,6,6,6,rep(11,5)), c(0,0,6,6,rep(10,5),11),
c(0,0,0,6,rep(10,4),12,12),c(rep(0,4),9,9,9,11,11,11), c(rep(0,5),8,8,8,13,13))
ultimo1<-rbind(rep(2,15),c(5,5,6,6,9,9,6,6,11,11,11,12,9,7,6),
c(8,8,4,9,6,4,9,9,7,4,5,8,6,10,9),c(3,3,8,3,3,7,11,12,4,9,8,4,12,14,13),
c(7,6,10,5,10,5,4,3,9,7,3,6,3,3,15),c(5,5,5,9,7,rep(6,9),7),
c(8,9,8,6,10,rep(11,5),12,11,12,13,15),c(9,8,9,9,rep(0,11)),
c(rep(0,4),9,10,10,10,10,13,11,15,14,12,11),c(8,8,8,9,9,10,10,18,11,rep(0,6)),
c(rep(0,9),rep(10,6)),c(9,11,11,12,12,13,13,14,14,rep(0,6)),
c(rep(0,9),17,14,17,14,17,20),c(4,4,rep(12,10),18,21,19),
c(10,10,rep(0,5),14,15,rep(0,4),16,0),c(0,0,rep(17,5),0,0,20,22,13,27,0,20),
c(0,16,17,0,0,12,12,12,15,16,0,0,0,18,20),c(12,0,0,14,14,rep(0,5),14,14,14,0,0),
c(6,6,6,7,rep(0,11)),c(rep(0,4),17,17,0,0,0,17,0,0,17,0,0),
c(rep(0,6),18,15,21,0,22,23,0,25,24),c(rep(16,9),rep(0,6)),
c(rep(0,9),rep(22,6)),c(13,rep(0,5),13,20,23,20,19,18,20,20,21),
c(0,16,17,18,19,20,rep(0,9)),c(15,15,15,15,15,15,rep(0,9)),
c(rep(0,6),21,22,rep(12,4),21,28,12))
repite <- c(2,4,6,8,10,3,6,9,4,8,5,10,6,7,8,9,10)
names(repite)<-1:17
kr<-2:10
pk <- c(1,6,9,11,14:19)
pk1<- c(1,6,10,14,17,19,22,24,26)
primeros <- list(1,1,1,1,1,c(1,2),c(1,3),c(1,2),c(1,2,4),c(1,2,5),c(1,2,3,5),
rbind(c(1,3,4,5),c(1,3,4,7)),c(1,2,3,4,7), c(1,2,3,4,5,8), c(1,2,3,4,5,7,9),
c(1,2,3,4,5,6,8,10),c(1,2,3,4,5,6,7,10,11))
primeros1 <- list(1,1,1,1,1, c(1,2),c(1,3),rbind(c(1,3),c(1,4)),rbind(c(1,2,4),c(1,2,7)),
rbind(c(1,2,6),c(1,3,13)),c(1,3,8,9),rbind(c(1,3,6,7),c(1,2,5,15)),rbind(c(1,2,3,7,10),c(1,3,9,10,13)),
rbind(c(1,2,3,4,9,13),c(1,2,3,5,9,12),c(1,3,8,9,11,12)),rbind(c(1,2,3,4,7,9,12),c(1,2,3,5,9,12,17)),
rbind(c(1,2,3,4,7,9,12,16),c(1,3,4,5,6,9,12,13)),rbind(c(1,2,4,5,8,9,10,11,13),c(1,2,3,4,7,9,13,16,20)))
#--------------------
B<-cbind(c(rep(2,5),3,3,3,4,4,5,5,6,7,8,9,10),c(2,4,6,8,10,3,6,9,4,8,5,10,6,7,8,9,10),
c(rep(1,11),2,rep(1,5)),c(rep(1,7),2,2,2,1,rep(2,6)))
nj <- r/k
r0<-0
inicial<-NULL
for (i in 1:nj) {
r0<-r0+k
if (ntr < 16 ) {
yp<-which(c(B[,1]==k & B[,2]==r0))
xp<-primeros[[yp]]
i1<-sum(B[1:yp,3])
if(!is.matrix(xp)) final <- ultimo[i1,ntr-5]
if( is.matrix(xp))
{i0<-sum(B[1:(yp-1),3])+1
m<-ultimo[i0:i1,ntr-5]
r1<-which(m !=0)
final<-m[r1]
xp<-xp[r1,] 
}
}

if ((ntr >= 16) & (ntr <=30) ) {
yp<-which(c(B[,1]==k & B[,2]==r0))
xp<-primeros1[[yp]]
i1<-sum(B[1:yp,4])
if(!is.matrix(xp)) final <- ultimo1[i1,ntr-15]
if( is.matrix(xp))
{i0<-sum(B[1:(yp-1),4])+1
m<-ultimo1[i0:i1,ntr-15]
r1<-which(m !=0)
final<-m[r1]
xp<-xp[r1,] 
}
}
inicial <- rbind(inicial,c(xp,final))
}

###############################
design<- NULL
BOOK<-NULL
for (irep in 1:nj) {
if(nj == 1)inicio <- inicial
if(nj>1)  inicio <- inicial[irep,]
###############################  
    c1 <- rep(0,ntr*k)
    dim(c1)<-c(ntr,k)
    c1[1,]<- inicio-1
    for (i in 2:ntr) {
    for (j in 1:k) {
    c1[i,j] <- c1[i-1,j]+1
    c1[i,j] <- c1[i,j]%%ntr
    }}
# randomize the elements from block, only if it is not row-col 
if( !rowcol) {
    for( i in 1:ntr) {
    nr<-1:k
    if(randomization)nr<-sample(1:k, replace=FALSE)
    c1[i,]<-c1[i,nr]
    }}
# randomize los bloques
#    nr<-sample(1:ntr, replace=FALSE)
    nr<-1:ntr
    if(randomization)nr<-sample(1:ntr, replace=FALSE)
    c1<-c1[nr,]
    c1<-t(c1)
    c1 <- as.numeric(c1) +1
    trt1<- trt[c1]
    book <- data.frame( irep, trt1)
    tr<-as.character(book[,2])
    dim(tr)<-c(k,ntr)
    design<-rbind(design,list(t(tr)))
    BOOK <-rbind(BOOK,book)
   }
    nr<-as.numeric(table(BOOK[,1]))
    nt<-length(nr)
    plots<-NULL
    for(i in 1:nt) plots<-c(plots,i*number+1:nr[i])
    block <- gl(nj*ntr, k)
    book <- data.frame(plots=plots,group=BOOK[,1],block=block,trt=BOOK[,2])
    names(book)[4] <- name.trt
    cat("\ncyclic design\n")
    cat("Generator block basic:\n")
    if (is.matrix(inicial)){
    for(i in 1:nj) cat(inicial[i,], "\n")
    } 
    else cat(inicial, "\n")
#   E <- (ntr - 1) * (r - 1)/((ntr - 1) * (r - 1) + r *(s - 1))
    cat("\nParameters\n===================")
    cat("\ntreatmeans :", ntr)
    cat("\nBlock size :", k)
    cat("\nReplication:", r, "\n")
    cat("\n")
#    cat("\nEfficiency factor\n(E )", E, "\n\n<<< Book >>>\n")
outdesign<-list(parameters=parameters,sketch=design,book=book)
return(outdesign)
}
