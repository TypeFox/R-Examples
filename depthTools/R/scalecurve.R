scalecurve<-function(x,y=NULL,xlab="p",ylab="A(p)",main="Scale curve",lwd=2,...)
{
#Inputs: 
         # x = Data matrix, where the rows are the observations and the columnns are the variables
         # y = optional vector, either numeric or factor, indicating the class of each observation in x
         # ... = graphical parameters (except xlab and ylab)
#Outputs: 
         # r = scale curve defined on previous grid; if y is provided, then r is a list with the scale curve of each class in each component

x<-as.matrix(x)
if (ncol(x)==1) x<-t(x)
if (length(y)==0)
{
n<-nrow(x);p<-ncol(x);    # size of the data matrix
t<-1:p;
arg<-matrix(0,1,n);       # dimension of the argument vector 
arg[1]<-1/n;              # first value of the vector arg
cont<-MBD(x,plotting=FALSE)$MBD          # calculate the modified band depth of each observations within the sample
I<-order(cont);           # order these observations from the one with lowest depth to the one with highest depth
xx<-x[I[n],];             # Deepest curve from the sample
r<-matrix(0,1,n);         # inizialize the the vector r as zero

for (j in 2:n)
{
   a<-x[I[n-j+1],]
   xx<-rbind(xx,a)        # curves that delimit the band 
   M <- apply(xx,2,max); m <- apply(xx,2,min)   # area of the band
   aa<-(t[p]-t[1])/(p-1);
   sM <- ( aa*sum(M[2:p]) + aa*sum(M[1:(p-1)]))/2
   sm <- (aa*sum(m[2:p]) + aa*sum(m[1:(p-1)]))/2
   r[j]<-sM-sm            # area of the band
   arg[j]<-j/n;           # proportion of curves that define the band 
}

arg <- c(0,arg); r <- c(0,r)
plot(arg,r,ty="l",xlab=xlab,ylab=ylab,main=main,lwd=lwd,...)
}
else
{
y<-as.matrix(y)
n<-nrow(x);p<-ncol(x)
if (nrow(y)!=n) {stop("Length of y mismatches the dimension of x")}
y <- as.factor(y)
r <- list(); arg<- list()
t<-1:p;
mx <- 0
for (class in 1:length(levels(y)))
  {
    n <- sum(y==levels(y)[class])
    xc <- x[y==levels(y)[class],]
    arg[[class]]<-matrix(0,1,n)    # dimension of the argument vector 
    arg[[class]][1]<-1/n           # first value of the vector arg
    cont<-MBD(xc, plotting=FALSE)$MBD              # calculate the modified band depth of each observations within the sample
    I<-order(cont)                 # order these observations from the one with lowest depth to the one with highest depth
    xx<-xc[I[n],]                  # Deepest curve from the sample
    r[[class]]<-matrix(0,1,n)      # inizialize the the vector r as zero

    for (j in 2:n)
     {
       a<-xc[I[n-j+1],]
       xx<-rbind(xx,a)             # curves that delimit the band 
       M <- apply(xx,2,max); m <- apply(xx,2,min); # area of the band
       aa<-(t[p]-t[1])/(p-1); 
       sM <- ( aa*sum(M[2:p]) + aa*sum(M[1:(p-1)]))/2
       sm <- (aa*sum(m[2:p]) + aa*sum(m[1:(p-1)]))/2
       r[[class]][j]<-sM-sm;  
       arg[[class]][j]<-j/n;       # proportion of curves that define the band 
     }

    arg[[class]] <- c(0,arg[[class]]); r[[class]] <- c(0,r[[class]])
    mx <- max(c(mx,r[[class]]))
  }  ##  end FOR class

plot(seq(0,1,0.1),c(0,rep(mx,10)),ty="n",xlab=xlab,ylab=ylab,main=main,lwd=lwd,...)
for (class in 1:length(levels(y)))
  {
    lines(arg[[class]],r[[class]], lwd=2, lty=class)
  }  ### end FOR class
legend("topleft",legend=levels(y),lty=1:length(levels(y)),lwd=2)
}  ## end ELSE
return(r)
}




