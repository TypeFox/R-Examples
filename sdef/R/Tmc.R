Tmc <-
function(iter=1000,output.ratio){
load(paste(output.ratio$dataname,".Rdata"))
if(output.ratio$pvalue==FALSE){
data=1-data
}

dim1=dim(data)[1]
lists = ncol(data)
l=length(output.ratio$Common)
Tmax = max(output.ratio$ratios,na.rm=TRUE)
Tmax.null = rep(NA,iter)
ratios.null = matrix(NA,l,iter)
sample = matrix(NA,dim1,lists)
    
sample[,1] <- data[,1]
for(k in 1:iter){
int = c()
L=matrix(0,l,lists)
data1 = matrix(NA,dim1,lists)
data1[,1] <- data[,1]

for(j in 2:lists){
sample[,j] = sample(data[,j])
data1[,j] = sample[,j]
}
    
threshold = output.ratio$h
for(i in 1:l){
temp = data1<=threshold[i]
for(j in 1:lists){
L[i,j] <- sum(temp[,j])
temp[temp[,j]==FALSE,j]<-0
temp[temp[,j]==TRUE,j]<-1
        }
        int[i] <- sum(apply(temp,1,sum)==lists)
}


expected = apply(L,1,prod)/(dim1)^(lists-1)
observed = int
ratios = matrix(0,l,1)

for(i in 1:l){
ratios[i,1] <- observed[i]/expected[i]
}
ratios.null[,k] <- ratios
ratios <- ratios[threshold>0]
Tmax.null[k] = max(ratios)
if(k%%10==0)cat(k,"of",iter,"completed","\n")
}

ID=seq(1,iter)
p=length(ID[Tmax.null>=Tmax])
pvalue<- p/iter

#########
XX=hist(Tmax.null,plot=FALSE)
    hist(Tmax.null,main=expression(paste("Distribution of ", T(h[max])," under independence")), xlab = "T", ylab = "", xaxt = "n", 
        cex.main = 0.7, xlim = c(min(Tmax.null), max(c(Tmax, 
            max(Tmax.null)))), yaxt = "n", cex.axis = 0.9)
    axis(1, at = seq(min(Tmax.null), max(c(Tmax, max(Tmax.null))), 
        length.out = 10), labels = round(seq(min(Tmax.null), 
        max(c(Tmax, max(Tmax.null))), length.out = 10), 2))
if(pvalue>0){
legend(x=XX$breaks[length(XX$breaks)/2],y=max(XX$counts)/1.5,legend=paste("P value =",pvalue),bty="n",cex=0.9)
abline(v=Tmax,lty=2)
}
if(pvalue==0){
        legend(x = XX$breaks[length(XX$breaks)/2], 
            y = max(XX$counts)/1.5, legend = paste("P value <", 1/iter), 
            bty = "n", cex = 0.9)
        abline(v = Tmax, lty = 2)
}
#########

postscript("MC p-value.ps")
XX=hist(Tmax.null,plot=FALSE)
    hist(Tmax.null,main=expression(paste("Distribution of ", T(h[max])," under independence")),xlab = "T", ylab = "", xaxt = "n", 
        cex.main = 0.7, xlim = c(min(Tmax.null), max(c(Tmax, 
            max(Tmax.null)))), yaxt = "n", cex.axis = 0.9)
    axis(1, at = seq(min(Tmax.null), max(c(Tmax, max(Tmax.null))), 
        length.out = 10), labels = round(seq(min(Tmax.null), 
        max(c(Tmax, max(Tmax.null))), length.out = 10), 2))
if(pvalue>0){
legend(x=XX$breaks[length(XX$breaks)/2],y=max(XX$counts)/1.5,legend=paste("P value =",pvalue),bty="n",cex=0.9)
abline(v=Tmax,lty=2)
dev.off()
return(list(pvalue=pvalue))
}
if(pvalue==0){
        legend(x = XX$breaks[length(XX$breaks)/2], 
            y = max(XX$counts)/1.5, legend = paste("P value <", 1/iter), 
            bty = "n", cex = 0.9)
        abline(v = Tmax, lty = 2)
dev.off()
return(noquote(paste("pvalue < ",1/iter)))
}

}

