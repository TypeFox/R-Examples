dplot <-
function(data, xlab="Variable x", ylab="Variable y", position=1, colors=TRUE, type="o", mean=TRUE){
       
means=function(data){
t=as.factor(data[,1])
d=data.frame(t,data[,-1])
s=split(data.frame(d[,-1]), d$t)
r=lapply(s, colMeans, na.rm=TRUE)
r=lapply(r, round,2)
rr=t(data.frame(r)); rr=data.frame(rr);rownames(rr)=NULL
treat=levels(t)
rr=data.frame(treat, rr); colnames(rr)=colnames(data)
return(rr)
}
	y=ylab
        x=xlab
        ti=type

rrr=means(data)
oi=ifelse(mean==TRUE,2,1)
l=list(data, rrr)
data=l[[oi]]


        d1=data[,1];d2=data[,-1]
        c1=1:length(d2)
        names=names(data)
	names=names[-1]
        c2=ifelse(colors==TRUE, 1,2)
        cor=list(c1,1)
        cor=cor[[c2]]
        t=list("top", "bottomright", "bottom", "bottomleft", "left", "topleft", "topright", "right", "center")
        position=position
        p=t[[position]]
        mat=matplot(d1, d2, type = ti, xlab = x, 
                ylab = y, col = cor, lty = c1, pch=c1)
        legend(p, names, bty = "n", col = cor, 
               lty =c1, pch=c1)

return(data)
    }
