`plotDilutionSeries` <-
function (D0=2, data.dilutes,sensible.min=5, sensible.max=1.e9,minimal.err=5, plot=T, title="") {
K = ncol(data.dilutes) # number of dilution steps in a dilution series
nsample = nrow(data.dilutes) # number of samples
D0 # this is a preset value in diluteion experiments, typical D0=10, 3, or 2.
a= max(5, min(data.dilutes, na.rm=T))
M= max(data.dilutes, na.rm =T) 
D=D0
S.y = as.numeric(data.dilutes[,-K]) # get all columns except the last 
S.x = as.numeric(data.dilutes[,-1])            # get all columns except the first
if(plot) {
    plot(S.x,S.y,pch='+',col='black',xlab='Signal at next dilution step', ylab='Signal', main=title)
    abline(0,1) # identity line
}
filter = ((S.x < sensible.min) | (S.y <sensible.min))
if(plot) {
    points(S.x[filter],S.y[filter], pch='+',col='red')
}
filter = ((S.x > sensible.max) | (S.y > sensible.max))
if(plot) {
    points(S.x[filter],S.y[filter], pch='o',col='red')
}
S.x1=S.x
S.x1[S.x1<max(0,sensible.min)] =max(0,sensible.min)
S.y1=S.y
S.y1[S.y1<max(0,sensible.min)] =max(0,sensible.min)
ratio= log(abs(S.y1)/abs(S.x1)) #log ratio
ratio.median= median(ratio, na.rm=T)
ratio.mad =mad(ratio,na.rm=T)
filter  =  abs(ratio - ratio.median)/3 > ratio.mad
if(plot) {
    points(S.x[filter],S.y[filter], pch='+',col='red')
}
S.x[filter] = NA
S.x[filter] = NA
S.x[S.x <sensible.min] =NA
S.y[S.y <sensible.min] =NA
S.x[S.x >sensible.max] =NA
S.y[S.y >sensible.max] =NA
fit = nls(S.y ~ a +1/((1/(S.x -a) -c)/D+c),start=list(a=a,D=D,c=1/M),algorithm="port", lower=list(minimal.err,1,0),weights=1/(minimal.err+abs(S.x)))
summary(fit)
a=summary(fit)$parameter[1]
D=summary(fit)$parameter[2]
c=summary(fit)$parameter[3]
d.a= summary(fit)$parameter[4]
d.D= summary(fit)$parameter[5]
d.c= summary(fit)$parameter[6]
M=a+1/summary(fit)$parameter[3]
S.pred = a + 1/((1/(sort(S.x) -a) -1/(M-a))/D+1/(M-a))
if(plot)  {
    points(sort(S.x), S.pred, pch='+', col='blue', type='l',lwd=2)
    #mtext(paste('a=', format(a,digits=3),' M =', format(M,digits=3),' Dilution factor',format(D,digits=3),'\n\n'))
    title(main=paste('a=', format(a,digits=3),' M =', format(M,digits=3),' Dilution factor',format(D,digits=3),'\n\n')
         ,sub=title) 
}
list(D=D,c=c, a=a, d.D=d.D, d.c=d.c,d.a=d.a, M=M)
} #end of function 

