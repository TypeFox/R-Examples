#this draw the coefficient plots for single stage coefficients
plot.linearcl<-function(x,index=NULL,names=NULL,ylab='std coefficients',xlab='',col='gray',...){
  p=length(x$beta)
  if (is.null(index)) index=1:p
  if (is.null(names)) names=paste('V',index,sep='')
names=factor(names,levels=names)
co=x$beta[index]/sqrt(sum(x$beta[index]^2))
method=rep('Olearn',p)
data1=data.frame(co,method,names)
f<-ggplot(data1, aes(names, co))+
  ylab(ylab) +
  xlab(xlab) +
  #scale_x_discrete(labels=names1)
  geom_bar(stat="identity",width=0.5,fill=col) + 
  facet_wrap(~method)+
  coord_flip() +
  #scale_fill_grey() +
  theme_bw(base_size=11) 
  suppressWarnings(print(f))
}

plot.qlearn<-function(x,index=NULL,names=NULL,ylab='std coefficients',xlab='',col='gray',...){
  p=length(x$co)
  if (is.null(index)) index=(p/2+1):p
  if (is.null(names)) names=paste('V',index-p/2,sep='')
  names=as.factor(names)
  co=x$co[index]/sqrt(sum(x$co[index]^2))
  method=rep('Qlearning',p)
  data1=data.frame(co,method,names)
  f<-ggplot(data1, aes(names, co))+
    ylab(ylab) +
    xlab(xlab) +
    #scale_x_discrete(labels=names1)
    geom_bar(stat="identity",width=0.5,fill=col) + 
    facet_wrap(~method)+
    coord_flip() +
    #scale_fill_grey() +
    theme_bw(base_size=11) 
    suppressWarnings(print(f))
}


