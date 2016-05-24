interactive2D=function(data,cl){
  require("tkrplot")
  if (is.matrix(data) == FALSE || ncol(data) != 2)
     stop("data must be a matrix with two columns")
  if (is.vector(cl) == FALSE || nrow(data) != length(cl))
     stop("cl must be a vector with length equaling the number of rows of data")
  if (any(is.na(data)) || any(is.na(cl)))
     stop("no missing values are allowed")

  clnew <- unclass(as.factor(cl))
  nc <- max(clnew)
  if (nc != 2) 
     stop("only two classes allowed in cl")

  data_1 = subset(data,clnew==1)
  data_2 = subset(data,clnew==2) 
  lcd_1<-mlelcd(data_1)
  lcd_2<-mlelcd(data_2)
  A1 = hatA(lcd_1)
  A2 = hatA(lcd_2)

  px<-seq(min(data[,1]),max(data[,1]),(max(data[,1])-min(data[,1]))/59)
  py<-seq(min(data[,2]),max(data[,2]),(max(data[,2])-min(data[,2]))/59)
  pxx<-sort(rep(px,60))
  pyy<-rep(py,60)
  x<-matrix(c(pxx,pyy),ncol=2)


  lh1<-dslcd(x,lcd_1,A1) * dim(data_1)[1]
  lh2<-dslcd(x,lcd_2,A2) * dim(data_2)[1]
 
  tt<-tktoplevel()
  start<-0
  logratio=tclVar(start)

  plotdens=function(...){
    ratio=exp(as.numeric(tclvalue(logratio)))
    contour(px,py,t(matrix(lh1>lh2*ratio,60,60)),0.5,xlab="x", ylab="y",sub=paste("Risk ratio = ",format(ratio,digits=6)))
    points(x,col=ifelse(lh1>lh2*ratio, "red", "green"),pch=46)
    points(data,col=ifelse(clnew==1,"red", "green"))
  }

  img = tkrplot(tt,plotdens)
  densplot = function(...){
    tkrreplot(img)
  }

  scl=tkscale(tt, command=densplot, from=log(min(lh1/lh2))-1, to=log(max(lh1/lh2))+1, showvalue=FALSE, variable=logratio, resolution=0.01, orient='hori')

  label <- tklabel(tt, text="Slider for risk ratio on log-scale")
  tkpack(img,side='top')
  tkpack(label)
  tkpack(scl)
}
