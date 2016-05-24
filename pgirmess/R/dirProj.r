dirProj<-function(df,deg=TRUE){
   df=as.matrix(df)
   if (!is.matrix(df) | ncol(df)!=4) stop("df must be a four column data.frame or matrix")
   x1<-df[,1]
   y1<-df[,2]
   alpha<-df[,3]
   d<-df[,4]
   if (deg) alpha<-alpha*pi/180
   ifelse (deg >=0 & deg <=pi/2, x2<-x1+d*sin(alpha),ifelse(deg >pi/2 & deg <=pi,x2<-x1+cos(alpha-pi/2),ifelse(deg >pi & deg <=1.5*pi,x2<-x1-sin(alpha-pi),x2<-x1-cos(alpha-1.5*pi))))
   ifelse (deg >=0 & deg <=pi/2, y2<-y1+d*cos(alpha),ifelse(deg >pi/2 & deg <=pi,y2<-y1-sin(alpha-pi/2),ifelse(deg >pi & deg <=1.5*pi,y2<-y1-cos(alpha-pi),y2<-y1+sin(alpha-1.5*pi))))
   cbind(x2,y2)
}