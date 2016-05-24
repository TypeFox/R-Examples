histo1d<-function(dendat,binlkm,ala=NULL,yla=NULL,
pic=TRUE,brush=NULL,brushcol=c("blue"),col=NULL,border=NULL,
xlab="",ylab="",cex.lab=1,cex.axis=1,data=FALSE,
weights=rep(1,length(dendat)),normalization=TRUE,
height=NULL,subweights=NULL,graphplot=FALSE)
{
if (is.null(ala)) ala<-min(dendat)
if (is.null(yla)) yla<-max(dendat)
step<-(yla-ala)/binlkm
frekv<-matrix(0,binlkm,1)
value<-matrix(0,binlkm,1)
if (!is.null(brush)){
   cnum<-max(brush)
   shade <-matrix(0,binlkm,cnum)
}
if (!is.null(subweights)) taint<-matrix(0,binlkm,1)
n<-length(dendat)
for (i in 1:n){
   hava<-dendat[i]
   weight<-weights[i]
   ind<-min(binlkm,floor((hava-ala)/step)+1)
   frekv[ind]<-frekv[ind]+weight
   if ((!is.null(brush)) && (brush[i]>0)) 
              shade[ind,brush[i]]<-shade[ind,brush[i]]+1
   if (!is.null(subweights)) taint[ind]<-taint[ind]+n*subweights[i]
}
if (normalization) value<-frekv/(n*step) else value<-frekv
if (!is.null(brush)) shade<-shade/(n*step)
if ((normalization) && (!is.null(subweights))) taint<-taint/(n*step)

if (pic){
   if (is.null(height)) height<-max(value)
   plot(x="",y="",xlab=xlab,ylab=ylab,xlim=c(ala,yla),ylim=c(0,height),
   cex.lab=cex.lab,cex.axis=cex.axis)
   for (i in 1:binlkm){
          xala<-ala+(i-1)*step
          xyla<-xala+step
          y<-value[i]

          if (graphplot){
               if (i==1) yeka<-0 else yeka<-value[i-1]
               if (i==binlkm) ytok<-0 else ytok<-value[i]
               segments(xala,yeka,xala,ytok)
               segments(xala,ytok,xyla,ytok)
          } 
          else
          polygon(c(xala,xala,xyla,xyla),c(0,y,y,0),col=col,border=border)

          if (!is.null(brush)){
              y0<-0
              for (kk in 1:cnum){
                  y<-y0+shade[i,kk]
                  polygon(c(xala,xala,xyla,xyla),c(y0,y,y,y0),col=brushcol[kk])
                  y0<-y
              }
          }
          if (!is.null(subweights)){
              if (graphplot){
                 if (i==1) yeka<-0 else yeka<-taint[i-1]
                 if (i==binlkm) ytok<-0 else ytok<-taint[i]
                 segments(xala,yeka,xala,ytok,col=brushcol)
                 segments(xala,ytok,xyla,ytok,col=brushcol)
              } 
              else{
                 y<-taint[i]
                 polygon(c(xala,xala,xyla,xyla),c(0,y,y,0),col=brushcol)
              } 
          }
   }
}
if (data){
     return(list(frekv=frekv,ala=ala,step=step,value=value))
}
}



