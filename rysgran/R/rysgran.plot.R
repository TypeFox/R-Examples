rysgran.plot <-
function (x=NULL , y=NULL , data=NULL, output = "phi", lang="en-US", method="folk",
main = NULL, xlab = NULL, ylab = NULL, show.grid=TRUE, show.labels = FALSE, 
label.points = FALSE, pch = 19, col.labels = "black",
labels = NULL, col = "black", cex.labels = 1, cex.points = 1,
pos=1, z.cex.range=NULL, z=NULL, ...)

{
 tab1 <- as.data.frame (data)
 a<-0
 b<-0

tab<-gran.stats(data=tab1, output=output, method=method, verbal=FALSE, lang="en-US")

 if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")
 {
  title<-as.character(c("Bivariated Plot"))
  class<-"Verbal."
  if (x=="mean" | x=="Mean" | x=="M\u00E9dia" | x=="Media" | x=="m\u00E9dia" | x=="media")
  {
   a<-as.matrix(tab$Mean);colnames(a)<-c("Mean")
  }

  if (x=="Sorting" | x=="sorting" | x=="sort" | x=="Sort" | x=="Sele\u00E7\u00E3o" | x=="Sele\u00E7ao" | x=="Selecao" | x=="sele\u00E7\u00E3o" | x=="sele\u00E7ao" | x=="selecao" | x=="Sele" | x=="sele")
  {
   a<-as.matrix(tab$Sorting);colnames(a)<-c("Sorting")
  }

  if (x=="Skewness" | x=="skewness" | x=="Skew" | x=="skew" | x=="Assimetria" | x=="assimetria" | x=="Ass" | x=="ass")
  {
   a<-as.matrix(tab$Skewness);colnames(a)<-c("Skewness")
  }

  if (x=="Kurtosis" | x=="kurtosis" | x=="Kurt" | x=="kurt" | x=="Curtose" | x=="curtose" | x=="Curt" | x=="curt")
  {
   a<-as.matrix(tab$Kurtosis);colnames(a)<-c("Kurtosis")
  }

  if (y=="Mean" | y=="mean" | y=="M\u00E9dia" | y=="Media" | y=="m\u00E9dia" | y=="media")
  {
   b<-as.matrix(tab$Mean);colnames(b)<-c("Mean")
  }

  if (y=="Sorting" | y=="sorting" | y=="Sort" | y=="sort" | y=="Sele\u00E7\u00E3o" | y=="Sele\u00E7ao" | y=="Selecao" | y=="sele\u00E7\u00E3o" | y=="sele\u00E7ao" | y=="selecao" | y=="Sele" | y=="sele")
  {
   b<-as.matrix(tab$Sorting);colnames(b)<-c("Sorting")
  }

  if (y=="Skewness" | y=="skewness" | y=="Skew" | y=="skew" | y=="Assimetria" | y=="assimetria" | y=="Ass" | y=="ass")
  {
   b<-as.matrix(tab$Skewness);colnames(b)<-c("Skewness")
  }

  if (y=="Kurtosis" | y=="kurtosis" | y=="Kurt" | y=="kurt" | y=="Curtose" | y=="curtose" | y=="Curt" | y=="curt")
  {
   b<-as.matrix(tab$Kurtosis);colnames(b)<-c("Kurtosis") 
  }
 }

 if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")
 {
  title<-as.character(c("Gr\u00E1fico Bivariado"))
  class<-"Class."
  if (x=="mean" | x=="Mean" | x=="M\u00E9dia" | x=="Media" | x=="m\u00E9dia" | x=="media")
  {
   a<-as.matrix(tab$Mean);colnames(a)<-c("M\u00E9dia") 
  }

  if (x=="Sorting" | x=="sorting" | x=="sort" | x=="Sort" | x=="Sele\u00E7\u00E3o" | x=="Sele\u00E7ao" | x=="Selecao" | x=="sele\u00E7\u00E3o" | x=="sele\u00E7ao" | x=="selecao" | x=="Sele" | x=="sele")
  {
   a<-as.matrix(tab$Sorting);colnames(a)<-c("Sele\u00E7\u00E3o")
  }

  if (x=="Skewness" | x=="skewness" | x=="Skew" | x=="skew" | x=="Assimetria" | x=="assimetria" | x=="Ass" | x=="ass")
  {
   a<-as.matrix(tab$Skewness);colnames(a)<-c("Assimetria") 
  }

  if (x=="Kurtosis" | x=="kurtosis" | x=="Kurt" | x=="kurt" | x=="Curtose" | x=="curtose" | x=="Curt" | x=="curt")
  {
   a<-as.matrix(tab$Kurtosis);colnames(a)<-c("Curtose")
  }

  if (y=="Mean" | y=="mean" | y=="M\u00E9dia" | y=="Media" | y=="m\u00E9dia" | y=="media")
  {
   b<-as.matrix(tab$Mean);colnames(b)<-c("M\u00E9dia")
  }

  if (y=="Sorting" | y=="sorting" | y=="Sort" | y=="sort" | y=="Sele\u00E7\u00E3o" | y=="Sele\u00E7ao" | y=="Selecao" | y=="sele\u00E7\u00E3o" | y=="sele\u00E7ao" | y=="selecao" | y=="Sele" | y=="sele")
  {
   b<-as.matrix(tab$Sorting);colnames(b)<-c("Sele\u00E7\u00E3o")  
  }

  if (y=="Skewness" | y=="skewness" | y=="Skew" | y=="skew" | y=="Assimetria" | y=="assimetria" | y=="Ass" | y=="ass")
  {
   b<-as.matrix(tab$Skewness);colnames(b)<-c("Assimetria")
  }

  if (y=="Kurtosis" | y=="kurtosis" | y=="Kurt" | y=="kurt" | y=="Curtose" | y=="curtose" | y=="Curt" | y=="curt")
  {
   b<-as.matrix(tab$Kurtosis);colnames(b)<-c("Curtose") 
  }
 }

 c<-a-b 
 if (sum(c)==0) stop ("x and y cannot be equal")
 if (sum(a)==0) stop ("x must be a valid argument")
 if (sum(b)==0) stop ("y must be a valid argument")

 if (is.null(xlab)) xlab <- colnames(a)
 if (is.null(ylab)) ylab <- colnames(b)
 if (is.null(main)) main <- title
 if (label.points | show.labels)
 {
  if (is.null(labels)) labels<-row.names(tab) 
 }

 par(xpd=FALSE)

 if (show.labels)
 {
  label.points=FALSE
  z=NULL
  plot(b~a, xlab = xlab, ylab = ylab,main=main, type="n", pch = NA,...)
  text(a,b, labels=labels, col=col.labels, cex = cex.labels, pos=NULL)
 }

 if (label.points & is.null(z))
 {
  plot(b~a, xlab = xlab, ylab = ylab,main=main, type="p", pch = pch, cex = cex.points, col = col,...)
  text(a,b, labels=labels, pos=pos, col=col.labels, cex = cex.labels)
 }

 if(!is.null(z))
 {
  if(length(a) != length(z)) stop("z must be a vector with the same length of x")
  if(is.null(z.cex.range)) z.cex.range<-c(1,3)
  cex.bubbles <- TT.str(z, z.cex.range[1], z.cex.range[2])
  plot(b~a, xlab = xlab, ylab = ylab,main=main, type="n",...)
  points(a, b, pch = pch, col = col, type = "p", cex = cex.bubbles,...)
  text(a,b, labels=labels, pos=pos, col=col.labels, cex = cex.labels)
 }

 if(label.points==FALSE & show.labels==FALSE & is.null(z))
 {
  plot(b~a, xlab = xlab, ylab = ylab, main=main, type="p", 
  col = col, pch = pch, cex = cex.points,...)
 }
 
 
 if (show.grid==TRUE) 
 {
   if (output=="phi"){
     if (colnames(a)=="Mean" | colnames(a)=="M\u00E9dia"){
       classesx<-c(-8,-6,-2,-1,0,1,2,3,4,5,6,7,8,9)
       linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
       abline(v=linesx,lty=2)
     }
     if (colnames(b)=="Mean" | colnames(b)=="M\u00E9dia"){
       classesy<-c(-8,-6,-2,-1,0,1,2,3,4,5,6,7,8,9)
       linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
       abline(h=linesy,lty=2)
     }
     if (colnames(a)=="Sorting" | colnames(a)=="Sele\u00E7\u00E3o"){
       classesx<-c(0.35,0.5,0.7,1,2,4)
       linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
       abline(v=linesx,lty=2)
     }
     if (colnames(b)=="Sorting" | colnames(b)=="Sele\u00E7\u00E3o"){
       classesy<-c(0.35,0.5,0.7,1,2,4)
       linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
       abline(h=linesy,lty=2)
     }
     if (method=="moment"){
       if (colnames(a)=="Skewness" | colnames(a)=="Assimetria"){
         classesx<-c(-1.3,-0.43,0.43,1.3)
         linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
         abline(v=linesx,lty=2)
       }
       if (colnames(b)=="Skewness" | colnames(b)=="Assimetria"){
         classesy<-c(-1.3,-0.43,0.43,1.3)
         linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
         abline(h=linesy,lty=2)
       }
       if (colnames(a)=="Kurtosis" | colnames(a)=="Curtose"){
         classesx<-c(1.7,2.55,3.7,7.4,15)
         linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
         abline(v=linesx,lty=2)
       }
       if (colnames(b)=="Kurtosis" | colnames(b)=="Curtose"){
         classesy<-c(1.7,2.55,3.7,7.4,15)
         linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
         abline(h=linesy,lty=2)
       }
     }
     else{
       if (colnames(a)=="Skewness" | colnames(a)=="Assimetria"){
         classesx<-c(-1,-0.3,-0.1,0.1,0.3,1)
         linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
         abline(v=linesx,lty=2)
       }
       if (colnames(b)=="Skewness" | colnames(b)=="Assimetria"){
         classesy<-c(-1,-0.3,-0.1,0.1,0.3,1)
         linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
         abline(h=linesy,lty=2)
       }
       if (colnames(a)=="Kurtosis" | colnames(a)=="Curtose"){
         classesx<-c(0.67,0.9,1.11,1.5,3)
         linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
         abline(v=linesx,lty=2)
       }
       if (colnames(b)=="Kurtosis" | colnames(b)=="Curtose"){
         classesy<-c(0.67,0.9,1.11,1.5,3)
         linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
         abline(h=linesy,lty=2)
       }
     }
   }
   if (output=="metric"){
     if (colnames(a)=="Mean" | colnames(a)=="M\u00E9dia"){
       classesx<-c(2,4,8,16,31,63,125,250,500,1000,2000,4000,64000,256000)
       linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
       abline(v=linesx,lty=2)
     }
     if (colnames(b)=="Mean" | colnames(b)=="M\u00E9dia"){
       classesy<-c(2,4,8,16,31,63,125,250,500,1000,2000,4000,64000,256000)
       linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
       abline(h=linesy,lty=2)
     }
     if (colnames(a)=="Sorting" | colnames(a)=="Sele\u00E7\u00E3o"){
       classesx<-c(1.27,1.41,1.62,2,4,16)
       linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
       abline(v=linesx,lty=2)
     }
     if (colnames(b)=="Sorting" | colnames(b)=="Sele\u00E7\u00E3o"){
       classesy<-c(1.27,1.41,1.62,2,4,16)
       linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
       abline(h=linesy,lty=2)
     }
     if (method=="moment"){
       if (colnames(a)=="Skewness" | colnames(a)=="Assimetria"){
         classesx<-c(-1.3,-0.43,0.43,1.3)
         linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
         abline(v=linesx,lty=2)
       }
       if (colnames(b)=="Skewness" | colnames(b)=="Assimetria"){
         classesy<-c(-1.3,-0.43,0.43,1.3)
         linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
         abline(h=linesy,lty=2)
       }
       if (colnames(a)=="Kurtosis" | colnames(a)=="Curtose"){
         classesx<-c(1.7,2.55,3.7,7.4,15)
         linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
         abline(v=linesx,lty=2)
       }
       if (colnames(b)=="Kurtosis" | colnames(b)=="Curtose"){
         classesy<-c(1.7,2.55,3.7,7.4,15)
         linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
         abline(h=linesy,lty=2)
       }
     }
     else{
       if (colnames(a)=="Skewness" | colnames(a)=="Assimetria"){
         classesx<-c(-1,-0.3,-0.1,0.1,0.3,1)
         linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
         abline(v=linesx,lty=2)
       }
       if (colnames(b)=="Skewness" | colnames(b)=="Assimetria"){
         classesy<-c(-1,-0.3,-0.1,0.1,0.3,1)
         linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
         abline(h=linesy,lty=2)
       }
       if (colnames(a)=="Kurtosis" | colnames(a)=="Curtose"){
         classesx<-c(0.67,0.9,1.11,1.5,3)
         linesx<-classesx[which(classesx[] >= min(a) & classesx[] <= max(a))]
         abline(v=linesx,lty=2)
       }
       if (colnames(b)=="Kurtosis" | colnames(b)=="Curtose"){
         classesy<-c(0.67,0.9,1.11,1.5,3)
         linesy<-classesy[which(classesy[] >= min(b) & classesy[] <= max(b))]
         abline(h=linesy,lty=2)
       }
     }
   }
   
   if(length(linesx)!=0){
     medx<-numeric(length=(length(linesx)+1))
     for(i in 1:(length(linesx)+1)){
       if (i==1){
         medx[i]<-mean(c(linesx[1],min(a)))     
       }
       if (i==(length(linesx)+1)){
         medx[i]<-mean(c(linesx[length(linesx)],max(a)))
       }
       if (i>1 & i<(length(linesx)+1)){
         medx[i]<-mean(c(linesx[i-1],linesx[i]))
       }
     }
     mtext(LETTERS[1:(length(linesx)+1)], side=3, line=0.3,
           at=medx, adj=c(1,(rep(0.5,length(medx)-2)),0))
   }
   
   if (length(linesy)!=0){
     medy<-numeric(length=(length(linesy)+1))
     for(i in 1:(length(linesy)+1)){
       if (i==1){
         medy[i]<-mean(c(linesy[1],min(b)))     
       }
       if (i==(length(linesy)+1)){
         medy[i]<-mean(c(linesy[length(linesy)],max(b)))
       }
       if (i>1 & i<(length(linesy)+1)){
         medy[i]<-mean(c(linesy[i-1],linesy[i]))
       }
     }
     mtext(c(1:(length(linesy)+1)), side=4, line=0.3,
           at=medy, las=2, padj=c(1,(rep(0.5,length(medy)-2)),0))
     }
   
   if(length(linesx)==0){
     menor<-classesx[which(classesx[] <= min(a))]
     linesx<-c(menor[length(menor)],classesx[which(classesx[] >= max(a))][1])
     abline(v=linesx,lty=2)
     medx<-mean(a)
     mtext(LETTERS[1:(length(linesx)-1)], side=3, 
           line=0.3, at=medx)
   }
   
   if(length(linesy)==0){
     menor<-classesy[which(classesy[] <= min(b))]
     linesy<-c(menor[length(menor)],classesy[which(classesy[] >= max(b))][1])
     abline(h=linesy,lty=2)
     medy<-mean(b)
     mtext(c(1:(length(linesy)-1)), side=4, line=0.3,
           las=2, at=medy)
   }
 
 tab.leg<-matrix(nrow=max(length(medx),length(medy)),ncol=4)
 tab.leg<-as.data.frame(tab.leg)
 colnames(tab.leg)<-c(colnames(a),(paste(class,colnames(a),sep="")),colnames(b),(paste("Verbal.",colnames(b),sep="")))
 
 length(medx)<-nrow(tab.leg)
 length(medy)<-nrow(tab.leg)
 
 tab.leg[,1]<-medx
 tab.leg[,3]<-medy
 
 for(i in 1:nrow(tab.leg)){
   if(is.na(tab.leg[i,1]))next
   if(colnames(a)=="Mean" | colnames(a)=="M\u00E9dia"){
     tab.leg[i,2]<-class.mean(tab.leg[i,1], output=output, lang=lang)
     }
   if(colnames(a)=="Sorting" | colnames(a)=="Sele\u00E7\u00E3o"){
     tab.leg[i,2]<-class.sort(tab.leg[i,1], method=method, output=output, lang=lang)
     }
   if(colnames(a)=="Skewness" | colnames(a)=="Assimetria"){
     tab.leg[i,2]<-class.skew(tab.leg[i,1], method=method, output=output, lang=lang)
     }
   if(colnames(a)=="Kurtosis" | colnames(a)=="Curtose"){
     tab.leg[i,2]<-class.kurt(tab.leg[i,1], method=method, output=output, lang=lang)
     }
   }
 for(i in 1:nrow(tab.leg)){
   if(is.na(tab.leg[i,3]))next
   if(colnames(b)=="Mean" | colnames(b)=="M\u00E9dia"){
     tab.leg[i,4]<-class.mean(tab.leg[i,3], output=output, lang=lang)
   }
   if(colnames(b)=="Sorting" | colnames(b)=="Sele\u00E7\u00E3o"){
     tab.leg[i,4]<-class.sort(tab.leg[i,3], method=method, output=output, lang=lang)
   }
   if(colnames(b)=="Skewness" | colnames(b)=="Assimetria"){
     tab.leg[i,4]<-class.skew(tab.leg[i,3], method=method, output=output, lang=lang)
   }
   if(colnames(b)=="Kurtosis" | colnames(b)=="Curtose"){
     tab.leg[i,4]<-class.kurt(tab.leg[i,3], method=method, output=output, lang=lang)
   }
   }

 let<-LETTERS[1:(length(linesx)+1)]
 num<-(1:(length(linesy)+1))
 for (i in 1:nrow(tab.leg)){
   if(is.na(tab.leg[i,1]))next
   tab.leg[i,1]<-let[i]
 }
 for (i in 1:nrow(tab.leg)){
   if(is.na(tab.leg[i,3]))next
   tab.leg[i,3]<-num[i]
 }
 for (i in 1:nrow(tab.leg)){
   for (j in 1:ncol(tab.leg)){
     if (is.na(tab.leg[i,j])){
       tab.leg[i,j]<-""
     }
   }
 }
 return(tab.leg)
 }
}#fim da fun\u00E7ao

