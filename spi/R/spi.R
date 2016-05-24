spi <-
function(nargs,filename,id,fd,title,output,txlab,tylab) 
{

 options(max.print=9999999)

 dados1=read.table(filename, header = TRUE)

 ncolumn=nrow(dados1)

 for(m in 2:ncol(dados1))
 {
    if(m==2)
       dat_aux=dados1[,m]
    else
	 {
         dat_aux2=rbind(dat_aux,dados1[,m])
         dat_aux=dat_aux2
	 } 
 }

 dados=dat_aux

 nlin=nrow(dados)
 if (nargs<3){
    return("Error: very small number of arguments")
 }  else {

          if (nargs>7){
             return("Error: very large number of arguments") }

          else if (nargs==3){
             tit=paste("Standard Precipitation Index", sep="") 
             output=2
          }
          else if (nargs==4){
             output=2
          }
          else {
                tit=paste(title,sep="") }
          ndates= fd - id + 1

         }

 data=dados[,1:ncolumn]
 dates=seq(id,fd,1)
 quad=1:ncolumn
 dates_ts=seq(id,fd,1)

 for(i in 2:ncolumn)
 {
    if (i==2) {
       data_aux1=rbind( t(data[i-1,]),t(data[i,])) }
    else {
          data_aux=rbind( data_aux1, t(data[i,]) )
         }
 }
 data2=data_aux

 #Routine for Calculating SPI
 fit.cdf <- ecdf(data2)
 cdfs <- sapply(data,fit.cdf)
 SPI <- qnorm(cdfs)

 #SPI graph
 spi.breaks <- c(-2.4,-2,-1.6,-1.3,-0.8,-0.5,0.5,0.8,1.3,1.6,2,2.4)
 spi.cols <- colorRampPalette(c("darkred","red","yellow","white","green","blue","darkblue"),space="rgb")
 spi.plot <- matrix(unlist(SPI),nlin,ncolumn) #convert list to matrix for plotting
 spi.plot[(spi.plot==Inf)] <- 2.2 #necessary to remove Infs because ecdf is being used
 spi.plot[(spi.plot==-Inf)] <- -2.2 #necessary to remove Infs because ecdf is being used

 ##plot diagram
 if (nargs>2){
    if (output==1)
    {
       if (nargs==7)
       {
          filled.contour(dates_ts,quad,matrix(unlist(spi.plot[,1:ncolumn]),ndates,ncolumn),col=spi.cols(11),xlab=txlab,ylab=tylab,cex.lab=1.1,
               	       font.lab=2,levels=spi.breaks,key.title=tit, las=2, cex.axis=.7,
           	             plot.axes = { axis(1, dates) 
                                       axis(2, quad) } )
          title(main=title,cex.main=1.1)
       } else {
               filled.contour(dates_ts,quad,matrix(unlist(spi.plot[,1:ncolumn]),ndates,ncolumn),col=spi.cols(11),xlab="dates",ylab="variables",cex.lab=1.1,
               	            font.lab=2,levels=spi.breaks,key.title=tit, las=2, cex.axis=.7,
           	                  plot.axes = { axis(1, dates) 
                                            axis(2, quad) } )
               title(main=title,cex.main=1.1)
              }


    } else {
            m=data.frame(dates=dates,spi.plot[,1:ncolumn] )
            #print(spi.plot)
            print(m)
           }

  } else {
           m=data.frame(dates=dates,spi.plot[,1:ncolumn] )
           #print(spi.plot)
           print(m)
         }


} #End spi function

