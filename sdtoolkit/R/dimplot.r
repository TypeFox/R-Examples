
#Function for plotting the relative dimension restrictions

dimplot <- function(boxseq, boxind=1, alldims=FALSE, thetitle=NULL, incol="lightblue", longcol="black"){#,outcol="transparent"){

  outcol <- "transparent"
  
  relbounds <- boxseq[[boxind]]$relbox

  if(is.null(thetitle)){thetitle <- "Normalized dimension restrictions"}

  relbounds[1,relbounds[1,]<0] <- 0
  relbounds[2,relbounds[2,]>1] <- 1

  if(!alldims){
    span <- relbounds[2,]-relbounds[1,]
    relbounds <- relbounds[,(span<1)]
  }

  relo <- relbounds[1,]
  rehi <- relbounds[2,]
    
  firstnum  <- relo
  midnum <- rehi-relo
  hinum  <- 1-rehi
   
  bmat <- rbind(firstnum,midnum,hinum)
  colnames(bmat) <- colnames(relbounds)
  
  ndims <- ncol(bmat)
  
  fakemat <- matrix(rep(1,3*ndims),ncol=ndims)
  
  #try setting margins to account for potentially long varnames:  
  maxchar <- max(nchar(colnames(relbounds)))
  leftmar <- max(maxchar*.8, 4) #approximate charter widths as .8 height, use default if tiny
  
  par(mar=c(5,leftmar,4,2)+.1)
                                
#  if(topdown){ 
  #  par(mai=c(1,1,.5,.5))

    barplot(fakemat,horiz=TRUE, add=FALSE, col="transparent", beside = TRUE, width=1, space=c(0,3),border=c("transparent",longcol,"transparent"))

    barplot(bmat[,ndims:1],add=TRUE,col=c(outcol,incol,outcol),horiz=TRUE,border=FALSE,names.arg=rep("",ndims),axes=FALSE, main=thetitle, width=3,space=1)
#    barplot(fakemat,horiz=TRUE, add=TRUE, col="transparent", beside = TRUE, width=.33333, space=c(0,.6),border=c("transparent",incol,"transparent"))
    axis(1)
    axis(2,at = seq(4.5,by=6,length.out=ndims),labels=colnames(bmat[,ndims:1]),las=1)


}

#  
#  else{
#    barplot(bmat,col=c(outcol,incol,outcol),horiz=TRUE,border=FALSE)
#    barplot(fakemat,horiz=TRUE, add=TRUE, col="transparent", beside = TRUE, width=.33333, space=c(0,.6),border=c("transparent",incol,"transparent"))
#  }
#    
#}
