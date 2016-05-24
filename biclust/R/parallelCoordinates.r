parallelCoordinates=function(x, bicResult, number,plotBoth=FALSE, plotcol=TRUE, compare=TRUE,info=F,bothlab=c("Rows","Columns"), order = FALSE, order2 = 0, ylab="Value",col=1,...)
{
    bicRows=which(bicResult@RowxNumber[,number])
    bicCols=which(bicResult@NumberxCol[number,])
    if(order)
    {
        if(order2<1)
        {
            bicRows <- bicRows[order(rowSums(x[bicRows,bicCols]))]
            bicCols <- bicCols[order(colSums(x[bicRows,bicCols]))]
        }
        else {
            bicRows <- bicRows[order(bicResult@info$probsR[[number]])]
            bicCols <- bicCols[order(bicResult@info$probsC[[number]])]
        }

    }




  if(plotBoth)
    {
    op<-par(mfrow=c(2,1))
    if(compare)
      {
      mat<-x[bicRows,]
      matplot(mat,type='l',lty=1,col=gray(0.7),axes=F,xlab=bothlab[1],ylab=ylab,...)
      axis(1,at=1:nrow(mat),labels=rownames(mat),...)
      axis(2,...)
      matplot(mat[,bicCols],type='l',lty=1,lwd=2,add=T,col=col,...)

      mat<-t(x[,bicCols])
      matplot(mat,type='l',lty=1,col=gray(0.7),axes=F,xlab=bothlab[2],ylab=ylab,...)
      axis(1,at=1:nrow(mat),labels=rownames(mat),...)
      axis(2,...)
      matplot(mat[,bicRows],type='l',lty=1,lwd=2,add=T,col=col,...)
      }

    else
      {
      matplot(x[bicRows,bicCols],type='l',lty=1,lwd=2,axes=F,xlab=bothlab[1],ylab=ylab,col=col,...)
      axis(1,at=1:nrow(x[bicRows,bicCols]),labels=rownames(x[bicRows,bicCols]),...)
      axis(2,...)

      matplot(t(x[bicRows,bicCols]),type='l',lty=1,lwd=2,axes=F,xlab=bothlab[1],ylab=ylab,col=col,...)
      axis(1,at=1:nrow(t(x[bicRows,bicCols])),labels=rownames(t(x[bicRows,bicCols])),...)
      axis(2,...)
      }
    par(op)
    }
  else
    {
    if(plotcol)
      {
      if(compare)
        {
        mat<-x[bicRows,]
        matplot(mat,type='l',lty=1,col=gray(0.7),axes=F,xlab=bothlab[1],ylab=ylab,...)
        axis(1,at=1:nrow(mat),labels=rownames(mat),...)
        axis(2,...)
        matplot(mat[,bicCols],type='l',lty=1,lwd=2,add=T,col=col,...)

        }
      else
        {
        matplot(x[bicRows,bicCols],type='l',lty=1,lwd=2,axes=F,xlab=bothlab[1],ylab=ylab,col=col,...)
        axis(1,at=1:nrow(x[bicRows,bicCols]),labels=rownames(x[bicRows,bicCols]),...)
        axis(2,...)
        }
      }
    else
      {
      if(compare)
        {
        mat<-t(x[,bicCols])
        matplot(mat,type='l',lty=1,col=gray(0.7),axes=F,xlab=bothlab[2],ylab=ylab,...)
        axis(1,at=1:nrow(mat),labels=rownames(mat),...)
        axis(2,...)
        matplot(mat[,bicRows],type='l',lty=1,lwd=2,add=T,col=col,...)
        }
      else
        {
        matplot(t(x[bicRows,bicCols]),type='l',lty=1,lwd=2,axes=F,xlab=bothlab[1],ylab=ylab,col=col,...)
        axis(1,at=1:nrow(t(x[bicRows,bicCols])),labels=rownames(t(x[bicRows,bicCols])),...)
        axis(2,...)
        }
      }
    }
  if(info)
    {
    title(main=paste("Bicluster",number,"\n(rows=", length(bicRows),";", "columns=",length(bicCols),")",sep=" "))
    }


}




