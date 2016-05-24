NormalizationEffect <-
function(MaxOrder,SingleEffect,TwoEffect,ThreeEffect,
                                FourEffect,FiveEffect,SaveFileName=""){
  # Normalization of SingleEffect, TwoEffect, ThreeEffect, FourEffect and FiveEffect
  #
  # input
  #    MaxOrder: the specified maximum order, must be setted as 1,2,3,4 or 5.
  #    SingleEffect: There are 2 columns. The first column saves all SNPs, and the 
  #                  second column saves their corresponding effects. Descending 
  #                  save according to their effects.                   
  #    TwoEffect: There are 3 columns. The first 2 columns save all 2-SNP combinations,
  #               and the last column saves their corresponding effects. Descending 
  #               save according to their interaction effects.
  #    ThreeEffect: There are 4 columns. The first 3 columns save all 3-SNP combinations,
  #                 and the last column saves their corresponding effects. Descending 
  #                 save according to their interaction effects.
  #    FourEffect: There are 5 columns. The first 4 columns save all 4-SNP combinations,
  #                and the last column saves their corresponding effects. Descending 
  #                save according to their interaction effects.
  #    FiveEffect: There are 6 columns. The first 5 columns save all 5-SNP combinations,
  #                and the last column saves their corresponding effects. Descending
  #                save according to their interaction effects.
  #    SaveFileName: basic file name for saving figure.
  #
  # output
  #    SingleEffect: The normalization of SingleEffect
  #    TwoEffect: The normalization of TwoEffect
  #    ThreeEffect: The normalization of ThreeEffect
  #    FourEffect: The normalization of FourEffect
  #    FiveEffect: The normalization of FiveEffect
  #
  # Junliang Shang
  # 4.6/2014
  
  if(MaxOrder==2){
    meanSingle <- mean(SingleEffect[,2])
    meanTwo <- mean(TwoEffect[,3])
    
    rowSingle <- nrow(SingleEffect)
    rowTwo <- nrow(TwoEffect)
    
    x <- seq(1, rowSingle+rowTwo)
    y <- c(SingleEffect[,2],TwoEffect[,3])
    
    SingleEffect[,2] <- SingleEffect[,2]*meanTwo/meanSingle
    
    z <- c(SingleEffect[,2],TwoEffect[,3])
    
    NmeanSingle <- mean(SingleEffect[,2])
    NmeanTwo <- meanTwo
    
    # plot
    x.pch <- rep(5,length(x))
    x.pch[1:rowSingle] <- 1
    
    x.col <- rep("red",length(x))
    x.col[1:rowSingle] <- "green"
    
    dev.new()
    set.seed(8)
    layout(matrix(1:2, 2, 1))
    
    plot(log(x),y,pch=x.pch,col=x.col,xlab=NA,ylab="Original Effect")
    lines(log(1:rowSingle),rep(meanSingle,rowSingle),col="green")
    lines(log((rowSingle+1):(rowSingle+rowTwo)),rep(meanTwo,rowTwo),col="red")
    
    plot(log(x),z,pch=x.pch,col=x.col,xlab=NA,ylab="Normalized Effect")
    lines(log(1:rowSingle),rep(NmeanSingle,rowSingle),col="green")
    lines(log((rowSingle+1):(rowSingle+rowTwo)),rep(NmeanTwo,rowTwo),col="red")
    # savePlot(filename=paste(SaveFileName,"_N",sep=""),type="tiff")
  }
  
  if(MaxOrder==3){
    
    meanSingle <- mean(SingleEffect[,2])
    meanTwo <- mean(TwoEffect[,3])
    meanThree <- mean(ThreeEffect[,4])
    
    rowSingle <- nrow(SingleEffect)
    rowTwo <- nrow(TwoEffect)
    rowThree <- nrow(ThreeEffect)
    
    x <- seq(1, rowSingle+rowTwo+rowThree)
    y <- c(SingleEffect[,2],TwoEffect[,3],ThreeEffect[,4])
    
    SingleEffect[,2] <- SingleEffect[,2]*meanTwo/meanSingle
    ThreeEffect[,4] <- ThreeEffect[,4]*meanTwo/meanThree
    
    z <- c(SingleEffect[,2],TwoEffect[,3],ThreeEffect[,4])
    
    NmeanSingle <- mean(SingleEffect[,2])
    NmeanTwo <- meanTwo
    NmeanThree <- mean(ThreeEffect[,4])
    
    # plot
    x.pch <- rep(5,length(x))
    x.pch[1:rowSingle] <- 1
    x.pch[(rowSingle+1):(rowSingle+rowTwo)] <- 2
    
    x.col <- rep("red",length(x))
    x.col[1:rowSingle] <- "green"
    x.col[(rowSingle+1):(rowSingle+rowTwo)] <- "blue"
    
    dev.new()
    set.seed(8)
    layout(matrix(1:2, 2, 1))
    
    plot(log(x),y,pch=x.pch,col=x.col,xlab=NA,ylab="Original Effect")
    lines(log(1:rowSingle),rep(meanSingle,rowSingle),col="green")
    lines(log((rowSingle+1):(rowSingle+rowTwo)),rep(meanTwo,rowTwo),col="blue")
    lines(log((rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)),
          rep(meanThree,rowThree),col="red")
    
    plot(log(x),z,pch=x.pch,col=x.col,xlab=NA,ylab="Normalized Effect")
    lines(log(1:rowSingle),rep(NmeanSingle,rowSingle),col="green")
    lines(log((rowSingle+1):(rowSingle+rowTwo)),rep(NmeanTwo,rowTwo),col="blue")
    lines(log((rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)),
          rep(NmeanThree,rowThree),col="red") 
    # savePlot(filename=paste(SaveFileName,"_N",sep=""),type="tiff")
  }
  
  if(MaxOrder==4){
    
    meanSingle <- mean(SingleEffect[,2])
    meanTwo <- mean(TwoEffect[,3])
    meanThree <- mean(ThreeEffect[,4])
    meanFour <- mean(FourEffect[,5])
    
    rowSingle <- nrow(SingleEffect)
    rowTwo <- nrow(TwoEffect)
    rowThree <- nrow(ThreeEffect)
    rowFour <- nrow(FourEffect)
    
    x <- seq(1, rowSingle+rowTwo+rowThree+rowFour)
    y <- c(SingleEffect[,2],TwoEffect[,3],ThreeEffect[,4],FourEffect[,5])
    
    SingleEffect[,2] <- SingleEffect[,2]*meanTwo/meanSingle
    ThreeEffect[,4] <- ThreeEffect[,4]*meanTwo/meanThree
    FourEffect[,5] <- FourEffect[,5]*meanTwo/meanFour
    
    z <- c(SingleEffect[,2],TwoEffect[,3],ThreeEffect[,4],FourEffect[,5])
    
    NmeanSingle <- mean(SingleEffect[,2])
    NmeanTwo <- meanTwo
    NmeanThree <- mean(ThreeEffect[,4])
    NmeanFour <- mean(FourEffect[,5])
    
    # plot
    x.pch <- rep(5,length(x))
    x.pch[1:rowSingle] <- 1
    x.pch[(rowSingle+1):(rowSingle+rowTwo)] <- 2
    x.pch[(rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)] <- 3
    
    x.col <- rep("red",length(x))
    x.col[1:rowSingle] <- "green"
    x.col[(rowSingle+1):(rowSingle+rowTwo)] <- "blue"
    x.col[(rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)] <- "cyan"
    
    dev.new()
    set.seed(8)
    layout(matrix(1:2, 2, 1))
    
    plot(log(x),y,pch=x.pch,col=x.col,xlab=NA,ylab="Original Effect")
    lines(log(1:rowSingle),rep(meanSingle,rowSingle),col="green")
    lines(log((rowSingle+1):(rowSingle+rowTwo)),rep(meanTwo,rowTwo),col="blue")
    lines(log((rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)),
          rep(meanThree,rowThree),col="cyan")
    lines(log((rowSingle+rowTwo+rowThree+1):(rowSingle+rowTwo+rowThree+rowFour)),
          rep(meanFour,rowFour),col="red")
    
    plot(log(x),z,pch=x.pch,col=x.col,xlab=NA,ylab="Normalized Effect")
    lines(log(1:rowSingle),rep(NmeanSingle,rowSingle),col="green")
    lines(log((rowSingle+1):(rowSingle+rowTwo)),rep(NmeanTwo,rowTwo),col="blue")
    lines(log((rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)),
          rep(NmeanThree,rowThree),col="cyan")   
    lines(log((rowSingle+rowTwo+rowThree+1):(rowSingle+rowTwo+rowThree+rowFour)),
          rep(NmeanFour,rowFour),col="red") 
    # savePlot(filename=paste(SaveFileName,"_N",sep=""),type="tiff")
  }
  
  if(MaxOrder==5){
    
    meanSingle <- mean(SingleEffect[,2])
    meanTwo <- mean(TwoEffect[,3])
    meanThree <- mean(ThreeEffect[,4])
    meanFour <- mean(FourEffect[,5])
    meanFive <- mean(FiveEffect[,6])
    
    rowSingle <- nrow(SingleEffect)
    rowTwo <- nrow(TwoEffect)
    rowThree <- nrow(ThreeEffect)
    rowFour <- nrow(FourEffect)
    rowFive <- nrow(FiveEffect)
    
    x <- seq(1, rowSingle+rowTwo+rowThree+rowFour+rowFive)
    y <- c(SingleEffect[,2],TwoEffect[,3],ThreeEffect[,4],FourEffect[,5],FiveEffect[,6])
    
    SingleEffect[,2] <- SingleEffect[,2]*meanTwo/meanSingle
    ThreeEffect[,4] <- ThreeEffect[,4]*meanTwo/meanThree
    FourEffect[,5] <- FourEffect[,5]*meanTwo/meanFour
    FiveEffect[,6] <- FiveEffect[,6]*meanTwo/meanFive
    
    z <- c(SingleEffect[,2],TwoEffect[,3],ThreeEffect[,4],FourEffect[,5],FiveEffect[,6])
    
    NmeanSingle <- mean(SingleEffect[,2])
    NmeanTwo <- meanTwo
    NmeanThree <- mean(ThreeEffect[,4])
    NmeanFour <- mean(FourEffect[,5])
    NmeanFive <- mean(FiveEffect[,6])
    
    # plot
    x.pch <- rep(5,length(x))
    x.pch[1:rowSingle] <- 1
    x.pch[(rowSingle+1):(rowSingle+rowTwo)] <- 2
    x.pch[(rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)] <- 3
    x.pch[(rowSingle+rowTwo+rowThree+1):(rowSingle+rowTwo+rowThree+rowFour)] <- 4
    
    x.col <- rep("red",length(x))
    x.col[1:rowSingle] <- "green"
    x.col[(rowSingle+1):(rowSingle+rowTwo)] <- "blue"
    x.col[(rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)] <- "cyan"
    x.col[(rowSingle+rowTwo+rowThree+1):(rowSingle+rowTwo+rowThree+rowFour)] <- "black"
    
    dev.new()
    set.seed(8)
    layout(matrix(1:2, 2, 1))
    
    plot(log(x),y,pch=x.pch,col=x.col,xlab=NA,ylab="Original Effect")
    lines(log(1:rowSingle),rep(meanSingle,rowSingle),col="green")
    lines(log((rowSingle+1):(rowSingle+rowTwo)),rep(meanTwo,rowTwo),col="blue")
    lines(log((rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)),
          rep(meanThree,rowThree),col="cyan")
    lines(log((rowSingle+rowTwo+rowThree+1):(rowSingle+rowTwo+rowThree+rowFour)),
          rep(meanFour,rowFour),col="black")
    lines(log((rowSingle+rowTwo+rowThree+rowFour+1):(rowSingle+rowTwo+rowThree+rowFour+rowFive)),
          rep(meanFive,rowFive),col="red")
    
    plot(log(x),z,pch=x.pch,col=x.col,xlab=NA,ylab="Normalized Effect")
    lines(log(1:rowSingle),rep(NmeanSingle,rowSingle),col="green")
    lines(log((rowSingle+1):(rowSingle+rowTwo)),rep(NmeanTwo,rowTwo),col="blue")
    lines(log((rowSingle+rowTwo+1):(rowSingle+rowTwo+rowThree)),
          rep(NmeanThree,rowThree),col="cyan")   
    lines(log((rowSingle+rowTwo+rowThree+1):(rowSingle+rowTwo+rowThree+rowFour)),
          rep(NmeanFour,rowFour),col="black") 
    lines(log((rowSingle+rowTwo+rowThree+rowFour+1):(rowSingle+rowTwo+rowThree+rowFour+rowFive)),
          rep(NmeanFive,rowFive),col="red") 
    # savePlot(filename=paste(SaveFileName,"_N",sep=""),type="tiff")
  }
  
  list(SingleEffect=SingleEffect,TwoEffect=TwoEffect,ThreeEffect=ThreeEffect,
       FourEffect=FourEffect,FiveEffect=FiveEffect)
}
