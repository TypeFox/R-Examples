heatmapBC <- function(x,bicResult,number=0,local=FALSE,order=FALSE, axes = FALSE,outside = FALSE, zlim = c(min(x), max(x)), ...){
    n <- bicResult@Number


    if(is.numeric(number)){
        if(length(number)==1){
              if(number==0){
                  bicRows <- numeric()
                  bicCols <- numeric()
                  xl <- numeric(n)
                  xr <- numeric(n)
                  yt <- numeric(n)
                  yb <- numeric(n)
                  xlo <- numeric(n)
                  yto <- numeric(n)
                  res <- list()

                  for(i in 1:n){
                      res <- heatorder(x, bicResult, bicRows, bicCols, order, i, n, i+1)
                      bicRows <- res[[1]]
                      bicCols <- res[[2]]
                      xl[i] <- res[[3]]
                      xr[i] <- res[[4]]
                      yt[i] <- res[[5]]
                      yb[i] <- res[[6]]
                      xlo[i] <- res[[7]]
                      yto[i] <- res[[8]]

                  }

                  bicRows <- c(bicRows, which(!(1:dim(x)[1] %in% bicRows)))
                  bicCols <- c(bicCols, which(!(1:dim(x)[2] %in% bicCols)))

                  image(t(x[rev(bicRows),bicCols]), axes=axes, zlim=zlim, x = 1:length(bicCols), y = 1:length(bicRows), ...)

                  rect(xleft = xl[1]+0.5, xright = xr[1] + 0.5, ytop = length(bicRows) - yt[1] + 0.5, ybottom = length(bicRows) - yb[1] + 0.5)
                  for(i in 2:n){
                      rect(xleft = xl[i]- xlo[i-1]+0.5, xright = xr[i] + 0.5, ytop = length(bicRows) - yt[i] + yto[i-1] + 0.5, ybottom = length(bicRows) - yb[i] + 0.5)
                  }

                  if(outside){
                      print("Hallo")
                      overR <- which(bicResult@RowxNumber[,1])
                      overC <- which(bicResult@NumberxCol[1,])
                      res <- list()
                      for(i in 2:n){
                          res <- printrect(x, bicResult, overR, overC, i, i-1, order, bicRows, bicCols, xl, xlo, yt, yto, xr, yb)
                          overR <- res[[1]]
                          overC <- res[[2]]
                      }
                  }


              }
              else {
                  bicRows=which(bicResult@RowxNumber[,number])
                  yb <- length(bicRows)
                  bicCols=which(bicResult@NumberxCol[number,])
                  xr <- length(bicCols)

                  if(order)
                  {
                      bicRows <- bicRows[order(rowSums(x[bicRows,bicCols]))]
                      bicCols <- bicCols[order(colSums(x[bicRows,bicCols]))]
                  }

                  if(!local){
                      bicRows <- c(bicRows, which(!(1:dim(x)[1] %in% bicRows)))
                      bicCols <- c(bicCols, which(!(1:dim(x)[2] %in% bicCols)))
                  }
                  image(t(x[rev(bicRows),bicCols]), axes=axes, zlim=zlim, x = 1:length(bicCols), y = 1:length(bicRows), ...)

                  if(!local){
                      rect(xleft = 0.5, xright = xr + 0.5, ytop = length(bicRows)+0.5, ybottom = length(bicRows) - yb + 0.5)
                  }

              }

          }
        if(length(number)>1){
                  bicRows <- numeric()
                  bicCols <- numeric()
                  xl <- numeric(n)
                  xr <- numeric(n)
                  yt <- numeric(n)
                  yb <- numeric(n)
                  xlo <- numeric(n)
                  yto <- numeric(n)
                  res <- list()

                  for(i in 1:length(number)){
                      res <- heatorder(x, bicResult, bicRows, bicCols, order, number[i], number[length(number)], number[min(i+1,length(number))])
                      bicRows <- res[[1]]
                      bicCols <- res[[2]]
                      xl[i] <- res[[3]]
                      xr[i] <- res[[4]]
                      yt[i] <- res[[5]]
                      yb[i] <- res[[6]]
                      xlo[i] <- res[[7]]
                      yto[i] <- res[[8]]
                  }
                  image(t(x[rev(bicRows),bicCols]), axes=axes, zlim=zlim, x = 1:length(bicCols), y = 1:length(bicRows),...)

                  rect(xleft = xl[1]+0.5, xright = xr[1] + 0.5, ytop = length(bicRows) - yt[1] + 0.5, ybottom = length(bicRows) - yb[1] + 0.5)
                  for(i in 2:length(number)){
                      rect(xleft = xl[i]- xlo[i-1]+0.5, xright = xr[i] + 0.5, ytop = length(bicRows) - yt[i] + yto[i-1] + 0.5, ybottom = length(bicRows) - yb[i] + 0.5)
                  }




        }


     }else{
        image(t(x),axes=axes,zlim=zlim,...)
    }





}

printrect <- function(x, bicResult, overR1, overC1, number, before, order, bicRows, bicCols, xl, xlo, yt, yto, xr, yb){
    overR <- which(bicResult@RowxNumber[,number])
    overC <- which(bicResult@NumberxCol[number,])

    overR <- overR[(overR %in% overR1)]
    overC <- overC[(overC %in% overC1)]
    R1 <- c(overR[!(overR %in% overR1)],overR1)
    C1 <- c(overC[(overC %in% overC1)], overC1)

    if(order){
        print("order")
        print(overR)
        if(length(overR)>0){
            for (j in 1:length(overR)){
                a<-which(overR[j] == bicRows)
                rect(xleft = xl[number]- xlo[before]+0.5, xright = xr[number] + 0.5, ytop = length(bicRows) - a  + 0.5, ybottom = length(bicRows) - a - 0.5)
            }
        }
        if(length(overC)>0){
            for (j in 1:length(overC)){
                b <- which(overC[j] == bicCols)
                rect(xleft = b+0.5, xright = b-0.5, ytop = yt[number]  + 0.5, ybottom = length(bicRows) - yb[number] + 0.5)
            }
        }

    } else {
        print("nichtorder")
        if(length(overR)>0){
            overR <- overR[!(overR %in% which(bicResult@RowxNumber[,before]))]
            print(overR)
            for (j in 1:length(overR)){
                a <- which(overR[j] == bicRows)
                rect(xleft = xl[number]- xlo[before]+0.5, xright = xr[number] + 0.5, ytop = length(bicRows) - a  + 0.5, ybottom = length(bicRows) - a - 0.5)
            }
        }
        if(length(overC)>0){
            overC <- overC[!(overC %in% which(bicResult@NumberxCol[before,]))]
            for (j in 1:length(overC)){
                b <-which(overC[j] == bicCols)
                rect(xleft = b+0.5, xright = b-0.5, ytop = yt[number] + yto[before] + 0.5, ybottom = length(bicRows) - yb[number] + 0.5)
            }
        }
    }

list(overR=R1, overC=C1)
}



heatorder <- function(x, bicResult, bicRows1, bicCols1, order, number, end, bicnext){
    xl <- length(bicCols1)
    yt <- length(bicRows1)

    bicRows <- which(bicResult@RowxNumber[,number])
    bicCols <- which(bicResult@NumberxCol[number,])
    if(order)
    {
        bicRows <- bicRows[order(rowSums(x[bicRows,bicCols]))]
        bicCols <- bicCols[order(colSums(x[bicRows,bicCols]))]
    }

    bicRows <- bicRows[!(bicRows %in% bicRows1)]
    bicCols <- bicCols[!(bicCols %in% bicCols1)]

    if(!(number == end) && !order){

            bicRows2 <- which(bicResult@RowxNumber[,bicnext])
            bicRows <- c(bicRows[!(bicRows %in% bicRows2)],
                          bicRows[ (bicRows %in% bicRows2)])
            yto <- sum((bicRows %in% bicRows2))
            print(paste("yto",yto))
            bicCols2 <- which(bicResult@NumberxCol[bicnext,])
            bicCols <- c(bicCols[!(bicCols %in% bicCols2)],
                          bicCols[ (bicCols %in% bicCols2)])

            xlo <- sum((bicCols %in% bicCols2))
                  print(paste("xlo",xlo))
        }else {
            yto <- 0
            xlo <- 0

        }

    list(bicRows=c(bicRows1,bicRows), bicCols=c(bicCols1, bicCols), xl=xl, xr=length(c(bicCols1,bicCols)), yt=yt, yb=length(c(bicRows1,bicRows)), xlo=xlo, yto=yto , bicCols, bicRows)
}




