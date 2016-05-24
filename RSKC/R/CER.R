CER <-
  function(ind, true.ind,nob=length(ind))
{
  if (length(ind) != length(true.ind))
    stop("length(ind) != length(true.ind)",length(ind) ,"!=", length(true.ind))
  ind.num <- ind
  true.ind.num <- true.ind
  if (is.character(ind)){
    uniInd <- unique(ind)
    for (i in 1 : length(uniInd))
      {
        ind.num[ind == uniInd[i]] <- i
      }
  }

  if (is.character(true.ind) ){
    uniInd <- unique(true.ind)
    for (i in 1 : length(uniInd))
      {
        true.ind.num[true.ind == uniInd[i]] <- i
      }
    
  }

  rand <- randIndex(as.numeric(ind.num),as.numeric(true.ind.num),correct=FALSE)
  names(rand) <- NULL
  return(1 - rand)
  ##return(sum(abs(one(true.ind)-one(ind)))/choose(nob,2))
}

## one <-
## function(index){
## 	on <- NULL 
## 	c <- combn(index,2)
## 	c <- t(c)
## 	on <- 1*(c[,1]==c[,2])
## 	return(on)
## 	}

norm1 <-
function(y){sum(abs(y))}

norm2 <-
function(y){sqrt(sum(y^2))}


## function specific for the opt digits
## generate bitmap of given observation
showbitmap <-function(index)
  {
    ## data(bitmapMat) ## lazyloading
    ## data(bitmapLab) ## lazyloading
    Nbit=32
    for (iindex in 1 : length(index))
      {
        indivindex<- index[iindex]
        for (ibit in 1 : Nbit)
          {
            cat("\n")
            cat(bitmapMat[[indivindex]][ibit])
          }
        cat("\n obs=",indivindex," true digit=",bitmapLab[indivindex]," \n")
      }
  }


## function declaration

showDigit <- function(index,cex.main=1)
  {
    ## data(DutchUtility) ## lazyloading
    ## 4. DutchUtility-pix: 240 pixel averages in 2 x 3 windows; 
    ## 16 by 15
    ncols = 15 ## replace to 15 
    nrows = 16 ## replace to 16
    labels <- rep(0:9,each=200)
    plot(NA,xlim=c(0,ncols),ylim=c(0,nrows),axes=FALSE,xlab="",ylab="",cex.main=cex.main,
         main=paste("observation",index," True digit",labels[index],sep=""))
    abline(h=0:ncols,v=0:nrows,lty=2,col="gray70")
    axis(1,(1:ncols)-0.5,1:ncols,lty=0,cex=0.5)
    axis(2,(1:nrows)-0.5,nrows:1,padj=1,lty=0,cex=0.5)
    ##
    cols <- gray.colors(n=6,start=0.9,end=0.3)
    for (iy in 1 : nrows)
      {
        for (ix in 1 : ncols)
          {
            ## each row vector of a matrix DutchUtility contains the bitmap X; 
            ## DutchUtility[1,] = c(X[1,],X[2,],...) 
            Pickedcolor <- cols[DutchUtility[index,ix+ncols*(iy-1)]]
            polygon(x=c(ix-1,ix-1,ix,ix),
                    y=c(nrows-iy,nrows-iy+1,nrows-iy+1,nrows-iy),
                    col=Pickedcolor,
                    border=FALSE)
          }
      }
  }


      
## sensitivity the other way around for the first example
## sensitivity for the second example both way

Sensitivity <- function (label1, label2
                         ##, Alpha=FALSE, which.Alpha="label1"
                         )
  {
    ## Given two partitions, label1 and label2, compute its "sensitivity"
    ## max{ (label1[label2==k])/sum(label2==k) }, and the corresponding label of label1 that achieves max 
    uni1 <- sort(unique(label1))
    uni2 <- sort(unique(label2))

    ## if (length(uni1) >= 27 & Alpha & which.Alpha=="label1")
    ##   stop("The alphabet option is valid only if the number of clusters from label1 is less than 27")
    ## if (length(uni2) >= 27 & Alpha & which.Alpha=="label2")
    ##   stop("The alphabet option is valid only if the number of clusters from label2 is less than 27")
    
    ## ## change the labels of clusters in which.Alpha
    ## if (Alpha & which.Alpha=="label1"){
    ##   temp <- label1
    ##   for (iuniL in 1 : length(uni1))
    ##     {
    ##       temp[label1==uni1[iuniL]] <- letters[iuniL]
    ##     }
    ##   label1 <- temp
    ##   uni1 <- sort(unique(label1))
      
    ## }else if (Alpha & which.Alpha=="label2"){
    ##   temp <- label2
    ##   for (iuniL in 1 : length(uni1))
    ##     {
    ##       temp[label2==uni2[iuniL]] <- letters[iuniL]
    ##     }
    ##   label2 <- temp
    ##   uni2 <- sort(unique(label2))
    ## }
    
    K <- length(uni2)
    senst <- correspondClass <- rep(NA, K)
    ## this will return table both columns and rows are in sorted order.
    tbl <- table(label1, label2)
    
    trueTot <- colSums(tbl)
    ## scale each column of the table matrix by its column sum
    ##  colSums(prMat)= 1,1,...
    prMat <- scale(tbl, scale = trueTot, center = FALSE)
    senst <- apply(prMat, 2, max)*100

    correspondClass <- rownames(prMat)[apply(prMat, 2, which.max)]
    senst <- sprintf("%1.0f",senst)
    
    re <- data.frame(rbind(senst, correspondClass))
    
    names(re) <- uni2
    rownames(re) <- c("Sensitivity. (%)", "Class label by label1.")
    return(list(prob=re,table=tbl,marginal=prMat))
  }
