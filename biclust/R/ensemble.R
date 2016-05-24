#### Jaccard Similarity Matrix Function
jaccard2<-function(Rows, Cols)
{
  le<-dim(Rows)[2]
  jaccardmat <- matrix(0,nrow=le,ncol=le)
    for(i in 1:(le-1))
    {
      for(j in (i+1):le)
      {
        #cat(paste("i=",i,"j=",j,"\n"))
        alle1<-Rows[,i] %*% t(Cols[i,])
        alle2<-Rows[,j] %*% t(Cols[j,])
        alle<-alle1 + alle2
        jaccardmat[i,j]<- sum(alle>1)/sum(alle>0)
      }
    }
     jaccardmat
}


### Correlation Similarity Approach

corappr<-function(Rows, Cols, method="both")
{
  le<-dim(Rows)[2]
  corapprmat <- matrix(0,nrow=le,ncol=le)
  rowcor <- cor(Rows, method="pearson")
  colcor <- cor(t(Cols), method="pearson")
    for(i in 1:(le-1))
    {
      for(j in (i+1):le)
      {
        if(method=="both")
        {
            corapprmat[i,j] <- (rowcor[i,j] + colcor[i,j])/2
        }
        else
        {
            if(method=="max")
            {
                corapprmat[i,j] <- min(rowcor[i,j],colcor[i,j])
            }
            else
            {
                corapprmat[i,j] <- max(rowcor[i,j],colcor[i,j])
            }
        }
      }
    }
  corapprmat <- (corapprmat +1) /2
  corapprmat
}



### Grid for Bicluster Algorithms
#Plaid
plaid.grid <- function(method="BCPlaid",cluster="b", fit.model = y ~ m + a + b, background = TRUE, background.layer = NA, background.df = 1, row.release = c(0.5, 0.6, 0.7), col.release = c(0.5, 0.6, 0.7), shuffle = 3, back.fit = 0, max.layers = 20, iter.startup = 5, iter.layer = 10, verbose = FALSE)
{
    resr <- expand.grid(method = method, cluster = cluster, fit.model = deparse(fit.model), background = background, background.layer = background.layer, background.df = background.df, row.release = row.release, col.release = col.release, shuffle = shuffle, back.fit = back.fit, max.layers = max.layers, iter.startup = iter.startup, iter.layer = iter.layer, verbose = verbose, stringsAsFactors=FALSE)

    res <- list()

    for(i in 1:dim(resr)[1])
    {
        res[[i]] <- as.list(resr[i,])
        res[[i]][[3]] <- as.formula(res[[i]][[3]])
    }
    res
}



### Grid for Bimax Function
bimax.grid <- function(method="BCBimax", minr=c(10,11), minc=c(10,11), number=10)
{
  resr <- expand.grid(method = method, minr = minr, minc = minc, number = number, stringsAsFactors = FALSE)

    res <- list()

    for(i in 1:dim(resr)[1])
    {
        res[[i]] <- as.list(resr[i,])
    }
    res
}

### Grid for Xmotifs
xmotif.grid <- function(method="BCXmotifs", ns=200, nd=200,sd=3:6,alpha=c(0.04,0.05, 0.06),number=10)
{
  resr <- expand.grid(method = method, ns=ns,nd=nd,sd=sd,alpha=alpha,number=number, stringsAsFactors = FALSE)

    res <- list()

    for(i in 1:dim(resr)[1])
    {
        res[[i]] <- as.list(resr[i,])
    }
    res
}

### Grid for CC
cc.grid <- function(method="BCCC", delta=c(0.8, 0.9, 1.0), alpha=c(0.8, 0.9, 1.0), number=10)
{
  resr <- expand.grid(method = method, delta = delta, alpha = alpha, number=number , stringsAsFactors = FALSE)

    res <- list()

    for(i in 1:dim(resr)[1])
    {
        res[[i]] <- as.list(resr[i,])
    }
    res
}

### Grid for Spectral
spectral.grid <- function(method="BCSpectral", normalization="log",numberOfEigenvalues=5, minr=c(10,11), minc=c(10,11), withinVar=c(0.1, 0.5, 1))
{
  resr <- expand.grid(method = method, normalization=normalization, numberOfEigenvalues=numberOfEigenvalues, minr = minr, minc=minc, withinVar = withinVar, stringsAsFactors = FALSE)

    res <- list()

    for(i in 1:dim(resr)[1])
    {
        res[[i]] <- as.list(resr[i,])
    }
    res
}



### Ensemble method for Biclustering
ensemble <- function(x, confs, rep = 1, maxNum = 5, similar = jaccard2, thr = 0.8, simthr =0.7, subs = c(1,1), bootstrap = FALSE, support = 0, combine=firstcome, ...)
{
    MYCALL <- match.call()
    le <- length(confs)
    dims <- dim(x)
    bicRow <- matrix(0, dims[1], rep*le*maxNum)
    bicCol <- matrix(0, rep*le*maxNum, dims[2])
    z <- 1
    ### Durch mcapply ersetzen
    for(j in 1:rep)
    {
      sub1 <- sample(1:dims[1], subs[1]*dims[1], replace = bootstrap)
      sub2 <- sample(1:dims[2], subs[2]*dims[2], replace = bootstrap)

      for(i in 1:length(confs))
      {
        res <- do.call("biclust", c(list(x[sub1,sub2]),confs[[i]]))
        ind <- min(res@Number,maxNum)
        if(ind > 0)
        {
          z1 <- z + ind -1
          bicRow[unique(sub1),z:z1] <- res@RowxNumber[!duplicated(sub1), 1:min(res@Number,maxNum)]
          bicCol[z:z1,unique(sub2)] <- res@NumberxCol[1:min(res@Number,maxNum),!duplicated(sub2)]
          z <- z + ind
        }
      }
    }

    bicRow <- bicRow[, 1:(z-1)]
    bicCol <- bicCol[1:(z-1), ]

    com <- combine(bicRow, bicCol, similar=similar, thr = thr, ...)

    ind <- com$ind
    number <- com$number
    counter <- com$counter

    RowxNumber <- matrix(0, dim(x)[1], number)
    NumberxCol <- matrix(0, number, dim(x)[2])


    for(i in 1:number)
    {
        if(counter[i]>1)
        {
            Row <- rowSums(bicRow[,ind[[i]]])
            RowxNumber[,i] <- Row/max(Row)
            Col <- colSums(bicCol[ind[[i]],])
            NumberxCol[i,] <- Col/max(Col)
        }
        else
        {
            Row <- bicRow[,ind[[i]]]
            RowxNumber[,i] <- Row/max(Row)
            Col <- bicCol[ind[[i]],]
            NumberxCol[i,] <- Col/max(Col)
        }


    }

    support <- support * le * rep
    print("Support:")
    print(support)

    if(length(counter)>1)
    {
        counter <- sort(counter, decreasing=TRUE, index.return=TRUE)
        RowxNumber <- RowxNumber[,counter$ix][,counter$x>=support]
        NumberxCol <- NumberxCol[counter$ix,][counter$x>=support,]
        number <- sum(counter$x>=support)
        print("Number of Bicluster:")
        print(counter$x)
    }
    else
    {
        number <- sum(counter>=support)
        print("Number of Bicluster:")
        print(counter)
    }

    if(number==1)
    {
        return(BiclustResult(c(Call=MYCALL,as.list(MYCALL)), matrix(RowxNumber>simthr, ncol=number), matrix(NumberxCol>simthr, nrow=number), number, list(Rowvalues=RowxNumber,Colvalues=NumberxCol, Counts = counter[1])))
    }
    else
    {
        if(number==0)
        {
            return(BiclustResult(c(Call=MYCALL,as.list(MYCALL)) ,matrix(NA,dim(x)[1],1),matrix(NA,1,dim(x)[2]), 0, list(Rowvalues=RowxNumber,Colvalues=NumberxCol, Counts = counter$x)))
        }
        else
        {
            return(BiclustResult(c(Call=MYCALL,as.list(MYCALL)), RowxNumber>simthr, NumberxCol>simthr, number, list(Rowvalues=RowxNumber,Colvalues=NumberxCol, Counts = counter$x)))
        }

    }

}



hcl <- function(bicRow, bicCol, similar=jaccard2, thr = 0.8, ...)
{
      sim <- similar(bicRow, bicCol)
      sim <- sim + t(sim)
      diag(sim) <- 1
      hc <- hclust(as.dist(-sim+1), ...)
      plot(hc)
      hcres <-cutree(hc, h=1-thr)
      hcres <- as.numeric(hcres)
      #print(hcres)
      number <- max(hcres)
      counter <- rep(0,number)
      ind <- list()
      for(i in 1:number)
      {
        ind[[i]] <- hcres==i
        counter[i] <- sum(ind[[i]])
      }
  return(list(ind=ind, counter=counter, number=number))
}

firstcome <- function(bicRow, bicCol, similar=jaccard2, thr = 0.8)
{
 sim <- similar(bicRow, bicCol)
 index <- which(apply(sim,2,max) < thr)
 number <- length(index)
 sim <- sim + t(sim)
 diag(sim) <- 1
 counter <- rep(0,number)
 ind <- list()
 for(i in 1:number)
 {
   ind[[i]] <- sim[,index[i]]>thr
   counter[i] <- sum(ind[[i]])
 }
  return(list(ind=ind, counter=counter, number=number))
}



biggest <- function(bicRow, bicCol, similar=jaccard2, thr = 0.8)
{
 sim <- similar(bicRow, bicCol)
 sim <- sim + t(sim)
 diag(sim) <- 1
 thrsim <- sim > thr
 ind <- list()
 number<-0
 index <- rep(TRUE,dim(sim)[1])
 indexf <- rep(FALSE,dim(sim)[1])
 counter <- c()
 while(sum(index)>1)
 {
     number <- number + 1
     ind[[number]] <- indexf
     ind1 <- as.numeric(which.max(colSums(thrsim[index,index]))[1])
     ind[[number]][index][thrsim[index,index][ind1,]] <- TRUE
     index[ind[[number]]] <- FALSE
     counter <- c(counter,sum(ind[[number]]))
 }
 if(sum(index)==1)
 {
    number <- number + 1
    ind[[number]] <- index
    counter <- c(counter,sum(ind[[number]]))
 }
  return(list(ind=ind, counter=counter, number=number))
}


qtbiggest <- function(bicRow, bicCol, similar=jaccard2, thr = 0.8)
{
 sim <- similar(bicRow, bicCol)
 sim <- sim + t(sim)
 diag(sim) <- 1
 thrsim <- sim > thr
 ind <- list()
 number<-0
 index <- rep(TRUE,dim(sim)[1])
 indexf <- rep(FALSE,dim(sim)[1])
 counter <- c()
 while(sum(index)>1)
 {
     number <- number + 1
     ind[[number]] <- indexf
     ind2 <- colSums(thrsim[index,index])
     ind1 <- sample(1:length(ind2),1,prob=ind2)
     ind[[number]][index][thrsim[index,index][ind1,]] <- TRUE
     index[ind[[number]]] <- FALSE
     counter <- c(counter,sum(ind[[number]]))
 }
 if(sum(index)==1)
 {
    number <- number + 1
    ind[[number]] <- index
    counter <- c(counter,sum(ind[[number]]))
 }
  return(list(ind=ind, counter=counter, number=number))
}








