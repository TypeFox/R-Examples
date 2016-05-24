#
#  Copyright (C) 2007 Sebastian Kaiser
#
#  Prelic, A.; Bleuler, S.; Zimmermann, P.; Wille, A. & , P. A systematic comparison and evaluation of biclustering methods for gene expression data Bioinformatics, Oxford Univ Press, 2006
#
# Calling the C Code from the bimax.c file.

cbimax<- function(logicalmatrix,minr=2,minc=2,number=100,er=0)
   .C("bimax",
   as.integer(logicalmatrix),
   as.integer(nrow(logicalmatrix)),
   as.integer(ncol(logicalmatrix)),
   as.integer(minr),
   as.integer(minc),
   as.integer(matrix(0,nrow=nrow(logicalmatrix),ncol=number)),
   as.integer(matrix(0,nrow=number,ncol=ncol(logicalmatrix))),
   as.integer(vector(mode="integer",length=nrow(logicalmatrix)+ncol(logicalmatrix))),
   as.integer(number),
   as.integer(er))

#### Bimax function for one time running

bimaxbiclust<- function(logicalmatrix,...){
MYCALL<-match.call()
flush.console()
ausgabe<-cbimax(logicalmatrix,...)
#if (ausgabe[[10]]==1)
#{warning("Too many biclusters found change number or minimal dimension of expected bicluster!")
#}
ausgabe[[6]] <- as.logical(ausgabe[[6]])
ausgabe[[7]] <- as.logical(ausgabe[[7]])
RowxNumber<-matrix(ausgabe[[6]],nrow=nrow(logicalmatrix),ncol=ausgabe[[9]])
NumberxCol<-matrix(ausgabe[[7]],nrow=ausgabe[[9]],ncol=ncol(logicalmatrix))
anzahl<-colSums(RowxNumber)
anzahl2<-rowSums(NumberxCol)
anzahl_ges<-anzahl+anzahl2
anzahl_id<-anzahl_ges>0
Number<-sum(anzahl_id)
if(Number==1){
RowxNumber <-matrix(RowxNumber[,anzahl_id],ncol=1)
NumberxCol <-matrix(NumberxCol[anzahl_id,],nrow=1)
}
if(Number>1)
{
RowxNumber<- RowxNumber[,anzahl_id]
NumberxCol <-NumberxCol[anzahl_id,]
}
return(BiclustResult(as.list(MYCALL),RowxNumber,NumberxCol,Number,list(0)))
}


#### Bimax function for repeated running without overlapping


repbimaxbiclust<- function(logicalmatrix,minr=2,minc=2,number=30,maxc=12)
{
  RowxNumber<-matrix(FALSE,nrow=nrow(logicalmatrix),ncol=number)
  NumberxCol<-matrix(FALSE,nrow=number,ncol=ncol(logicalmatrix))
  daten<-logicalmatrix
  datenrows<-rep(TRUE,nrow(logicalmatrix))
  forstop <- TRUE

  for(j in 1:number)
    {
    res1<-0
    res2<-1
    i<-maxc
    k<-0
    while(res1==0 & i>0)
      {
      i<-i-1
      res_bimax <- bimaxbiclust(daten[datenrows,], minr=minr, minc=i, number=100)
      res1<-res_bimax@Number
      #print(i)
      }
    while(res2>0 & i>minc)
      {
      resbic<-res_bimax
      res_bimax <- bimaxbiclust(daten[datenrows,], minr=minr+k, minc=i, number=100)
      k<-k+1
      res2<-res_bimax@Number
      #print(bicluster28_bimax)
      }
    #print(j)
    #print(resbic)
    if(i>minc)
      {
      forstop <- FALSE
      ind<-which.max(colSums(resbic@RowxNumber))
      RowxNumber[datenrows,j]<-resbic@RowxNumber[,ind]
      NumberxCol[j,]<-resbic@NumberxCol[ind,]
      datenrows[datenrows][resbic@RowxNumber[,ind]]<-FALSE
      }
    else
      {
      break
      }

    }
  if(i>minc)
    {

    RowxNumber <-matrix(RowxNumber[,1:j],ncol=j)
    NumberxCol <-matrix(NumberxCol[1:j,],nrow=j)
    bimaxbic<-BiclustResult(resbic@Parameters,RowxNumber,NumberxCol,j,list())
    }
  else
    {
    if(forstop)
      {
      bimaxbic<-BiclustResult(res_bimax@Parameters,matrix(NA,1,1),matrix(NA,1,1),0,list())
      }
    else
      {
      RowxNumber <-matrix(RowxNumber[,1:(j-1)],ncol=(j-1))
      NumberxCol <-matrix(NumberxCol[1:(j-1),],nrow=(j-1))

      bimaxbic<-BiclustResult(resbic@Parameters,RowxNumber,NumberxCol,(j-1),list())
      }
    }
return(bimaxbic)
}


maxbimaxbiclust <- function(logicalmatrix,minr=2,minc=2,number=5, backfit=2, n2=30)
{
  RowxNumber<-matrix(FALSE,nrow=nrow(logicalmatrix),ncol=number)
  NumberxCol<-matrix(FALSE,nrow=number,ncol=ncol(logicalmatrix))
  daten<-logicalmatrix
  datenrows<-rep(TRUE,nrow(logicalmatrix))
  forstop <- TRUE
  sums <- rowSums(logicalmatrix)

  for(j in 1:number)
  {
    print(j)
    i <- minc-1
    k <- 0
    k2 <- FALSE
    size <- 0
    RUN <- TRUE
    FOUND <- FALSE
    i_found <- minc
    while(RUN)
    {
      #print(i)
      i<-i+1
      scountb <- sums >= i
      if(sum(scountb)<max(floor(size/i),minr))
      {
          res1 =0
      }
      else
      {
          res_bimax <- bimaxbiclust(daten[(datenrows&scountb),], minr=max(floor(size/i),minr), minc=i, number=n2)
          res1<-res_bimax@Number
      }

      if(res1==0)
      {
          k <- k + 1
      }
      else
      {
          FOUND <- TRUE
          forstop <- FALSE
          cS <- colSums(res_bimax@RowxNumber)
          rS <- rowSums(res_bimax@NumberxCol)
          sizeb <- max(si <- cS * rS)

          if(sizeb > size)
          {
              #print("da")
              k <- 0
              k2 <- FALSE
              ind <- which.max(si)[1]
              resbic <- res_bimax
              size <- sizeb
              scount <- scountb
              i_found <- i
          }
          else
          {
              k <- k + 1
          }

      }
      #print(k)
      #print(k2)
      if(k>=backfit)
      {
          if(k2)
          {
              RUN <- FALSE
          }
          else
          {
              k2 <- TRUE
              i <- max(i_found - k - 1, minc - 1)
              k <- 0
          }
      }
    }

    if(FOUND)
    {
      RowxNumber[(datenrows&scount),j] <- resbic@RowxNumber[,ind]
      NumberxCol[j,] <- resbic@NumberxCol[ind,]
      datenrows[(datenrows&scount)][resbic@RowxNumber[,ind]]  <- FALSE
    }
    else
    {
      break
    }

  }


  if(FOUND)
    {

    RowxNumber <-matrix(RowxNumber[,1:j],ncol=j)
    NumberxCol <-matrix(NumberxCol[1:j,],nrow=j)
    bimaxbic<-BiclustResult(resbic@Parameters,RowxNumber,NumberxCol,j,list())
    }
  else
    {
    if(forstop)
      {
      bimaxbic<-BiclustResult(res_bimax@Parameters,matrix(NA,1,1),matrix(NA,1,1),0,list())
      }
    else
      {
      RowxNumber <-matrix(RowxNumber[,1:(j-1)],ncol=(j-1))
      NumberxCol <-matrix(NumberxCol[1:(j-1),],nrow=(j-1))

      bimaxbic<-BiclustResult(resbic@Parameters,RowxNumber,NumberxCol,(j-1),list())
      }
    }
return(bimaxbic)
}


