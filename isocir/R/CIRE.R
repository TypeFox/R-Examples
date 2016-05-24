"CIRE"<-function(data, groups=c(1:nrow(data)), circular=TRUE){

# DATE WRITTEN: 04 Mar 2010          LAST REVISED:  20 Feb 2012
# AUTHOR: Sandra Barragan based on the SAS routines written by Miguel A. Fernandez.
# DESCRIPTION: Computes the algorithm to obtain the CIRE.
# REFERENCE: Rueda, C., Fernandez, M. A. and Shyamal, P. (2009). 
#             Estimation of parameters subject to order restrictions 
#             on a circle with application to estimation of phase angles of cell-cycle genes. 
#             JASA.
# SEE ALSO: cond.test, CIRE.


##################################################
##  auxiliary functions needed for obCIRE:
##################################################

 cirmean<-function (data) {
    data <- data[!is.na(data)]
    A <- sum(sin(data))
    B <- sum(cos(data))
    if ((A >= 0) && (B > 0)) {
        result <- atan(A/B)
    }
    if (B < 0) {
        result <- atan(A/B) + pi
    }
    if ((A < 0) && (B > 0)) {
        result <- atan(A/B) + (2 * pi)
    }
    return(result)
}

mmeans <- function(data){
   means <- matrix(0, nrow=length(data), ncol=length(data))
   for(l in 1:length(data)){
    for(m in l:length(data)){
     means[l,m]<-cirmean(data[l:m])
    }
   }
   means<-means+t(means)
   diagonal<-diag(diag(means),nrow=nrow(means),ncol=ncol(means))
   means<-means-(diagonal/2)
   return(means)
   # the result is a matrix nxn where n=length(data)
}

cirPAVA <- function(data){
   means <- mmeans(data)
   data1 <- data
   problems <- ifelse(data[1:length(data) - 1] - data[2:length(data)] > 0, 1, 0)
   while(sum(problems) != 0){
    pos <- list()
    k <- 1
    for(i in 1:length(problems)){
     if(problems[i] == 1){
      if((i != 1) && problems[i - 1] == 0){pos[[k]] <- vector()}
      if(i == 1){pos[[k]] <- vector()}
      pos[[k]] <- c(pos[[k]], c(1:length(problems))[i])
     }
     if((problems[i] == 0) && (i != 1)){
      if((problems[i - 1] == 1) && (i != length(problems))){
       k <- k + 1
      }
     }
    } # close for(i)
    for(j in 1:length(pos)){
     a <- min(pos[[j]])
     a <- min(c(1:length(data))[data == data[a]])
     b <- max(pos[[j]])+1
     b <- max(c(1:length(data))[data==data[b]])
     data[a:b] <- means[a, b]
    } # close for(j)
    problems <- ifelse(data[1:length(data) - 1] - data[2:length(data)] > 0, 1, 0)
   } # close while
   SCE <- sum(1 - cos(data1 - data))
   return(list(data = data, SCE = SCE))
}

funCIRE <- function(data, means, mrl){

   data1 <- data
   q <- length(data)
   lA1.PAVA <- list()
   lA2.PAVA <- list()
   lA3.PAVA <- list()
   resultlist<-.Call("funCIRE", data, means, mrl, PACKAGE="isocir")
   lA1.PAVA<-resultlist[[1]]
   lA2.PAVA<-resultlist[[2]]
   lA3.PAVA<-resultlist[[3]]
  
# #applying PAVA 
#    lA1.PAVA<-list()
#    SCE1<-array()
#    for(i in 1:length(lA1)){
#      if(!is.null(lA1[[i]][1])){ 
#      lA1.PAVA[[i]]<-cirPAVA(lA1[[i]][1,])$data 
#      SCE1[i]<-sce(data[1:i],lA1.PAVA[[i]],mrl[1:i])
#      lA1[[i]]<-unique(lA1[[i]], MARGIN=1)
#      if(nrow(lA1[[i]])>1){  
#        for(n in 2:nrow(lA1[[i]])){
#          result1<-cirPAVA(lA1[[i]][n,])$data
#          auxSCE<-sce(data[1:i],result1,mrl[1:i])
#          if(result1[length(result1)]<(pi/2)){ # condition (5)
#            if(auxSCE<SCE1[i]){
#              lA1.PAVA[[i]]<-result1
#              SCE1[i]<-auxSCE
#            }# close if<SCE
#          }# close if condition (5)
#        }# close for solutions
#      }# close if solutions
#      if(lA1.PAVA[[i]][length(lA1.PAVA[[i]])]>(pi/2)){lA1.PAVA[[i]]<-rep(10000,length(lA1.PAVA[[i]]))}
#      }# close check NA
#      if(is.null(lA1[[i]][1])){lA1.PAVA[[i]]<-NA}
#    }# close for(lA1) 
#    
#    lA2.PAVA<-list()
#    SCE2<-matrix(ncol=length(lA2), nrow=length(lA2))
#    for(i in 1:length(lA2)){
#      lA2.PAVA[[i]]<-list()
#      for(j in i:length(lA2[[i]])){
#        if(!is.null(lA2[[i]][[j]][1])){
#          means2<-means[i:j,i:j]
#          lA2.PAVA[[i]][[j]]<-cirPAVA(lA2[[i]][[j]][1,])$data
#          SCE2[i,j]<-sce(data[i:j],lA2.PAVA[[i]][[j]],mrl[i:j])
#          lA2[[i]][[j]]<-unique(lA2[[i]][[j]], MARGIN=1)
#          if(nrow(lA2[[i]][[j]])>1){
#            for(n in 2:nrow(lA2[[i]][[j]])){
#              result2<-cirPAVA(lA2[[i]][[j]][n,])$data
#              auxSCE<-sce(data[i:j],result2,mrl[i:j])
#              if(result2[1]>(pi/2)&&result2[length(result2)]<(3*pi/2)){ # conticion (5)
#                if(auxSCE<SCE2[i,j]){
#                  lA2.PAVA[[i]][[j]]<-result2
#                  SCE2[i,j]<-auxSCE
#                }# close if<SCE
#              }# close if conticion (5)
#            }# close for solutions
#          }# close if solutions
#          if(lA2.PAVA[[i]][[j]][length(lA2.PAVA[[i]][[j]])]>(3*pi/2)||lA2.PAVA[[i]][[j]][1]<(pi/2)){
#            lA2.PAVA[[i]][[j]]<-rep(10000,length(lA2.PAVA[[i]][[j]]))
#          }
#        }# close check null
#        if(is.null(lA2[[i]][[j]][1])){lA2.PAVA[[i]][[j]]<-NA}
#      }# close for(j)
#    }# close for(i)
# 
#    lA3.PAVA<-list()
#    SCE3<-array()
#    for(i in 1:length(lA3)){
#    if(!is.null(lA3[[i]][1])){
#      means3<-means[i:nrow(means),i:ncol(means)]
#      lA3.PAVA[[i]]<-cirPAVA(lA3[[i]][1,])$data
#      SCE3[i]<-sce(data[i:q],lA3.PAVA[[i]],mrl[i:q])
#      lA3[[i]]<-unique(lA3[[i]], MARGIN=1)
#      if(nrow(lA3[[i]])>1){
#        for(n in 1:nrow(lA3[[i]])){
#        result3<-cirPAVA(lA3[[i]][n,])$data
#        auxSCE<-sce(data[i:q],result3,mrl[i:q])
#        if(result3[1]>(3*pi/2)){  # conticion (5)
#          if(auxSCE<SCE3[i]){
#            lA3.PAVA[[i]]<-result3
#            SCE3[i]<-auxSCE
#          }# close if<SCE
#        }# close if conticion (5)
#        }# close for (n)
#      }# close if solutions
#      if(lA3.PAVA[[i]][1]<(3*pi/2)){lA3.PAVA[[i]]<-rep(10000,length(lA3.PAVA[[i]]))}
#    }# close check NA
#    if(is.null(lA3[[i]][1])){lA3.PAVA[[i]]<-NA}
#    }# close for(i)
#    
   
   fis.tilde<-matrix(ncol=length(data), nrow=(length(data)+1)*length(data))
   z<-1
   #fis.tilde[z,]<-c(data,0)
   #z<-z+1
   fis.tilde[z,]<-lA3.PAVA[[1]]
   z<-z+1
   for(j in 2:length(data)){
     fis.tilde[z,]<-c(lA2.PAVA[[1]][[j-1]],lA3.PAVA[[j]])
     z<-z+1     
   }
   fis.tilde[z,]<-lA2.PAVA[[1]][[length(data)]]
   z<-z+1
   for(i in 2:length(data)){
     for(j in i:length(data)){
       if(j!=length(data)){
         fis.tilde[z,]<-c(lA1.PAVA[[i-1]],lA2.PAVA[[i]][[j]],lA3.PAVA[[j+1]])
         z<-z+1       
       }
       if(j==length(data)){
         fis.tilde[z,]<-c(lA1.PAVA[[i-1]],lA2.PAVA[[i]][[j]])
         z<-z+1         
       }
     }     
   }
   for(i in 1:(length(data)-1)){
     fis.tilde[z,]<-c(lA1.PAVA[[i]],lA3.PAVA[[i+1]])
     z<-z+1
   }
   fis.tilde[z,]<-c(lA1.PAVA[[length(data)]])
   fis.tilde<-fis.tilde[1:z,]
   addition<-apply(fis.tilde,1,sum)
   fis.tilde<-fis.tilde[addition<10000,]
   fis.tilde<-unique(fis.tilde,MARGIN=1)
   SCE<-apply(fis.tilde,1,sce,data,mrl)
   CIRE<-fis.tilde[which.min(SCE),]
   SCE<-min(SCE)
   rm(fis.tilde,lA1.PAVA,lA2.PAVA,lA3.PAVA)
return(list(CIRE=CIRE,SCE=SCE))
}# end function funCIRE

####################

##################################################
##  control of arguments:
##################################################

if(missing(data)){stop("Data cannot be misssing")}
if(all(c(!is.matrix(data),!is.vector(data),!is.data.frame(data),
         !(is.circular(data)&&is.numeric(data))))){
      stop("The argument data must be a matrix or a vector")}
if(!is.matrix(data)){
  if(is.vector(data)){data<-cbind(data)}
  if(is.data.frame(data)){data<-as.matrix(data)}
}
#if(!is.circular(data)){data <- suppressWarnings(as.circular(data))}
if(length(data) == 1){stop("Data cannot be of length 1")}

means <- apply(data,1,cirmean)
resultant <- mrl(data)

n <- as.vector(table(groups))
b<-as.numeric(names(table(groups)))
if(length(b)!=b[length(b)]){stop("There is some missing level")}
k <- length(n)
listdata <- list()
length(listdata) <- k
for(i in 1:length(groups)){
  if(!is.null(listdata[[groups[i]]])){listdata[[groups[i]]] <- c(listdata[[groups[i]]],means[i])}
  if(is.null(listdata[[groups[i]]])){listdata[[groups[i]]] <- means[i]}
}
ordresult <- groups
names(ordresult) <- resultant
sortresult <- as.numeric(names(sort(ordresult)))


# all necessary data:
ldata <- listdata
rm(listdata)
mrl <- sortresult
num <- n
nr <- prod(factorial(num))
nc <- sum(num)
n <- factorial(num)


# create mdata, the matrix with all the possibilities
# dim>1 just in the case of the partial order
indexes<-matrix(ncol=nc,nrow=nr)
mdata <- matrix(ncol=nc,nrow=nr)
aux <- 1
for(p in 1:length(ldata)){
 l <- length(ldata[[p]])
 if(l == 1){
  mdata[, aux] <- rep(ldata[[p]], nr)
  indexes[, aux] <- rep(aux, nr)
  aux <- aux + 1 
  }# end if level of size 1
 if(l > 1){
  iaux <- c(aux:(aux + l - 1))
  permu <- permn(ldata[[p]])
  ipermu <- permn(iaux)
  lp <- length(permu)
  nrep <- nr / lp
  if(p == 1){
   maux <- matrix(nrow = lp, ncol = l)
   miaux <- matrix(nrow = lp, ncol = l)
   for(d in 1:lp){
    maux[d,] <- permu[[d]]
    miaux[d,] <- ipermu[[d]]
    }
   for(c in 1:nrep){
    mdata[(n[1]*(c - 1) + 1):(c*n[1]), 1:l] <- maux
    indexes[(n[1]*(c - 1) + 1):(c*n[1]), 1:l] <- miaux
    }
   aux <- aux+l
   } # close if p=1
  if(p != 1){
   repe <- prod(n[1:(p - 1)])
   maux <- matrix(nrow = repe, ncol = l)
   miaux <- matrix(nrow = repe, ncol = l)
   for(a in 1:(nr / repe)){
    aux1 <- a%%n[p]
    if(aux1 == 0){aux1 <- n[p]}
    repe <- prod(n[1:(p - 1)])
    maux <- matrix(nrow = repe, ncol = l)
    miaux <- matrix(nrow = repe, ncol = l)
    for(d in 1:repe){
     maux[d,] <- permu[[aux1]]
     miaux[d,] <- ipermu[[aux1]]
     }
    if(a > n[p]){aux1 <- aux1 + 1}
    mdata[(((a - 1)*repe) + 1):(a*repe), aux:(aux + l - 1)] <- maux
    indexes[(((a - 1)*repe) + 1):(a*repe), aux:(aux + l - 1)] <- miaux
    } # close for a
   aux <- aux + l
   } # close if p!=1
  } # close if l>1
 } # close p



CIRE <- matrix(ncol = ncol(mdata), nrow = nrow(mdata))
SCE <- array()

# calculate the CIRE for each row of mdata
for(m in 1:nrow(mdata)){
 data <- mdata[m,]
 mrl2 <- mrl
 ind <- indexes[m,]
 means <- mmeans(data)
 resultfunCIRE <- funCIRE(data, means, mrl2)
 CIRE[m,] <- resultfunCIRE$CIRE
 SCE[m] <- resultfunCIRE$SCE
 rm(resultfunCIRE)
 # In case of circular order:
 if(circular == TRUE){
  aux <- 1
  for (i in 2:length(data)){
   data <- c(data[2:length(data)], data[1])
   mrl2 <- c(mrl2[2:length(mrl)], mrl2[1])
   ind <- c(ind[2:length(ind)], ind[1])
   means <- mmeans(data)
   resultfunCIRE <- funCIRE(data, means, mrl2)
   # If the SCE is better, we take that new result
   if(resultfunCIRE$SCE < SCE[m]){
    CIRE[m,] <- resultfunCIRE$CIRE
    SCE[m] <- resultfunCIRE$SCE
    rm(resultfunCIRE)
    indexes[m,] <- ind
    aux <- i
    } # close if<SCE
   } # close for(i)
  } # close circular=TRUE
 }# close for(m)

rm(mdata, mmeans)

CIRESfinals <- CIRE[c(1:length(SCE))[SCE == min(SCE)],]
if(is.matrix(CIRESfinals)){eCIRE <- as.vector(unique(CIRESfinals))}
if(!is.matrix(CIRESfinals)){eCIRE <- CIRESfinals}
minSCE <- unique(SCE[SCE == min(SCE)])
finalindexes <- indexes[c(1:length(SCE))[SCE == min(SCE)],]
if(is.matrix(finalindexes)){names(eCIRE) <- finalindexes[1,]}
if(!is.matrix(finalindexes)){names(eCIRE) <- finalindexes}

rm(indexes, CIRESfinals)

sorted <- as.numeric(names(eCIRE))
names(sorted) <- eCIRE
sorted <- as.numeric(names(sort(sorted)))
cumform <- cumsum(num)
lCIRE <- list()
length(lCIRE) <- length(cumform)
lCIRE[[1]] <- sorted[1:cumform[1]]
if(length(cumform) > 1){
 for(i in 2:length(cumform)){
  lCIRE[[i]] <- sorted[((cumform[(i - 1)]) + 1):cumform[i]] 
  }
 }

lexit <- list(cirmeans = ldata, SCE = minSCE, CIRE = lCIRE)

class(lexit) <- "isocir"
return(lexit)

} ###### end function CIRE





