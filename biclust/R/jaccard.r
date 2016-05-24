#jaccardind_old<-function(bicres1,bicres2){
#  alle1<-bicres1@RowxNumber %*% bicres1@NumberxCol
#  alle2<-bicres2@RowxNumber %*% bicres2@NumberxCol
#  alle<-alle1 + alle2
#  loalle<-alle>0
#  loalle1<-alle1>0
#  loalle2<-alle2>0
#  res<- (sum(loalle1)+sum(loalle2)-sum(loalle))/sum(loalle)
#
#  res
#
#
#}

jaccardind <- function(bicres1, bicres2){
	jaccard1(bicres1, bicres2)/max(jaccard1(bicres1, bicres1),jaccard1(bicres2, bicres2))
}


jaccard1<-function(bicres1,bicres2){
	le1 <- bicres1@Number
	le2 <- bicres2@Number
	jacvec <- c()
	
	for (i in 1:le1) {
		jacvec2 <- c()
		for (j in 1:le2) {
			
			# Taking into account the case when NUMBER=1 and the output is in WRONG FORMAT
			#	Note on 'Wrong Format': This occurs only when the expected BC for the algorithm is >1, but only 1 BC is discovered. (NOT when you are asking an algorithm for 1 BC)
			#		Why? as.matrix() of a vector results into a matrix with 1 column always
			if(le1==1 & (dim(bicres1@NumberxCol)[2]==1)  ){      
				res1.NumberxCol <- bicres1@NumberxCol[,1]
			}
			else{
				res1.NumberxCol <- bicres1@NumberxCol[i,]
			}
			
			if(le2==1 & (dim(bicres2@NumberxCol)[2]==1)){
				res2.NumberxCol <- bicres2@NumberxCol[,1]
				
			}
			else{
				res2.NumberxCol <- bicres2@NumberxCol[j,]
			}
			
			
			
			alle1 <- bicres1@RowxNumber[, i] %*% t(res1.NumberxCol)
			alle2 <- bicres2@RowxNumber[, j] %*% t(res2.NumberxCol)
			alle <- alle1 + alle2
			loalle <- alle > 0
			loalle1 <- alle1 > 0
			loalle2 <- alle2 > 0
			jacvec2 <- c(jacvec2,((sum(loalle1) + sum(loalle2) -sum(loalle))/sum(loalle))) 
		}
		jacvec <- c(jacvec, sum(jacvec2,na.rm=TRUE))
	}
	
	res <- sum(jacvec)/max(le1, le2)
	return(res)
}


#jaccardind_vec<-function(bicres1,bicres2){
#le1<-bicres1@Number
#le2<-bicres2@Number
#jacvec<-matrix(0,le1,le2)
#
#for (i in 1:le1)
# {
#  for (j in 1:le2)
#  {
#    alle1<-bicres1@RowxNumber[,i] %*% t(bicres1@NumberxCol[i,])
#    alle2<-bicres2@RowxNumber[,j] %*% t(bicres2@NumberxCol[j,])
#    alle<-alle1 + alle2
#    loalle<-alle>0
#    loalle1<-alle1>0
#    loalle2<-alle2>0
#    jacvec[i,j] <- (sum(loalle1)+sum(loalle2)-sum(loalle))/sum(loalle)
#  }
#}
#jacvec
#
#}



