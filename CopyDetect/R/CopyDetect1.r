CopyDetect1 <- function(data,item.par=NULL,pair) { #Start Function

#########################################################################################################

	# data,  N x n response data
		# "data" should be a data frame object
		# N is number of examinee and n is number of items
		# Binary responses 

	# item.par, n x 3 item parameter matrix, in the logistic metric
		#first column is item discrimination parameters
		#second column is item difficulty parameters
		#third column is item guessing parameters

		#If the user does not supply an item parameter matrix, item parameters
		#are estimated using ltm-package in R based on 2 PL IRT Model.

	# pair, a vector of length two. 

		#the first element of the vector is the row number of the suspected copier examinee in 
		 #in the data file
		#the second element of the vector is the row number of the source examinee for the 
		 #suspected copier examinee in the data file.
	
	#CopyDetect computes w, GBT, S2, S1, K, K1, and K2 indices internally, and print the information
		#regarding the likelihood of agreement between suspected copier and source examinees'
		#response vector.

#########################################################################################################


	row.names <- c("Examinee1");for(i in 2:nrow(data)){row.names <- c(row.names,paste("Examinee",i,sep="")) }
	rownames(data) <- row.names



#Check item parameter file, if supplied check if that is appropriate, if not supplied
	#use ltm package to estimate 2PL item parameters


	if(is.null(item.par)==TRUE) { item.par <- est(data, model = "2PL", engine = "ltm")$est
                                  }


#Other checks for inputs

		if(ncol(item.par)!=3){ stop("The item parameter file should have three columns. Check your input matrix")}
		if(nrow(item.par)!=ncol(data)) {stop("The number of rows in item parameter file should be equal to the number of columns in the response data file. Check your input matrices")}
		if(is.data.frame(data)!=TRUE)  {stop("The input response data file is not a data frame object. Use as.data.frame() first to make your response data file a data frame object")}
		if(length(pair)!=2){stop("The input pair should be a vector of length two. The first element indicates the row number in the response data file for the suspected copier examinee, and the second element indicates the row number in the response data file for the suspected source examinee")}
                if(pair[1]>nrow(data) | pair[2] > nrow(data)) { stop("The elements of the input vector pair can not be bigger than the number of rows in the response data file")}
		
#Internal function to compute probability correct for a binary IRT model

		irtprob <- function(param,theta,d=1){ 
 				   return( param[,3] + (1 - param[,3])/( 1 + exp(-d*param[,1]*(theta-param[,2]))))
            		    }

######################             Computing w index       #########################################################

	omega <- function(form,ip,pa) { #start internal function
	
		theta.est1 <- mlebme(form[pa[1],],ip)[1]       
 
		common.missing <- which(is.na(form[pa[1],])==TRUE & is.na(form[pa[2],])==TRUE)
            form[pa[1],common.missing]=0
            form[pa[2],common.missing]=0
	
		obs.match <- length(which(form[pa[1],]==form[pa[2],]))


	      prob.cor <- cbind(irtprob(param=ip,theta=theta.est1),1-irtprob(param=ip,theta=theta.est1))

			colnames(prob.cor) <- c(paste("Probability Correct for Examinee ",pa[1],sep=""),
							paste("Probability Incorrect for Examinee ",pa[1],sep=""))		
			row.names <- c("Item 1");for(i in 2:ncol(form)){row.names <- c(row.names,paste("Item ",i,sep="")) }
			rownames(prob.cor) <- row.names

				pvec <- c(prob.cor[which(form[pa[2],]==1),1],
                                 (1-prob.cor[which(form[pa[2],]==0),1]),
					   (1-prob.cor[which(is.na(form[pa[2],])==TRUE),1])
					   )

		exp.match <- sum(pvec)
		sd.match  <- sqrt(sum(pvec*(1-pvec)))

		w.value <- (obs.match-exp.match)/sd.match 
		p.value <- pnorm(w.value,0,1,lower.tail=FALSE)
		
		return(list(ability1=theta.est1,obs.match=obs.match,
				exp.match=exp.match,
				sd.match=sd.match,
				W.value=w.value,
			 	p.value=p.value))
	}#end internal function

	w.index <- omega(form=data,ip=item.par,pa=pair)


##################                     Computing GBT index              ##################################################


	GBT <- function(form,ip,pa) { #start internal function
	
		theta.est1 <- mlebme(form[pa[1],],ip)[1]  
		theta.est2 <- mlebme(form[pa[2],],ip)[1]       

		common.missing <- which(is.na(form[pa[1],])==TRUE & is.na(form[pa[2],])==TRUE)
            form[pa[1],common.missing]=0
            form[pa[2],common.missing]=0

		obs.match <- length(which(form[pa[1],]==form[pa[2],]))                                       

	      prob.cor <- cbind(irtprob(param=ip,theta=theta.est1),irtprob(param=ip,theta=theta.est2)) #probability correct
			
			colnames(prob.cor) <- c(paste("Examinee ",pa[1],sep=""),
							paste("Examinee ",pa[2],sep=""))	
			row.names <- c("Item 1");for(i in 2:ncol(form)){row.names <- c(row.names,paste("Item ",i,sep="")) }
			rownames(prob.cor) <- row.names

		prob.incor <- 1-prob.cor #probability incorrect

			colnames(prob.incor) <- c(paste("Examinee ",pa[1],sep=""),
							paste("Examinee ",pa[2],sep=""))	
			rownames(prob.incor) <- row.names


		Pi <- (prob.cor[,1]*prob.cor[,2])+(prob.incor[,1]*prob.incor[,2]) #probability of matching

  			Qi <-1-Pi
			   Cpi <- cumprod(Pi)
			     I <- length(Pi)
	
      			M <- matrix(1,(I + 1),I)
		      	M[1,] <- cumprod(Qi)
			       for(o in 1:I) {
				   M[(o+1),o]<- Cpi[o]
			       }

	      	for(m in 2:(I+1)) {
				for(o in 2:I) {
					if(m <= o)
						M[m,o] <- Qi[o]* M[m,(o-1)]+Pi[o]* M[(m-1),(o-1)]
					   else M[m,o] <- M[m,o]
				}
			}

	GBT.p.value <- sum(M[(obs.match+1):(I+1),I])

		matchings <- c("Prob. of 0 Match","Prob. of 1 Match");for(i in 2:ncol(form)){matchings<- c(matchings,paste("Prob. of ",i," Matches",sep="")) }
		prob.dist.match <- as.data.frame(cbind(matchings,round(M[,I],6)))

		return(list(ability1=theta.est1,
				ability2=theta.est2,
				prob.cor = prob.cor,
				prob.match=Pi,
				exact.prob.dist=prob.dist.match,
				p.value=GBT.p.value
			))
	
	}#end internal function


	GBT.index <- GBT(form=data,ip=item.par,pa=pair)


#####################                    Computing K-index                 #####################################################

	k <- function(form,pa) { #start internal function


			wc  <- ncol(form)-rowSums(form[pa[1],],na.rm=TRUE)                    #number-incorrect score for copier
			ws  <- ncol(form)-rowSums(form[pa[2],],na.rm=TRUE)                    #number-incorrect score for source

			
	            m   <- length(which(form[pa[1],]==0 & form[pa[2],]==0))+
                         length(which(is.na(form[pa[1],])==TRUE & is.na(form[pa[2],])==TRUE))  #number of identical incorrect responses between two vectors

			subgroup <- which((ncol(form)-rowSums(form,na.rm=TRUE))==wc)

				if(length(subgroup)!=0) {

				smatrix <- as.data.frame(matrix(rep(as.matrix(form[pa[2],]),length(subgroup)),nrow=length(subgroup),byrow=TRUE))

				emp.agg <- rowSums((form[subgroup,]==0)&(smatrix==0),na.rm=TRUE)+
                                   rowSums((is.na(form[subgroup,])==TRUE & is.na(smatrix)==TRUE),na.rm=TRUE)
				p=mean(emp.agg)/ws 

				 } else p=NA

		if(is.na(p)!=TRUE) { k.index=1-pbinom(m-1,ws,p) } else k.index=NA

		return(list(
				subgroups=subgroup,
				emp.agg=emp.agg,
				k.index=k.index
			))

	
	}#end internal function

	if(rowSums(data[pair[2],],na.rm=TRUE)!=ncol(data)) { k.index <- k(form=data,pa=pair) } else k.index <- NULL


##############                       Computing K variants                   #######################################

	ks12 <- function(form,pa) { #start internal function

		options(warn=-1)
		subgroups <- vector("list",ncol(form)+1)

			for(j in 1:(ncol(form)+1)){
				subgroups[[j]] <- which((ncol(form)-rowSums(form,na.rm=TRUE))==j-1)
			}

	
			wc  <- ncol(form)-rowSums(form[pa[1],],na.rm=TRUE)
			qc  <- wc/ncol(form)
			ws  <- ncol(form)-rowSums(form[pa[2],],na.rm=TRUE)
	            m   <- length(which(form[pa[1],]==0 & form[pa[2],]==0))+
				 length(which(is.na(form[pa[1],])==TRUE & is.na(form[pa[2],])==TRUE))
			cm <- which(form[pa[1],]==1 & form[pa[2],]==1)

			pr     <- c()
			prob   <- matrix(nrow=(ncol(form)+1),ncol=ncol(form))
			weight <- matrix(nrow=(ncol(form)+1),ncol=ncol(form))
			pj <- c()

			g=.2
			d2=-(1+g)/g

				for(j in 1:(ncol(form)+1)){

					if(length(subgroups[[j]])!=0) {

						smatrix <- as.data.frame(matrix(rep(as.matrix(form[pa[2],]),length(subgroups[[j]])),nrow=length(subgroups[[j]]),byrow=TRUE))
						emp.agg <- rowSums((form[subgroups[[j]],]==0)&(smatrix==0),na.rm=TRUE)+
                                               rowSums((is.na(form[subgroups[[j]],])==TRUE & is.na(smatrix)==TRUE),na.rm=TRUE)
						pr[j]=mean(emp.agg)/ws 
						cor.match <- ((form[subgroups[[j]],]==1 & smatrix==1)*1)
						cor.match[which(is.na(cor.match)==TRUE)]=0
						prob[j,] <- colMeans(cor.match)
						weight[j,] <- (((1+g)/(1-g))*exp(1))^(prob[j,]*d2)
						pj[j] <-  mean(cor.match%*%t(t(weight[j,])))

					 } else 
				
					if(length(subgroups[[j]])==0) {
							pr[j]=NA 
							pj[j]=NA
					}
				}


			Qrs  <- (0:ncol(form))/ncol(form)
			Qrs2 <- Qrs^2
			Qrs3 <- 0:ncol(form)
			mm <- ceiling(sum(weight[wc+1,cm]))+m
				
				pred1 <- predict(lm(pr~1+Qrs))
				pred.1 <- c()
				for(i in 1:(ncol(form)+1)){ 
					if(length(which(as.character(1:(ncol(form)+1))[i]==names(pred1)))!=0){
						pred.1[i]=pred1[which(as.character(1:(ncol(form)+1))[i]==names(pred1))]} else pred.1[i]=NA
                        }
				pred2 <- predict(lm(pr~1+Qrs+Qrs2))
				pred.2 <- c()
				for(i in 1:(ncol(form)+1)){ 
					if(length(which(as.character(1:(ncol(form)+1))[i]==names(pred2)))!=0){
						pred.2[i]=pred2[which(as.character(1:(ncol(form)+1))[i]==names(pred2))]} else pred.2[i]=NA
                        }

				pred3 <- exp(predict(glm(ws*pr ~ Qrs3 ,family=poisson())))
				pred.3 <- c()
				for(i in 1:(ncol(form)+1)){ 
					if(length(which(as.character(1:(ncol(form)+1))[i]==names(pred3)))!=0){
						pred.3[i]=pred3[which(as.character(1:(ncol(form)+1))[i]==names(pred3))]} else pred.3[i]=NA
                        }
				pred4 <- exp(predict(glm(ws*pr+ceiling(pj) ~ Qrs3 ,family=poisson())))
				pred.4 <- c()
				for(i in 1:(ncol(form)+1)){ 
					if(length(which(as.character(1:(ncol(form)+1))[i]==names(pred4)))!=0){
						pred.4[i]=pred4[which(as.character(1:(ncol(form)+1))[i]==names(pred4))]} else pred.4[i]=NA
                        }


				p1 <- pred.1[which(Qrs==qc)]
				p2 <- pred.2[which(Qrs==qc)]
				s1 <- pred.3[which(Qrs==qc)]
				s2 <- pred.4[which(Qrs==qc)]


			if(is.na(p1)!=TRUE & p1>=1)  { p1=.999 };if(is.na(p1)!=TRUE & p1<=0) { p1=.001 }
			if(is.na(p2)!=TRUE & p2>=1)  { p2=.999 };if(is.na(p2)!=TRUE & p2<=0) { p2=.001 }
			if(is.na(s1)!=TRUE & s1>=ws) { s1=ws }  ;if(is.na(s2)!=TRUE & s2>=ncol(form)) { s2=ncol(form) }  ;

			if(is.na(p1)!=TRUE) { k1.index=1-pbinom(m-1,ws,p1) } else k1.index=NA
			if(is.na(p2)!=TRUE) { k2.index=1-pbinom(m-1,ws,p2) } else k2.index=NA
			if(is.na(s1)!=TRUE) { s1.index=(1-ppois(m-1,s1)) - (1 - ppois(ws,s1))       } else s1.index=NA
			if(is.na(s2)!=TRUE) { s2.index=(1-ppois(mm-1,s2)) - (1 - ppois(ncol(form),s2))      } else s2.index=NA


		return(list(mean.iden.incorrect = pr*ws,
				weighted.iden.correct=pj,
				subgroups=subgroups,
				pred1 = pred.1*ws,
				pred2 = pred.2*ws,
				pred3 = pred.3,
				pred4 = pred.4,
				K1.index =k1.index,
				K2.index =k2.index,
				S1.index =s1.index,
				S2.index =s2.index))

	}#end internal function

	if(rowSums(data[pair[2],],na.rm=TRUE)!=ncol(data)) { k.variants <- ks12(form=data,pa=pair) } else k.variants <- NULL

	outCD <- list(data=data,suspected.pair=pair,W.index=w.index,GBT.index=GBT.index,
                    K.index=k.index,K.variants=k.variants)
	class(outCD) <- "CopyDetect"
	return(outCD)


}#end CopyDetect1 function

















