CopyDetect2 <- function(data,item.par=NULL,pair,options,key=NULL) { #Start Function

#########################################################################################################

	# data,  N x n response data
		# "data" should be a data frame object, the variables should be all character variables
		# N is number of examinee and n is number of items
		# Nominal responses 

	# item.par, n x 10 item parameter matrix, in the logistic metric
		#first r columns are discrimination parameters, where r is the number of options in the test
		#second r columns are intercept parameters
		#the key is inferred from item parameter matrix
		
	# pair, a vector of length two. 

		#the first element of the vector is the row number of the suspected copier examinee in 
		 #in the data file
		#the second element of the vector is the row number of the source examinee for the 
		 #suspected copier examinee in the data file.


	# One of the incorrect options are randomly assigned to missing data. 
	
	#CopyDetect computes w, GBT, S2, S1, K, K1, and K2 indices internally, and print the information
		#regarding the likelihood of agreement between suspected copier and source examinees'
		#response vector.

#########################################################################################################

for(i in 1:ncol(data)) {
				
			if(is.character(data[,i])!=TRUE) {

				stop("Variable ",i," is not a character variable. All variables (columns) in the data file
has to be character variable. Use as.character() to transform your variables")
			}
		}


for(i in 1:ncol(data)){

	if(sum(options==sort(unique(data[,1]))[which(sort(unique(data[,1]))!="NA")])!=length(options)){
	stop("Variable ",i, " includes a character element other than what has been specified as an input for
possible options")
	}
}


r <- ncol(item.par)/2

	if(r!=length(options)){
                            stop("The number of columns in the item parameter matrix does not match with the 
number of options in the data file. Check your input files.")
                            cat("        Item parameter file should have", length(options),"* 2 columns","\n")
                           }

#Key

if(is.null(item.par)==TRUE & is.null(key)==TRUE){
stop("Either key responses or item parameters must be provided.")
}



if(is.null(item.par)!=TRUE){
	key2 <- c()
	for(i in 1:ncol(data)) { key2[i]=options[which.max(item.par[i,1:r])]}
}

if(is.null(item.par)!=TRUE){
	
	if(sum(key==key2)!=length(key)){

		 stop("Check your input for key. It does not match with the item parameter matrix.
In Nominal Response Model, the response alternative with the highest slope parameter is
the key response.")
	}
}

if(is.null(key)==TRUE){ key <- key2}

			

			
	row.names <- c("Examinee1");for(i in 2:nrow(data)){row.names <- c(row.names,paste("Examinee",i,sep="")) }
	rownames(data) <- row.names



#Check item parameter file is appropriate

	
		if(ncol(item.par)!=length(options)*2)          {
                                            stop("The item parameter file should have ",length(options),"*2 columns Check your input matrix")
                                           }
		if(nrow(item.par)!=ncol(data)) {
                                            stop("The number of rows in item parameter file should be equal to the number of 
columns in the response data file. Check your input matrices")
							 }
		if(is.data.frame(data)!=TRUE)  {
							  stop("The input response data file is not a data frame object. Use as.data.frame() 
first to make your response data file a data frame object")
							 }

		if(length(pair)!=2)		 {
							  stop("The input pair should be a vector of length two. The first element indicates 
the row number in the response data file for the suspected copier examinee, and the second element indicates the row number in 
the response data file for the suspected source examinee")
							 }

		

		if(pair[1]>nrow(data) | pair[2] > nrow(data)) { stop("The elements of the input vector pair can not be bigger than 
the number of rows in the response data file")
									    }


#

scored.data <- as.data.frame(matrix(nrow=nrow(data),ncol=ncol(data)))
for(i in 1:nrow(scored.data)){ scored.data[i,] <- as.numeric((data[i,]==key)*1)}


		
#Internal function to compute probability matrix for an ability level of "x"

		irtprob <- function(ability,item.param) {  
		
			prob <- matrix(nrow=nrow(item.param),ncol=ncol(item.param)/2)

			for(i in 1:nrow(prob)){
				ps <- c()
				for(j in 1:ncol(prob)){ps[j]=exp((item.param[i,j]*ability)+item.param[i,j+ncol(prob)])
							    }
				prob[i,]=ps/sum(ps)
                   }
		prob
		}

#Using scored data, estimate ability levels needed for GBT and w index

ipar.dic <- est(scored.data, model = "2PL", engine = "ltm")$est
abilities <- mlebme(scored.data,ipar.dic)[,1]


######################             Computing w index       #########################################################

	omega <- function(form,thetas,ip,pa,options) { #start internal function

		key <- c()
		for(i in 1:ncol(data)) { key[i]=options[which.max(ip[i,1:r])]}

		theta.est1 <- thetas[pa[1]]

		obs.match <- length(which(form[pa[1],]==form[pa[2],]))                                       

	      probabilities <- irtprob(ability=theta.est1,item.param=ip)

			colnames(probabilities) <- options
			row.names <- c("Item 1");for(i in 2:ncol(form)){row.names <- c(row.names,paste("Item ",i,sep="")) }
			rownames(probabilities) <- row.names

			miss.items <- which(form[pa[2],]=="NA")
                          if(length(miss.items)==0) { miss.items <- which(is.na(form[pa[2], ]) == TRUE)}

			for(i in miss.items){ form[pa[2], i] = options[which(options != key[i])][which.max(probabilities[i,which(options != key[i])])] }

				pvec <- c()
				for(i in 1:ncol(form)){ pvec[i]=probabilities[i,which(options==form[pa[2],i])] }
                                        
		exp.match <- sum(pvec)
		sd.match  <- sqrt(sum(pvec*(1-pvec)))

		w.value <- (obs.match-exp.match)/sd.match 
		p.value <- pnorm(w.value,0,1,lower.tail=FALSE)
		
		return(list(exp.match=exp.match,obs.match=obs.match,
				sd.match=sd.match,
				W.value=w.value,
			 	p.value=p.value))
	}#end internal function

	if(is.null(item.par)!=TRUE) {
			w.index <- omega(form=data,thetas=abilities,ip=item.par,pa=pair,options=options)
		} else w.index <- NULL


##################                     Computing GBT index              ##################################################


	GBT <- function(form,thetas,ip,pa) { #start internal function
	
		theta.est1 <- thetas[pa[1]]  
		theta.est2 <- thetas[pa[2]]

		obs.match <- length(which(form[pa[1],]==form[pa[2],]))                                       

	      probabilities1 <- irtprob(ability=theta.est1,item.param=ip)

			colnames(probabilities1) <- options
			row.names <- c("Item 1");for(i in 2:ncol(form)){row.names <- c(row.names,paste("Item ",i,sep="")) }
			rownames(probabilities1) <- row.names

		probabilities2 <- irtprob(ability=theta.est2,item.param=ip)

			colnames(probabilities2) <- options
			row.names <- c("Item 1");for(i in 2:ncol(form)){row.names <- c(row.names,paste("Item ",i,sep="")) }
			rownames(probabilities2) <- row.names

		Pi <- c()
		for(i in 1:ncol(form)){ Pi[i]=sum(probabilities1[i,]*probabilities2[i,])}

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

		return(list(probabilities1 = probabilities1,probabilities2 = probabilities2,
				prob.match=Pi,
				exact.prob.dist=prob.dist.match,
				p.value=GBT.p.value
			))
	
	}#end internal function

	if(is.null(item.par)!=TRUE) {

		GBT.index <- GBT(form=data,thetas=abilities,ip=item.par,pa=pair)
	} else GBT.index <- NULL 


#####################                    Computing K-index                 #####################################################

	k <- function(form,form2,pa) { #start internal function


			wc  <- ncol(form2)-rowSums(form2[pa[1],])                    #number-incorrect score for copier
			ws  <- ncol(form2)-rowSums(form2[pa[2],])                    #number-incorrect score for source

			incorrect.items <- which(form2[pa[2],]==0)
			m <- length(which(form[pa[1],incorrect.items]==form[pa[2],incorrect.items]))
			
			subgroup <- which((ncol(form2)-rowSums(form2))==wc)

				if(length(subgroup)!=0) {

				incorrect.items <- which(form2[pa[2],]==0)
				smatrix <- as.data.frame(matrix(rep(as.matrix(form[pa[2],incorrect.items]),length(subgroup)),nrow=length(subgroup),byrow=TRUE))
				emp.agg <- rowSums(form[subgroup,incorrect.items]==smatrix)
				p=mean(emp.agg)/ws 

				 } else p=NA

		if(is.na(p)!=TRUE) { k.index=1-pbinom(m-1,ws,p) } else k.index=NA

		return(list(
				subgroups=subgroup,
				emp.agg=emp.agg,
				k.index=k.index
			))

	
	}#end internal function

	if(rowSums(scored.data[pair[2],])!=ncol(data)) { k.index <- k(form=data,form2=scored.data,pa=pair) } else k.index <- NULL


##############                       Computing K variants                   #######################################

	ks12 <- function(form,form2,pa) { #start internal function

		options(warn=-1)
		subgroups <- vector("list",ncol(form)+1)

			for(j in 1:(ncol(form)+1)){
				subgroups[[j]] <- which((ncol(form2)-rowSums(form2))==j-1)
			}

	
			wc  <- ncol(form2)-rowSums(form2[pa[1],])
			qc  <- wc/ncol(form)
			ws  <- ncol(form2)-rowSums(form2[pa[2],])
			incorrect.items <- which(form2[pa[2],]==0)
			m <- length(which(form[pa[1],incorrect.items]==form[pa[2],incorrect.items]))
                  cm <- which(form2[pa[1],]==1 & form2[pa[2],]==1)

			pr     <- c()
			prob   <- matrix(nrow=(ncol(form)+1),ncol=ncol(form))
			weight <- matrix(nrow=(ncol(form)+1),ncol=ncol(form))
			pj <- c()

			g=1/length(options)
			d2=-(1+g)/g

				for(j in 1:(ncol(form)+1)){

					if(length(subgroups[[j]])!=0) {

						incorrect.items <- which(form2[pa[2],]==0)
						smatrix1 <- as.data.frame(matrix(rep(as.matrix(form[pa[2],incorrect.items]),length(subgroups[[j]])),nrow=length(subgroups[[j]]),byrow=TRUE))
						smatrix2 <- as.data.frame(matrix(rep(as.matrix(form2[pa[2],]),length(subgroups[[j]])),nrow=length(subgroups[[j]]),byrow=TRUE))
						emp.agg <- rowSums(form[subgroups[[j]],incorrect.items]==smatrix1)
						pr[j]=mean(emp.agg)/ws 
						prob[j,] <- colMeans((form2[subgroups[[j]],]==1)&(smatrix2==1))
						weight[j,] <- (((1+g)/(1-g))*exp(1))^(prob[j,]*d2)
						pj[j] <-  mean(((form2[subgroups[[j]],]==1 & smatrix2==1)*1)%*%t(t(weight[j,])))

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
			if(is.na(s1)!=TRUE & s1>=ws) { s1=ws }  ;if(is.na(s1)!=TRUE & s2>=ncol(form)) { s2=ncol(form) }  ;

			if(is.na(p1)!=TRUE) { k1.index=1-pbinom(m-1,ws,p1) } else k1.index=NA
			if(is.na(p2)!=TRUE) { k2.index=1-pbinom(m-1,ws,p2) } else k2.index=NA
			if(is.na(s1)!=TRUE) { s1.index= (1-ppois(m-1,s1)) - (1 - ppois(ws,s1))       } else s1.index=NA
			if(is.na(s2)!=TRUE) { s2.index= (1-ppois(mm-1,s2)) - (1 - ppois(ncol(form),s2))     } else s2.index=NA


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

	if(rowSums(scored.data[pair[2],])!=ncol(data)) { k.variants <- ks12(form=data,form2=scored.data,pa=pair) } else k.variants <- NULL

	outCD <- list(data=data, key=key, 
                    scored.data=scored.data,
                    theta.par=abilities,
                    suspected.pair=pair,
                    W.index=w.index,GBT.index=GBT.index,
                    K.index=k.index,K.variants=k.variants)
	class(outCD) <- "CopyDetect2"
	return(outCD)


}#end CopyDetect2 function

















