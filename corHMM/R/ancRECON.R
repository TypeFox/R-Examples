#RECONSTRUCTION OF ANCESTRAL STATES

#written by Jeremy M. Beaulieu and Jeffrey C. Oliver

ancRECON <- function(phy, data, p, method=c("joint", "marginal", "scaled"), hrm=FALSE, rate.cat, ntraits=NULL, charnum=NULL, rate.mat=NULL, model=c("ER", "SYM", "ARD"), root.p=NULL){
	
	#Note: Does not like zero branches at the tips. Here I extend these branches by just a bit:
	phy$edge.length[phy$edge.length<=1e-5]=1e-5
	
	if(hrm==FALSE){
		if(ntraits==1){
			data.sort<-data.frame(data[,charnum+1],data[,charnum+1],row.names=data[,1])
		}
		if(ntraits==2){
			data.sort<-data.frame(data[,2], data[,3], row.names=data[,1])
		}
		if(ntraits==3){
			data.sort<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
		}
	}
	else{
		data.sort <- data.frame(data[,2], data[,2],row.names=data[,1])
	}
	data.sort<-data.sort[phy$tip.label,]
	
	#Some initial values for use later
	k=2
	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	#Builds the rate matrix based on the specified rate.cat. Not exactly the best way
	#to go about this, but it is the best I can do for now -- it works, so what me worry?
	if(hrm==TRUE){
		if(is.null(rate.mat)){
			rate<-rate.mat.maker(hrm=TRUE,rate.cat=rate.cat)
			rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
		}
		else{
			rate<-rate.mat
			rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
		}
		#Makes a matrix of tip states and empty cells corresponding 
		#to ancestral nodes during the optimization process.	
		x <- data.sort[,1]
		TIPS <- 1:nb.tip
		
		for(i in 1:nb.tip){
			if(is.na(x[i])){x[i]=2}
		}
		if (rate.cat == 1){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			TIPS <- 1:nb.tip
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,1]=1}
				if(x[i]==1){liks[i,2]=1}
				if(x[i]==2){liks[i,1:2]=1}
			}
		}
		if (rate.cat == 2){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3)]=1}
				if(x[i]==1){liks[i,c(2,4)]=1}
				if(x[i]==2){liks[i,1:4]=1}
			}
		}
		if (rate.cat == 3){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5)]=1}
				if(x[i]==1){liks[i,c(2,4,6)]=1}
				if(x[i]==2){liks[i,1:6]=1}
			}
		}
		if (rate.cat == 4){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5,7)]=1}
				if(x[i]==1){liks[i,c(2,4,6,8)]=1}
				if(x[i]==2){liks[i,1:8]=1}
			}
		}
		if (rate.cat == 5){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5,7,9)]=1}
				if(x[i]==1){liks[i,c(2,4,6,8,10)]=1}
				if(x[i]==2){liks[i,1:10]=1}
			}
		}
		if (rate.cat == 6){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5,7,9,11)]=1}
				if(x[i]==1){liks[i,c(2,4,6,8,10,12)]=1}
				if(x[i]==2){liks[i,1:12]=1}
			}
		}
		if (rate.cat == 7){
			liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
			for(i in 1:nb.tip){
				if(x[i]==0){liks[i,c(1,3,5,7,9,11,13)]=1}
				if(x[i]==1){liks[i,c(2,4,6,8,10,12,14)]=1}
				if(x[i]==2){liks[i,1:14]=1}
			}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
		tranQ <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if(hrm==FALSE){
		#Imported from Jeffs rayDISC -- will clean up later, but for now, it works fine:
		if(ntraits==1){
			k <- 1
			factored <- factorData(data.sort) # was acting on data, not data.sort
			nl <- ncol(factored)
			obj <- NULL
			nb.tip<-length(phy$tip.label)
			nb.node <- phy$Nnode
			
			if(is.null(rate.mat)){
				rate<-rate.mat.maker(hrm=FALSE,ntraits=ntraits,nstates=nl,model=model)
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}
			else{
				rate<-rate.mat
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}
			
			stateTable <- NULL # will hold 0s and 1s for likelihoods of each state at tip
			for(column in 1:nl){
				stateTable <- cbind(stateTable,factored[,column])
			}
			colnames(stateTable) <- colnames(factored)
			
			ancestral <- matrix(0,nb.node,nl) # all likelihoods at ancestral nodes will be 0
			liks <- rbind(stateTable,ancestral) # combine tip likelihoods & ancestral likelihoods
			rownames(liks) <- NULL
		}
		if(ntraits==2){
			k=2
			nl=2
			if(is.null(rate.mat)){
				rate<-rate.mat.maker(hrm=FALSE,ntraits=ntraits,model=model)
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}
			else{
				rate<-rate.mat
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}
			x<-data.sort[,1]
			y<-data.sort[,2]
			
			liks <- matrix(0, nb.tip + nb.node, nl^k)
			TIPS <- 1:nb.tip
			for(i in 1:nb.tip){
				if(is.na(x[i])){
					x[i]=2
				}
				if(is.na(y[i])){
					y[i]=2
				}
			}
			for(i in 1:nb.tip){
				if(x[i]==0 & y[i]==0){liks[i,1]=1}
				if(x[i]==0 & y[i]==1){liks[i,2]=1}
				if(x[i]==1 & y[i]==0){liks[i,3]=1}
				if(x[i]==1 & y[i]==1){liks[i,4]=1}
				if(x[i]==2 & y[i]==0){liks[i,c(1,3)]=1}
				if(x[i]==2 & y[i]==1){liks[i,c(2,4)]=1}
				if(x[i]==0 & y[i]==2){liks[i,c(1,2)]=1}
				if(x[i]==1 & y[i]==2){liks[i,c(3,4)]=1}
				if(x[i]==2 & y[i]==2){liks[i,1:4]=1}
			}
		}
		if(ntraits==3){
			k=3
			nl=2
			if(is.null(rate.mat)){
				rate<-rate.mat.maker(hrm=FALSE,ntraits=ntraits,model=model)
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}
			else{
				rate<-rate.mat
				rate[is.na(rate)]<-max(rate,na.rm=TRUE)+1
			}			
			x<-data.sort[,1]
			y<-data.sort[,2]
			z<-data.sort[,3]
			
			liks <- matrix(0, nb.tip + nb.node, nl^k)
			TIPS <- 1:nb.tip
			for(i in 1:nb.tip){
				if(is.na(x[i])){
					x[i]=2
				}
				if(is.na(y[i])){
					y[i]=2
				}
				if(is.na(z[i])){
					z[i]=2
				}
			}
			for(i in 1:nb.tip){
				if(x[i]==0 & y[i]==0 & z[i]==0){liks[i,1]=1}
				if(x[i]==1 & y[i]==0 & z[i]==0){liks[i,2]=1}
				if(x[i]==0 & y[i]==1 & z[i]==0){liks[i,3]=1}
				if(x[i]==0 & y[i]==0 & z[i]==1){liks[i,4]=1}
				if(x[i]==1 & y[i]==1 & z[i]==0){liks[i,5]=1}
				if(x[i]==1 & y[i]==0 & z[i]==1){liks[i,6]=1}
				if(x[i]==0 & y[i]==1 & z[i]==1){liks[i,7]=1}
				if(x[i]==1 & y[i]==1 & z[i]==1){liks[i,8]=1}
				#If x is ambiguous but the rest are not:
				if(x[i]==2 & y[i]==0 & z[i]==0){liks[i,c(1,2)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==0){liks[i,c(3,5)]=1}
				if(x[i]==2 & y[i]==0 & z[i]==1){liks[i,c(4,6)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==1){liks[i,c(7,8)]=1}
				#If y is ambiguous but the rest are not:
				if(x[i]==0 & y[i]==2 & z[i]==0){liks[i,c(1,3)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==0){liks[i,c(2,5)]=1}
				if(x[i]==0 & y[i]==2 & z[i]==1){liks[i,c(4,7)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==1){liks[i,c(6,8)]=1}
				#If z is ambiguous but the rest are not:
				if(x[i]==0 & y[i]==0 & z[i]==2){liks[i,c(1,4)]=1}
				if(x[i]==0 & y[i]==1 & z[i]==2){liks[i,c(3,7)]=1}
				if(x[i]==1 & y[i]==0 & z[i]==2){liks[i,c(2,6)]=1}
				if(x[i]==1 & y[i]==1 & z[i]==2){liks[i,c(5,8)]=1}
				#If x and y is ambiguous but z is not:
				if(x[i]==2 & y[i]==2 & z[i]==0){liks[i,c(1,2,3,5)]=1}
				if(x[i]==2 & y[i]==2 & z[i]==1){liks[i,c(4,6,7,8)]=1}
				#If x and z is ambiguous but y is not:
				if(x[i]==2 & y[i]==0 & z[i]==2){liks[i,c(1,2,4,6)]=1}
				if(x[i]==2 & y[i]==1 & z[i]==2){liks[i,c(3,5,7,8)]=1}
				#If y and z is ambiguous but x is not:
				if(x[i]==0 & y[i]==2 & z[i]==2){liks[i,c(1,3,4,7)]=1}
				if(x[i]==1 & y[i]==2 & z[i]==2){liks[i,c(2,5,6,8)]=1}
				#All states are ambiguous:			
				if(x[i]==2 & y[i]==2 & z[i]==2){liks[i,1:8]=1}
			}
		}
		Q <- matrix(0, nl^k, nl^k)
		tranQ <- matrix(0, nl^k, nl^k)
	}
	
    p[p==0] = exp(-21)
	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)
	phy <- reorder(phy, "pruningwise")
	TIPS <- 1:nb.tip
	anc <- unique(phy$edge[,1])
	
	if(method=="joint"){
		lik.states<-numeric(nb.tip + nb.node)
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		for (i  in seq(from = 1, length.out = nb.node)) {
			#The ancestral node at row i is called focal:
			focal <- anc[i]
			#Get descendant information of focal:
			desRows<-which(phy$edge[,1]==focal)
			#Get node information for each descendant:
			desNodes<-phy$edge[desRows,2]
			#Initiates a loop to check if any nodes are tips:
			for (desIndex in sequence(length(desRows))){
				#If a tip calculate C_y(i) for the tips and stores in liks matrix:
				if(any(desNodes[desIndex]==phy$edge[,1])==FALSE){
					if(hrm==TRUE){
						v <- c(rep(1, k*rate.cat))
					}
					if(hrm==FALSE){
						v <- c(rep(1, nl^k))
					}
					Pij <- expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77"))
					v = v * liks[desNodes[desIndex],]
					for(i in 1:dim(Pij)[1]){
						L <- Pij[i,] * v
						liks[desNodes[desIndex],i] <- max(L)
						comp[desNodes[desIndex],i] <- which.max(L==max(L))[1]
					}					
				}
			}
			#Collects t_z, or the branch subtending focal:
			tz<-phy$edge.length[which(phy$edge[,2] == focal)]	
			if(length(tz)==0){
				#The focal node is the root, calculate P_k:
				root.state=1
				for (desIndex in sequence(length(desRows))){
					#This is the basic marginal calculation:
					root.state <- root.state * liks[desNodes[desIndex],]
				}
                equil.root <- NULL
                for(i in 1:ncol(Q)){
                    posrows <- which(Q[,i] >= 0)
                    rowsum <- sum(Q[posrows,i])
                    poscols <- which(Q[i,] >= 0)
                    colsum <- sum(Q[i,poscols])
                    equil.root <- c(equil.root,rowsum/(rowsum+colsum))
                }
                if (is.null(root.p)){
                    flat.root = equil.root
                    k.rates <- 1/length(which(!is.na(equil.root)))
                    flat.root[!is.na(flat.root)] = k.rates
                    flat.root[is.na(flat.root)] = 0
                    liks[focal, ] <- root.state
                }
                else{
                    if(is.character(root.p)){
                        # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10), q01/(q01+q10), etc.
                        if(root.p == "yang"){
                            Q.tmp <- Q
                            diag(Q.tmp) = 0
                            root.p = colSums(Q.tmp) / sum(Q.tmp)
                            liks[focal, ] <- root.state * root.p
                        }else{
                            # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
                            root.p = root.state / sum(root.state)
                            liks[focal,] <- root.state * root.p
                        }
                    }
                    # root.p!==NULL will fix root probabilities based on user supplied vector:
                    else{
                        liks[focal, ] <- root.state * root.p
                    }
                }

			}
			#All other internal nodes, except the root:
			else{
				#Calculates P_ij(t_z):
				Pij <- expm(Q * tz, method=c("Ward77"))
				#Calculates L_z(i):
				if(hrm==TRUE){
					v<-c(rep(1, k*rate.cat))
				}
				if(hrm==FALSE){
					v<-c(rep(1, nl^k))
				}
				for (desIndex in sequence(length(desRows))){
					v = v * liks[desNodes[desIndex],]
				}
				for(i in 1:dim(Pij)[1]){
					L <- Pij[i,] * v
					liks[focal,i] <- max(L)
					comp[focal,i] <- which.max(L==max(L))[1]
				}
				sum.tot <- sum(liks[focal,])
				liks[focal,]=liks[focal,]/sum.tot
				
#				if(sum(liks[focal,])<1e-200){
					#Kicks in arbitrary precision calculations: 
#					library(Rmpfr)
#					liks <- mpfr(liks, 15)	
#				}
			}
		}
		root <- nb.tip + 1L
        lik.states[root] <- which(liks[root,]==max(liks[root,]))
		N <- dim(phy$edge)[1]
		for(i in N:1){
			anc <- phy$edge[i,1]
			des <- phy$edge[i,2]
			lik.states[des] <- comp[des,lik.states[anc]]
		}
		#For later use:
		#logl <- as.numeric(log(liks[root,lik.states[root]]))
		#Outputs likeliest tip states
		obj$lik.tip.states <- lik.states[TIPS]
		#Outputs likeliest node states
		obj$lik.anc.states <- lik.states[-TIPS]
	}
	if(method=="marginal"){
		#A temporary likelihood matrix so that the original does not get written over:
		liks.down<-liks
		#A transpose of Q for assessing probability of j to i, rather than i to j:
		tranQ<-t(Q)
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		#The first down-pass: The same algorithm as in the main function to calculate the conditional likelihood at each node:
		for (i  in seq(from = 1, length.out = nb.node)) {
			#the ancestral node at row i is called focal
			focal <- anc[i]
			#Get descendant information of focal
			desRows<-which(phy$edge[,1]==focal)
			desNodes<-phy$edge[desRows,2]
			v <- 1
			for (desIndex in sequence(length(desRows))){
				v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks.down[desNodes[desIndex],]
			}
			comp[focal] <- sum(v)
			liks.down[focal, ] <- v/comp[focal]
		}
		root <- nb.tip + 1L
		#Enter the root defined root probabilities if they are supplied by the user:
        equil.root <- NULL
        for(i in 1:ncol(Q)){
            posrows <- which(Q[,i] >= 0)
            rowsum <- sum(Q[posrows,i])
            poscols <- which(Q[i,] >= 0)
            colsum <- sum(Q[i,poscols])
            equil.root <- c(equil.root,rowsum/(rowsum+colsum))
        }
        if (is.null(root.p)){
            flat.root = equil.root
            k.rates <- 1/length(which(!is.na(equil.root)))
            flat.root[!is.na(flat.root)] = k.rates
            flat.root[is.na(flat.root)] = 0
            liks.down[root, ] <- flat.root * liks.down[root, ]
            liks.down[root, ] <- liks.down[root,] / sum(liks.down[root, ])
            root.p = flat.root
        }
        else{
            if(is.character(root.p)){
                # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10), q01/(q01+q10), etc.
                if(root.p == "yang"){
                    Q.tmp <- Q
                    diag(Q.tmp) = 0
                    root.p = colSums(Q.tmp) / sum(Q.tmp)
                    liks.down[root, ] <- root.p * liks.down[root, ]
                    liks.down[root, ] <- liks.down[root,] / sum(liks.down[root,])
                }else{
                    # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
                    root.p = liks.down[root,] / sum(liks.down[root,])
                    liks.down[root, ] <- root.p * liks.down[root, ]
                    liks.down[root, ] <- liks.down[root,] / sum(liks.down[root, ])
                }
            }
            # root.p!==NULL will fix root probabilities based on user supplied vector:
            else{
                liks.down[root, ] <- root.p * liks.down[root, ]
                liks.down[root, ] <- liks.down[root,] / sum(liks.down[root, ])
            }
        }
		#The up-pass
		liks.up<-liks
		states<-apply(liks,1,which.max)
		N <- dim(phy$edge)[1]
		comp <- numeric(nb.tip + nb.node)
		for(i in length(anc):1){
			focal <- anc[i]
			if(!focal==root){
				#Gets mother and sister information of focal:
				focalRow<-which(phy$edge[,2]==focal)
				motherRow<-which(phy$edge[,1]==phy$edge[focalRow,1])
				motherNode<-phy$edge[focalRow,1]
				desNodes<-phy$edge[motherRow,2]
				sisterNodes<-desNodes[(which(!desNodes==focal))]
				sisterRows<-which(phy$edge[,2]%in%sisterNodes==TRUE)
				#If the mother is not the root then you are calculating the probability of being in either state.
				#But note we are assessing the reverse transition, j to i, rather than i to j, so we transpose Q to carry out this calculation:
				if(motherNode!=root){
					v <- expm(tranQ * phy$edge.length[which(phy$edge[,2]==motherNode)], method=c("Ward77")) %*% liks.up[motherNode,]
				}
				#If the mother is the root then just use the marginal. This can also be the prior, which I think is the equilibrium frequency. 
				#But for now we are just going to use the marginal at the root -- it is unclear what Mesquite does.
				else{
					v <- root.p
				}
				#Now calculate the probability that each sister is in either state. Sister can be more than 1 when the node is a polytomy. 
				#This is essentially calculating the product of the mothers probability and the sisters probability:
				for (sisterIndex in sequence(length(sisterRows))){
					v <- v*expm(Q * phy$edge.length[sisterRows[sisterIndex]], method=c("Ward77")) %*% liks.down[sisterNodes[sisterIndex],]
				}
				comp[focal] <- sum(v)
				liks.up[focal,] <- v/comp[focal]
			}
		}
		#The final pass
		liks.final<-liks
		comp <- numeric(nb.tip + nb.node)
		#In this final pass, root is never encountered. But its OK, because root likelihoods are set after the loop:
		for (i in seq(from = 1, length.out = nb.node-1)) { 
			#the ancestral node at row i is called focal
			focal <- anc[i]
			focalRows<-which(phy$edge[,2]==focal)
			#Now you are assessing the change along the branch subtending the focal by multiplying the probability of 
			#everything at and above focal by the probability of the mother and all the sisters given time t:
			v <- liks.down[focal,]*expm(tranQ * phy$edge.length[focalRows], method=c("Ward77")) %*% liks.up[focal,]
			comp[focal] <- sum(v)
			liks.final[focal, ] <- v/comp[focal]
		}
		#Just add in the marginal at the root calculated on the original downpass or if supplied by the user:
		liks.final[root,] <- liks.down[root,] * root.p
		root.final <- liks.down[root,] * root.p
		comproot <- sum(root.final)
		liks.final[root,] <- root.final/comproot
		#Reports the probabilities for all internal nodes as well as tips:
		#Outputs likeliest tip states
		obj$lik.tip.states <- liks.final[TIPS,]
		#Outputs likeliest node states
		obj$lik.anc.states <- liks.final[-TIPS,]
	}	
	
	if(method=="scaled"){
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		#The same algorithm as in the main function. See comments in either corHMM.R, corDISC.R, or rayDISC.R for details:
		for (i  in seq(from = 1, length.out = nb.node)) {
			#the ancestral node at row i is called focal
			focal <- anc[i]
			#Get descendant information of focal
			desRows<-which(phy$edge[,1]==focal)
			desNodes<-phy$edge[desRows,2]
			v <- 1
			for (desIndex in sequence(length(desRows))){
				v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
			}
			comp[focal] <- sum(v)
			liks[focal, ] <- v/comp[focal]
		}
        equil.root <- NULL
        for(i in 1:ncol(Q)){
            posrows <- which(Q[,i] >= 0)
            rowsum <- sum(Q[posrows,i])
            poscols <- which(Q[i,] >= 0)
            colsum <- sum(Q[i,poscols])
            equil.root <- c(equil.root,rowsum/(rowsum+colsum))
        }
        if (is.null(root.p)){
            flat.root = equil.root
            k.rates <- 1/length(which(!is.na(equil.root)))
            flat.root[!is.na(flat.root)] = k.rates
            flat.root[is.na(flat.root)] = 0
            liks[root, ] <- flat.root * liks[root, ]
            liks[root, ] <- liks[root,] / sum(liks[root, ])
        }
        else{
            if(is.character(root.p)){
                # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10), q01/(q01+q10), etc.
                if(root.p == "yang"){
                    Q.tmp <- Q
                    diag(Q.tmp) = 0
                    root.p = colSums(Q.tmp) / sum(Q.tmp)
                    liks[root, ] <- root.p * liks[root, ]
                    liks[root, ] <- liks[root,] / sum(liks[root, ])
                }else{
                    # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:
                    root.p = liks[root,] / sum(liks[root,])
                    liks[root, ] <- root.p * liks[root, ]
                    liks[root, ] <- liks[root,] / sum(liks[root, ])
                }
            }
            # root.p!==NULL will fix root probabilities based on user supplied vector:
            else{
                liks[root, ] <- root.p * liks[root, ]
                liks[root, ] <- liks[root,] / sum(liks[root, ])
            }
        }
		#Reports the probabilities for all internal nodes as well as tips:
		obj$lik.tip.states <- liks[TIPS,]
		#Outputs likeliest node states
		obj$lik.anc.states <- liks[-TIPS,]
	}	
	obj
}

