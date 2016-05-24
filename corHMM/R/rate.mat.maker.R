#Rate matrix maker and manipulating functions

#written by Jeremy M. Beaulieu and Jeffrey C. Oliver

rate.mat.maker<-function(rate.cat, hrm=TRUE, ntraits=NULL, nstates=NULL, model=c("ER", "SYM", "ARD")){
	if(hrm==TRUE){
		k=2
		mat1 <- matrix(NA, k*rate.cat, k*rate.cat)
		mat2 <- matrix(NA, k*rate.cat, k*rate.cat)
		vec.tmp1<-rep(c(0,1),rate.cat)
		vec.tmp2<-rep(1:rate.cat, rep(2,rate.cat))-1
		
		for(i in 1:(k*rate.cat)){
			mat1[i,]<-abs(vec.tmp1-vec.tmp1[i])
			mat2[i,]<-abs(vec.tmp2-vec.tmp2[i])
		}
		matFINAL<-mat1+mat2
		rate.mat.index <- matrix(NA, k*rate.cat, k*rate.cat)
		np <- k + (rate.cat-1) * 6
		index<-matFINAL==1
		rate.mat.index[index] <- 1:np
		if (rate.cat == 1){
			rownames(rate.mat.index) <- c("(0)","(1)")
			colnames(rate.mat.index) <- c("(0)","(1)")			
		}
		if (rate.cat == 2){
			rownames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
			colnames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
		}
		if (rate.cat == 3){
			rownames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
			colnames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
		}
		if (rate.cat == 4){
			rownames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
			colnames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
		}
		if (rate.cat == 5){
			rownames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
			colnames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
		}
		if (rate.cat == 6){
			rownames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)")
			colnames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)")
		}
		if (rate.cat == 7){
			rownames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)","(0,R7)","(1,R7)")
			colnames(rate.mat.index) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)","(0,R6)","(1,R6)","(0,R7)","(1,R7)")
		}		
	}
	if(hrm==FALSE){
		k=ntraits
		nl=2
		if(ntraits==1){
			k <- 1
			nl <- nstates
			if (is.character(model)) {
				rate.mat.index <- matrix(NA, nl, nl) #An nl x nl matrix, filled with NA values
				tmp2 <- cbind(1:(nl^k), 1:(nl^k)) # For setting diagonals
				index<-matrix(TRUE,nl^k,nl^k)
				diag(index) <- FALSE
				if (model == "ER") {
					np <- 1 #np is the number of parameters in the rate matrix
					rate.mat.index[index] <- 1:np
				}
				if (model == "SYM") {
					np <- nl * (nl - 1)/2
					sel <- col(rate.mat.index) < row(rate.mat.index)
					#Use transpose of the rate category matrix to finish enumerating:
					rate.mat.index <- t(rate.mat.index)
					rate.mat.index[sel] <- 1:np
				}
				if (model == "ARD") {
					np <- nl*(nl - 1)
					rate.mat.index[index] <- 1:np
				}
			}
		}
		if(ntraits==2){
			#Hard-coded for now
			mat1<-matrix(,nl^k,nl^k)
			mat2<-matrix(,nl^k,nl^k)
			vec.tmp1<-c(0,0,1,1)
			vec.tmp2<-c(0,1,0,1)
			for(i in 1:(nl^k)){
				mat1[i,]<-abs(vec.tmp1-vec.tmp1[i])
				mat2[i,]<-abs(vec.tmp2-vec.tmp2[i])
			}
			matFINAL<-mat1+mat2
			
			if (is.character(model)) {
				rate.mat.index <- matrix(NA, nl^k, nl^k)
				if (model == "ER"){
					np <- 1
					index<-matFINAL==1
					rate.mat.index[index] <- 1:np
				}
				if (model == "SYM") {
					np <- 4
					index<-matFINAL==1
					rate.mat.index[index][c(1,2,4,6)] <- rate.mat.index[index][c(3,5,7,8)] <- 1:np
				}
				if (model == "ARD") {
					np <- 8
					index<-matFINAL==1
					rate.mat.index[index] <- 1:np
				}
			}
		}
		if(ntraits==3){
			#Hard-coded for now
			mat1<-matrix(,nl^k,nl^k)
			mat2<-matrix(,nl^k,nl^k)
			mat3<-matrix(,nl^k,nl^k)
			vec.tmp1<-c(0,1,0,0,1,1,0,1)
			vec.tmp2<-c(0,0,1,0,1,0,1,1)
			vec.tmp3<-c(0,0,0,1,0,1,1,1)
			
			for(i in 1:(nl^k)){
				mat1[i,]<-abs(vec.tmp1-vec.tmp1[i])
				mat2[i,]<-abs(vec.tmp2-vec.tmp2[i])
				mat3[i,]<-abs(vec.tmp3-vec.tmp3[i])
			}
			matFINAL<-mat1+mat2+mat3
			
			if (is.character(model)) {
				rate.mat.index <- matrix(NA, nl^k, nl^k)
				if (model == "ER"){
					np <- 1
					index<-matFINAL==1
					rate.mat.index[index] <- 1:np
				}
				if (model == "SYM") {
					np <- 12
					index<-matFINAL==1
					rate.mat.index[index][c(1,2,3,5,6,8,9,11,12,15,18,21)] <- rate.mat.index[index][c(4,7,10,13,16,14,19,17,20,22,23,24)] <- 1:np
				}
				if (model == "ARD") {
					np <- 24
					index<-matFINAL==1
					rate.mat.index[index] <- 1:np
				}
			}
		}
		if(ntraits==1){
			rownames(rate.mat.index) <- as.character(1:nl)
			colnames(rate.mat.index) <- as.character(1:nl)
		}
		if(ntraits==2){
			rownames(rate.mat.index) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
			colnames(rate.mat.index) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
		}
		if(ntraits==3){
			rownames(rate.mat.index) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
			colnames(rate.mat.index) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
		}			
	}

	return(rate.mat.index)
}

rate.par.drop <- function(rate.mat.index=NULL,drop.par=NULL){
	if(is.null(rate.mat.index)){
		cat("Rate matrix needed.  See mat.maker to create one.\n")
		return
	}
	if(is.null(drop.par)){
		cat("No parameters indicated to drop.  Original matrix returned.\n")
		return(rate.mat.index)
	}
	if(max(rate.mat.index,na.rm=TRUE) < max(drop.par,na.rm=TRUE)){
		cat("Some parameters selected for dropping were not in the original matrix.\n")
	}
	drop.par <- unique(drop.par) # in case parameters listed more than once in drop vector
	drop.par <- drop.par[order(drop.par)]
	max <- max(rate.mat.index,na.rm=TRUE)
	for(drop.which in 1:length(drop.par)){
		drop.locs <- which(rate.mat.index == drop.par[drop.which],arr.ind=TRUE)
		rate.mat.index[drop.locs] <- NA
	}
	max <- max - length(drop.par)
	exclude <- which(is.na(rate.mat.index))
	rate.mat.index[-exclude] <- 1:max
	
	return(rate.mat.index)
}

rate.par.eq <- function(rate.mat.index=NULL,eq.par=NULL){
	if(is.null(rate.mat.index)){
		cat("Rate matrix needed.  See mat.maker to create one.\n")
		return
	}
	if(is.null(drop) || length(eq.par) < 2){
		cat("Fewer than two parameters indicated to equalize.  Original matrix returned.\n")
		return(rate.mat.index)
	}
	too.big <- which(eq.par > max(rate.mat.index,na.rm=TRUE))
	if(length(too.big) > 0){
		cat("Some parameters selected for equalizing were not in the original matrix:\n")
		cat("Not in original rate.mat.index:",eq.par[too.big],"\n")
		cat("Original matrix returned.\n")
		return(rate.mat.index)
	}
	eq.par <- unique(eq.par)
	eq.par <- eq.par[order(eq.par)]	
	min <- min(eq.par) # rm.na unnecessary?

	# the decrement index will hold counters to decrement rate index
	dec.index <- matrix(0,length(rate.mat.index[,1]),length(rate.mat.index[1,]))
	for(eq.which in 2:length(eq.par)){
		to.eq <- which(rate.mat.index == eq.par[eq.which],arr.ind=TRUE)
		rate.mat.index[to.eq] <- min
	}
	# the decrement index will hold counters to decrement rate index
	dec.index <- matrix(0,length(rate.mat.index[,1]),length(rate.mat.index[1,]))
	for(eq.which in 2:length(eq.par)){
		to.dec <- which(rate.mat.index > eq.par[eq.which],arr.ind=TRUE) #greater than current decrementer
		dec.index[to.dec] <- dec.index[to.dec] + 1
	}
	rate.mat.index <- rate.mat.index - dec.index
		
	return(rate.mat.index)
}

