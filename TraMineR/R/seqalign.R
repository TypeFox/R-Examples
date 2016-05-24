seqalign <- function(seqdata, indices, indel=1, sm, with.missing = FALSE) {
	
	l1 <- c(seqlength(seqdata[indices[1],]))
	#print(str(l1))
	l2 <- c(seqlength(seqdata[indices[2],]))
	#seqdata <- seqnum(seqdata, with.missing=with.missing)
	seq1 <- seqdata[indices[1],]
	seq2 <- seqdata[indices[2],]
	allalphabet <- alphabet(seqdata)
	if(with.missing) {
		allalphabet <- c(allalphabet, attr(seqdata, "nr"))
	}
	rownames(sm) <- allalphabet
	colnames(sm) <- allalphabet
	#print(seq1)
	#print(seq2)
	leven <- matrix(0, nrow=(l1+1), ncol=(l2+1))
	## Init first indels row and col
	leven[,1] <- 0:l1*indel
	leven[1,] <- 0:l2*indel
	operation <- matrix("", nrow=(l1+1), ncol=(l2+1))
	operation[, 1] <- "D"
	operation[1,] <- "I"
	## Filling the matrix
    for(i in 2:(l1+1)) {
        for(j in 2:(l2+1)) {
            if (seq1[1,i-1] == seq2[1,j-1]) {
                cost <- 0
				## Comparison equal
            } else {
                cost <-  sm[as.character(seq1[1,i-1]),as.character(seq2[1,j-1])];
				##Substitution
            }
            leven[i,j] <- min(leven[i,j-1]   + indel,
							  leven[i-1,j]   + indel,
							  leven[i-1,j-1] + cost)
			op <- which.min(c(leven[i-1,j-1] + cost,
							  leven[i,j-1]   + indel,
							  leven[i-1,j]   + indel))
			operation[i,j] <- c("S","I", "D")[op]
			if(op==1 && cost==0) {
				operation[i,j] <- "E"
			}
        }
    }
	i <- l1+1
	j <- l2+1
	#print(operation)
	#print(leven)
	operationlist <- character(length=l1+l2)
	seq1c <- character(length=l1+l2)
	seq2c <- character(length=l1+l2)
	cost <- numeric(length=l1+l2)
	opi <- 1
	while(TRUE) {
		op <- operation[i,j]
		oldcost <- leven[i,j]
		# dest_i <- i - 1
		# dest_j <- j - 1
		# if (op=="I") {
			# dest_i <- i
		# } else if(op=="D") {
			# dest_j <- j
		# }
		operationlist[opi] <- op
		
		if(op=="S" || op=="E") {
			seq1c[opi] <- as.character(seq1[1,i-1])
			seq2c[opi] <- as.character(seq2[1,j-1])
			i <- i-1
			j <- j-1
		} else if (op=="I") {
			seq1c[opi] <- "-"
			seq2c[opi] <- as.character(seq2[1,j-1])
			j <- j-1
		} else if (op=="D") {
			seq2c[opi] <- "-"
			seq1c[opi] <- as.character(seq1[1,i-1])
			i <- i-1
		}
		cost[opi] <- oldcost-leven[i,j]
		opi <- opi+1
	#	print(c(i,j))
		#print(operation[i,j])
		if(i==1&&j==1){
			break
		}
	}
		
	
	
	revopi <- function(l){
		return(rev(l[1:(opi-1)]))
	}
	
	ret <- list()
	ret$operation <- revopi(operationlist)
	ret$seq1 <- revopi(seq1c)
	ret$seq2 <- revopi(seq2c)
	ret$cost <- revopi(cost)
	ret$opmatrix <- operation
	ret$costmatrix <- leven
	ret$stsseq <- seqdata[indices,]
	class(ret) <- "seqalign"
	return(ret)
}
print.seqalign <- function(x, digits=3, ...){
	formd <- function(n){
		return(format(n, digits=digits))
	}
	toprint <- matrix("", nrow=4, ncol=length(x$operation))
	toprint[1,] <- x$seq1
	toprint[2,] <- x$operation
	toprint[3,] <- formd(x$cost)
	toprint[4,] <- x$seq2
	rownames(toprint) <- c("Seq 1", "Operations", "Costs", "Seq 2")
	colnames(toprint) <- as.character(1:length(x$operation))
	print(toprint, quote=FALSE, ...)
	omdist <- sum(x$cost)
	cat("\nOM distance between sequences:", formd(omdist), "\n")
}
  #data(biofam)
  #biofam.seq <- seqdef(biofam, 10:25)
 # costs <- seqsubm(biofam.seq, method="TRATE")
  #costs <- seqsubm(biofam.seq, method="CONSTANT")
 # sa <- seqalign(biofam.seq, 1:2, indel=1, sm=costs)
  #sa <- seqalign(biofam.seq, c(1,5), indel=0.5, sm=costs)