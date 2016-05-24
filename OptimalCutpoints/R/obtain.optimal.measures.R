obtain.optimal.measures <-
function(value, measures.acc) {
	position <- which(measures.acc$cutoffs %in% value)
	
	Se.v <- measures.acc$Se[position,,drop = FALSE]
	Sp.v <- measures.acc$Sp[position,,drop = FALSE]	
	PPV.v <- measures.acc$PPV[position,,drop = FALSE]	
	NPV.v <- measures.acc$NPV[position,,drop = FALSE]	
	DLR.Positive.v <- measures.acc$DLR.Positive[position,,drop = FALSE]
	DLR.Negative.v <- measures.acc$DLR.Negative[position,,drop = FALSE]

	FP <- measures.acc$n$h*(1-Sp.v)
	FN <- measures.acc$n$d*(1-Se.v)
	if(ncol(Se.v) == 3) { #Confidence intervals
		FP[,c(2,3)] <- NA
		FN[,c(2,3)] <- NA
	}
	res <- list(cutoff = value, Se = Se.v, Sp = Sp.v, PPV = PPV.v, NPV = NPV.v, DLR.Positive = DLR.Positive.v, DLR.Negative = DLR.Negative.v, FP = FP, FN = FN)
	return(res)
}
