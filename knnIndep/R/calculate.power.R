calculate.power <-
function(vals,alpha=.95,comp=`>`){
	sapply(vals,function(v){
				## find quant% significance level cutoff
				cut = apply(sapply(v,function(lis){lis$null}),1,quantile,probs=alpha,na.rm=T)
				if(is.null(dim(cut))) dim(cut) = c(1,length(cut))
				## Now we estimate the power as the number of alternative statistics exceeding our estimated cutoffs
				vals = sapply(v,function(lis){lis$dep})
				exceeding = sapply(1:nrow(cut),function(r){sapply(1:nrow(vals),function(i){comp(vals[i,],cut[r,i])})},simplify="array")
				#return(colSums(exceeding)/length(v))
				return(apply(exceeding,c(2,3),sum,na.rm=T)/length(v))
			},simplify="array")
}
