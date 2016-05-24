aggr_data <- function(data,disp=FALSE,fort=FALSE){

    # find distinct response configurations and corresponding frequencies
    # does not work properly with missing data
    data = as.matrix(data)
    data_dis = data
    rim = nrow(data)    
    nc = ncol(data)
    if(fort){
    	out = .Fortran("aggrdata",as.double(data),as.integer(rim),as.integer(nc),ndis=as.integer(0),
               datadis=matrix(0,rim,nc),freq=as.integer(rep(0,rim)),
               label=as.integer(rep(0,rim)))
		data_dis = out$datadis[1:out$ndis,]; freq = out$freq[1:out$ndis]; label = out$label    	
    }else{
	    freq = rep(0,rim)
	    label = rep(0,rim)
	    j = 0
	    label0 = 1:rim
	    if(disp){
	        cat("------------|-------------|\n")
	   	    cat("  remaining |   average   |\n")
	   	    cat("------------|-------------|\n")
	   	}
	    while(rim>0){
	    	j = j+1
	    	data_dis[j,] = data[1,]
			D = t(matrix(data_dis[j,],nc,rim))
		    ind = which(rowSums(D==data)==nc)
		    label[label0[ind]] = j
			freq[j] = length(ind)    
			data = as.matrix(data[-ind,])
			label0 = label0[-ind]
		    rim = length(label0)
		    if(rim==1) data = t(data)
		    if(disp) cat(sprintf("%11g",c(rim,mean(freq[1:j]))),"\n",sep=" | ")    
	    } 
	   	if(disp) cat("------------|-------------|\n");
	    data_dis = as.matrix(data_dis[1:j,])
	    freq = freq[1:j]
	}
    # output
    out = list(data_dis=data_dis,freq=freq,label=label)
}