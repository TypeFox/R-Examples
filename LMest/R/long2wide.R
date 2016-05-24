long2wide <-
function(data,nameid,namet,colx,coly,aggr=T,full=999){
	
# preliminaries	
	data = as.matrix(data[,c(nameid,namet,colx,coly)])
	listid = data[,nameid]
	listid = sort(unique(listid))
	listt = data[,namet]
	listt = sort(unique(listt))
	n = length(listid)
	TT = length(listt)
	nx = length(colx)
	ny = length(coly)
# wide data
	data_wide = matrix(full,n,TT*(nx+ny))
	ind2 = c(colx,coly)
	for(j in 1:dim(data)[1]){
		i = which(data[j,nameid]==listid)
		t = which(data[j,namet]==listt) 
		ind = (t-1)*(nx+ny)+(1:(nx+ny))
		data_wide[i,ind] = data[j,3:(2+nx+ny)]
	}
	data_wide = cbind(listid,data_wide)
# put names
	namev = colnames(data[,-2])
	if(is.null(namev)){
		namev = "id"
		for(j in 1:nx) namev = c(namev,paste("x",1,sep=""))
		for(j in 1:ny) namev = c(namev,paste("y",1,sep=""))
	}
	namev_wide = namev[1]
	for(t in 1:TT) namev_wide = c(namev_wide,paste(namev[2:(nx+ny+1)],"_",listt[t],sep=""))
	colnames(data_wide) = namev_wide
# if aggregate data are required
	if(aggr){
		if(is.na(full)) data_wide[is.na(data_wide)]=999
		out = aggr_data(data_wide[,2:(1+TT*(nx+ny))])
		S = out$data_dis
		S = array(t(S),c(nx+ny,TT,length(out$freq)))
		S = aperm(S)
		if(is.na(full)) S[S==999]=NA
		XX = S[,,1:nx]
		YY = S[,,(nx+1):(nx+ny)]
		freq = out$freq
	}
	else {XX=NULL;YY=NULL;freq=NULL}
# final output
	out = list(listid=listid,listt=listt,data_wide=data_wide,XX=XX,YY=YY,freq=freq)
}
