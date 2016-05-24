pgdata <-
function(data,plot=NULL,type=NULL,...){

	if(is.null(plot))plot <- FALSE
	if(is.null(type))type <- c("barplot")
	res <- pgdata1(data,plot,type)
	return(res)

}



pgdata1 <-
function(data,plot=NULL,type=NULL,...){
  	n <- ncol(data)
  	rn <- nrow(data)
  	level <- NULL
  	for(i in 1:n) level <- union(data[,i],level)
  	id.na = which(level == "NA")
	if(length(id.na)==0) id.na = which(is.na(level))
  	if(length(id.na)==1)level <- union(level[-id.na],level[id.na])
  	num <- length(level)
  	ID <- matrix(0,n,num)
  	for(i in 1:n){
		miss <- 0
		id <- matrix(0,2,num)
		for(j in 1:num){
			if(!(is.na(level[j]))){ 
				id[1,j] <- length(which(data[,i]==level[j]))
				id[2,j] <- id[1,j]/rn
			}
			if(is.na(level[j])){
				miss <- j
			}
		}	
      	if (miss != 0){
			id[1,miss] <- rn - sum(id[1,])
			id[2,miss] <- id[1,miss]/rn
		}
		ID[i,] <- id[2,]
  	}
  	colnames(ID) <- level
	if(is.null(plot))plot <- FALSE
	if(is.null(type))type <- c("barplot")
      if(plot){
	if(type=="barplot")barplot(t(ID[,num]),col="red",ylab="proportion of missings")
	if(type=="stacked"){
		layout(matrix(c(1,2), nrow = 1), widths = c(0.95, 0.05))
		par(mar = c(5,3,5,1)+0.1)
		n = dim(ID)[2]
		barplot(t(ID),col=rainbow(n),legend.text=rownames(t(ID)),args.legend=list(x=3,y=-0.01))
	}
	if(type=="dist"){	
		missing.id <- missing.dist(data)	
		x = missing.id[,1]
		y = missing.id[,2]
		plot(x = x,y = y, ylim = rev(range(y)),pch=20,col='lightblue',main="Distribution of Missings", xlab="", ylab="")
	}
      }
	return(ID)
}
