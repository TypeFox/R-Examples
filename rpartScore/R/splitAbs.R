splitAbs<-function (y, wt, x, parms, continuous,minb){
    	n <- length(y)
    	no.classes <- length(table(y))
    	class.labs <- as.numeric(names(table((y))))
    	node<-xtabs(wt~y)/sum(wt)
	cum.node<-cumsum(node)
	root <- 0
    	for (j in 1:(no.classes - 1)) {
        	root <- root + 2*(class.labs[j+1]-class.labs[j])*cum.node[j]*(1-cum.node[j])
    	}

    	if (continuous) {

		if (no.classes==1) list(goodness = rep(0,n-1), direction = rep(-1,n-1))
		else {
        	ymat <- matrix(rep(0, length(x) * no.classes), ncol = no.classes)
        	for (j in 1:no.classes) {
            	ymat[which(y == class.labs[j]),j] <- wt[which(y==class.labs[j])]
        	}
        	left <- apply(ymat, 2, cumsum)
        	col.left <- apply(left, 1, sum)
        	row.total <- apply(ymat, 2, sum)
        	impure.l <- left/col.left
        	right <- matrix(rep(row.total, length(x)), ncol = no.classes, byrow = TRUE) - left
        	col.right <- apply(right, 1, sum)
        	impure.r <- right/col.right
        	cum.matrix.l <- t(apply(impure.l, 1, cumsum))
        	cum.matrix.r <- t(apply(impure.r, 1, cumsum))
        	left.imp <- 0
        	right.imp <- 0
        	for (j in 1:(no.classes - 1)) {
            	left.imp <- left.imp + 2*(class.labs[j+1]-class.labs[j])*cum.matrix.l[,j]*(1-cum.matrix.l[,j])
            	right.imp <- right.imp + 2*(class.labs[j+1]-class.labs[j])*cum.matrix.r[,j]*(1-cum.matrix.r[, j])
       	}
        	impure <- (root - col.left/sum(wt) * left.imp - col.right/sum(wt)*right.imp)
        	impure <- impure[-n]
        	direction <- rep(-1, (n - 1))
        	list(goodness = impure, direction = direction)
		}
    	} else {
		
		num.cat<-length(unique(x))
		if (no.classes==1) list(goodness = rep(0,num.cat-1), direction = unique(x))
		else {
		left<-expand.grid(rep(list(c(0,1)),(num.cat)))
		left<-left[1:(nrow(left)/2),]
		left<-left[-1,]


		
		colnames(left)<-unique(x)
		
		tabella<-xtabs(wt~y+x)
		nodo.l<-function(a) apply(as.table(tabella[,a==1]),1,sum)
	
		freq.l<-apply(left,1,nodo.l)
		freq.r<-as.numeric(xtabs(wt~y))-freq.l
		col.left<-apply(freq.l,2,sum)
		col.right<-apply(freq.r,2,sum)
		freq.cum.l<-apply(freq.l,2,cumsum)
		freq.cum.r<-apply(freq.r,2,cumsum)
		cum.matrix.l<- t(freq.cum.l)/col.left
		cum.matrix.r<- t(freq.cum.r)/col.right
		

		left.imp <- 0
        	right.imp <- 0
        	for (j in 1:(no.classes - 1)) {
            	left.imp <- left.imp + 2*(class.labs[j+1]-class.labs[j])*cum.matrix.l[,j]*(1-cum.matrix.l[,j])
            	right.imp <- right.imp + 2*(class.labs[j+1]-class.labs[j])*cum.matrix.r[,j]*(1-cum.matrix.r[, j])
       	}
        	impure <- (root - col.left/sum(wt) * left.imp - col.right/sum(wt)*right.imp)
		
		check.left<-col.left<minb
		check.right<-col.right<minb
		check.minbucket<-(check.left+check.right)==0
		
		if (sum(check.minbucket)>0){
			left1<-left[check.minbucket,]
			impure1<-impure[check.minbucket]
			split.ott<-which.max(impure1)
			pos.cat.left<-left1[split.ott,]
			pos.cat.right<-1-pos.cat.left
			cat.left<-unique(x)[pos.cat.left==1]
			cat.right<-unique(x)[pos.cat.right==1]
			direction<-c(cat.left,cat.right)
			goodness<-rep(max(impure1)*0.1,num.cat-1)
			goodness[sum(pos.cat.left)]<-max(impure1)
		} else {
			direction<-unique(x)
			goodness<-rep(0,num.cat-1)
		}
		list(goodness = goodness, direction = direction)
    	}
	}
}
