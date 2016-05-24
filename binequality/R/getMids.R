getMids <-
function(ID, hb, lb, ub, alpha_bound = 10/9){
	ID <- as.character(ID)
	mids.out <- c()
	c.out <- c()
	alpha.out <- c()
	ID.out <- c()
	hb.out <- c()
	for(id in unique(ID)){
		use.id <- which(ID == id & hb >0)
  		open.lb <- which(is.na(lb[use.id])==TRUE)
  		open.ub <- which(is.na(ub[use.id])==TRUE)
  		if(length(open.lb)==0&length(open.ub)==0){
    			mids <- (ub[use.id] + lb[use.id])/2
    			alpha <- NA
    			c <- NA
  		}#end if length(open.lb)==0&length(open.ub)==0
  		if(length(open.lb)>0){
    			stop('The code is not written to handle left censored data',' current ID is ', id,'\n','\n')
  		}
  		if(length(open.ub)>2){
    			stop('The code is not written to handle more than 1 right censored bin',' current ID is ', id,'\n','\n')
  		}
  		if(length(open.ub)==1){
    			mids<-rep(NA, length(ub[use.id]))
    			mids[-open.ub]<-(ub[use.id][-open.ub] + lb[use.id][-open.ub])/2
    			use<-which(hb[use.id]>0)
    			hb.use<-hb[use.id][use]
    			lb.use<-lb[use.id][use]
    			ub.use<-ub[use.id][use]
    			open.ub.use <- which(is.na(ub.use)==TRUE)
    
    			a.num<-log((hb.use[(open.ub.use-1)]+hb.use[open.ub.use])/hb.use[open.ub.use])
    			a.denom<-log(lb.use[open.ub.use]/lb.use[(open.ub.use-1)])
    			alpha<-a.num/a.denom
    
    			if(length(alpha_bound) == 0){
    				warning('alpha is unbounded', '\n')
    			}else{
    				alpha<-max(alpha, alpha_bound)
    			}
    			c = alpha/(alpha - 1)
    			mids[open.ub]<-lb[use.id][open.ub]*c
  			}#end if length(open.ub) == 1 
  		mids.out <- c(mids.out, mids)
		c.out <- c(c.out, c)
		alpha.out <- c(alpha.out, alpha)
		ID.out <- c(ID.out, ID[use.id])
		hb.out <- c(hb.out, hb[use.id])
  		}#end for id
	mids.return <- data.frame(ID.out,mids.out,hb.out)
	colnames(mids.return) <- c('ID', 'mids', 'hb')
	return(list('mids' = mids.return, 'c' = c.out, 'alpha' = alpha.out))
}