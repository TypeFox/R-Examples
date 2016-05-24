bos.cat2 <-
function(data,nodevec,col){
	#nodevec: indicator variable that takes the value 1 if a subject falls into that node and 0 if not
	#col= categorical predictor variable to be used for the split, is not in the dataset
	#f2=measure of effect size in multiple regression (see Cohen,1988)
	x<-col[nodevec==1]
	##Transformation according to mean on residuals (see Algorithm 2 of paper)
	y<-resid(lm(data))	
	y<-y[nodevec==1]
	meanx<-data.frame(sort(tapply(y,x,mean)))
	ordernames<-row.names(meanx)
	trx2<-factor(col,levels=ordernames)
    trx2<-as.numeric(trx2)
	object<-list(ordernames=ordernames,trx=trx2)
   }
