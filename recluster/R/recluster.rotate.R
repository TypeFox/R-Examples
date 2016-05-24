recluster.rotate<-function(table,m=FALSE,q=FALSE,flip="none",centre=TRUE){	
    table2<-table
    if(m&q){
	for (i in 1:nrow(table))
		{
		dist<-(abs(-m*table[i,1]+table[i,2]-q))/sqrt(m^2+1)
		if (table[i,2] < (m*table[i,1]+q)){dist<--dist}
		table2[i,2]<-dist
		m2<--1/m
		q2<-table[i,2]-m2*table[i,1]
		x1<-(q2-q)/(m-m2)
		y1<-m*x1+q
		table2[i,1]<-sqrt((0-x1)^2+(q-y1)^2)*sign(x1)
		}
		}	
	if(flip=="hor"){table2[,1]<--table2[,1]}
	if(flip=="ver"){table2[,2]<--table2[,2]}
	if(flip=="both"){table2[,1]<--table2[,1]
		table2[,2]<--table2[,2]}
	if(centre){table2[,1]<-table2[,1]-mean(table2[,1])
	table2[,2]<-table2[,2]-mean(table2[,2])}
    return (table2)
}