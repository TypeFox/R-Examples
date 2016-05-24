Aggregate <-
function(x,by) {


categoricals <- data.frame(x[sapply(x, is.factor)],x[sapply(x, is.character)])
dummies <- dummy.data.frame(data=categoricals,dummy.classes="ALL", sep='_')

dummies_sum <- aggregate(dummies,list(ID=by),sum)
dummies_tail <- aggregate(dummies,list(ID=by),tail,n=1)
ID <- dummies_tail$ID
dummies_sum$ID <- NULL
dummies_tail$ID <- NULL

names(dummies_sum) <- paste(names(dummies_sum),'_sum',sep='')
names(dummies_tail) <- paste(names(dummies_tail),'_last',sep='')
dummies <- data.frame(ID,dummies_sum,dummies_tail)

numerics <- x[sapply(x, is.numeric )]

numerics_sum <- aggregate(numerics,by=list(ID=by),sum)
numerics_mean <- aggregate(numerics,by=list(ID=by),mean)
numerics_var <- aggregate(numerics,by=list(ID=by),var)
ID <- numerics_sum$ID
numerics_sum$ID <- NULL
numerics_mean$ID <- NULL
numerics_var$ID <- NULL
names(numerics_sum) <- paste(names(numerics_sum),'_sum',sep='')
names(numerics_mean) <- paste(names(numerics_mean),'_mean',sep='')
names(numerics_var) <- paste(names(numerics_var),'_var',sep='')
numerics <- data.frame(numerics_sum,numerics_mean,numerics_var)

final <- data.frame(dummies,numerics)
final
}
