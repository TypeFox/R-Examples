# calculate entropy
# p can be count vector or probability vector, but not a vector of membership indicator
H=function (p, logbase=c("e","2")) { 
    logbase <- match.arg(logbase)    
    if (sum(p)!=1) p=p/sum(p) # if p is count, transform to probability
    p=p[p!=0] # remove zero entries
    if (logbase=="e") {
        sum(-p*log(p)) 
    } else if (logbase=="2") {
        sum(-p*log2(p)) 
    }
}

# mutual information
mutual.info=function(two.way.table, logbase=c("e","2")){
    H(rowSums(two.way.table), logbase=logbase) + H(colSums(two.way.table), logbase=logbase) - H(c(two.way.table), logbase=logbase)    
}
# test
#mutual.info(matrix(c(1,1,1,1),2,2))

cor.mixed <- function(x, ...) UseMethod("cor.mixed") 

cor.mixed.default=function(x, na.fun, method=c("pearson","spearman"), ...) {
    p=ncol(x)
    res=matrix(1,p,p, dimnames=list(colnames(x), colnames(x)))
    for (i in 2:p) {
        for (j in 1:(i-1)){
            res[i,j]<-res[j,i]<-cor.mixed.vector (x[,i], x[,j], na.fun, method) 
        }
    }
    res
}

cor.mixed.vector=function(x, y, na.fun, method=c("pearson","spearman"), ...) {
    method=match.arg(method)
    x=ifelse(na.fun(x), NA, x)
    y=ifelse(na.fun(y), NA, y)
    mean(is.na(x) & is.na(y)) + mean(!is.na(x) & !is.na(y)) * cor(x,y,method=method,use="p")
}

cor.mixed.formula=function(formula, data, na.fun, method=c("pearson","spearman"), ...) {
    vars=dimnames(attr(terms(formula),"factors"))[[1]]
    x=data[,vars[1]]; 
    y=data[,vars[2]]; 
    cor.mixed.vector(x,y,na.fun,method)
}


# information coefficient of correlation
info.cor = function(two.way.table) {
    I=mutual.info(two.way.table)
    sqrt(1-exp(-2*I))
}
## test: 
#info.cor(matrix(c(1,0,0,1),2,2))

yule.y=function(two.by.two.matrix) {
    (sqrt(two.by.two.matrix[1,1] * two.by.two.matrix[2,2]) - sqrt(two.by.two.matrix[1,2] * two.by.two.matrix[2,1])) / 
    (sqrt(two.by.two.matrix[1,1] * two.by.two.matrix[2,2]) + sqrt(two.by.two.matrix[1,2] * two.by.two.matrix[2,1])) 
}

kappa.cor=function(two.by.two.matrix, weight=c(1,1), maximum=FALSE) {
    total=sum(two.by.two.matrix)
    
#    E=matrix(c(
#        sum(two.by.two.matrix[1,])*sum(two.by.two.matrix[,1])/total,     
#        sum(two.by.two.matrix[1,])*sum(two.by.two.matrix[,2])/total,     
#        sum(two.by.two.matrix[2,])*sum(two.by.two.matrix[,1])/total,     
#        sum(two.by.two.matrix[2,])*sum(two.by.two.matrix[,2])/total
#    ), 2, 2) 
#    w=matrix(c(0,weight,0),2,2)
#    1-sum(two.by.two.matrix*w)/sum(E*w)
    
    f1=sum(two.by.two.matrix[1,])
    f2=sum(two.by.two.matrix[2,])
    g1=sum(two.by.two.matrix[,1])
    g2=sum(two.by.two.matrix[,2])
    p.e=(f1*g1+f2*g2)/total^2    
    p.0=(two.by.two.matrix[1,1]+two.by.two.matrix[2,2])/total
    if (maximum) {
        max.p.0=ifelse (which.min(c(f1,f2,g1,g2)) %in% c(1,4), (f1+g2)/total, (f2+g1)/total)
        (p.0-p.e)/(max.p.0-p.e)
    } else {
        (p.0-p.e)/(1-p.e)
    }
}

# L-measure, only works for two.by.two.matrix b/c of the way I.sup is computed
l.measure = function(two.by.two.matrix) {
    I=mutual.info(two.by.two.matrix)
    I.sup = mutual.info(matrix(c(1,0,0,1),2,2))
    sqrt(1-exp(-2*I / (1-I/I.sup) ))
}
## test: 
# l.measure(matrix(c(1,0,0,1),2,2))
