recluster.line<-function(mat,type="maxd",X1=NULL,X2=NULL){
res<-NULL
if (type=="maxd"){
    maxd<-max(dist(mat))
    matr<-as.matrix(dist(mat))
    res$first<-0
    res$second<-0
        for (i in 1:nrow(mat)){
            for (k in 1:nrow(mat)){
                if (matr[i,k]==maxd){
                    res$first<-i
                    res$second<-k
                    }
                }
             }
        x1<-mat[res$first,1]
        y1<-mat[res$first,2]
        x2<-mat[res$second,1]
        y2<-mat[res$second,2]
        res$m<-(y2-y1)/(x2-x1)
        res$q<-((x1*y1-x1*y2)/(x2-x1))-y1
        }
    if (type=="points"){
        x1<-mat[X1,1]
        y1<-mat[X1,2]
        x2<-mat[X2,1]
        y2<-mat[X2,2]
        res$first<-X1
        res$second<-X2
        res$m<-(y2-y1)/(x2-x1)
        res$q<-((x1*y1-x1*y2)/(x2-x1))-y1
    }
    if (type=="regression"){
        regr<-lm(mat[,X1]~mat[,X2])
        res$m<-regr$coefficients[2]
        res$q<-regr$coefficients[1]
    }
return(res)
}
