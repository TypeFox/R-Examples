dbind <-
function(a,b){
        out1<-cbind(a,matrix(0,nrow(a),ncol(b)))
        out2<-cbind(matrix(0,nrow(b),ncol(a)),b)
        out<-rbind(out1,out2)
        out
}
