lis <-
function( x ){
    N<-length(x);
    if(N==1){rr<-1}
    if(N>1){
        gr<-c(1:N)*0
        gr[1]<-1
        for ( i in 2:N ){
            gr[i] = 1;
            for ( j in 1:(i - 1) ){ 
                if ( x[i] > x[j] ) { gr[i] = max( gr[i], gr[j] + 1) } 
            }
        }
        rr<-max( gr )
    }
    result<-rr
}
