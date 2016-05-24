"permcont"<-function(Table){

# Giraudoux 13.7.2004
# return a random permutation of contingency table 
# n rows x 2 columns keeping the marginal totals
    
    # transform the contingency table in a 1 vector
    # of categories (cate) and a vector of 0/1 (eff)
    col1<-Table[,1]; col2<-Table[,2]
    totcol1<-sum(col1);totcol2<-sum(col2);totlig<-col1+col2
    for (i in 1:length(col1)){
        if (i==1){cate<-rep(i,totlig[i])} else {cate<-c(cate,rep(i,totlig[i]))}
        if (i==1)eff<-c(rep(1,col1[i]),rep(0,col2[i]))
        else eff<-c(eff,rep(1,col1[i]),rep(0,col2[i]))
    }

    # permute the 0/1 vector and return the tableau
    eff<-sample(eff)
    tt<-table(cate,eff)
    n1<-tt[,2];n2<-tt[,1]
    tt<-cbind(n1,n2)
    tt
}
