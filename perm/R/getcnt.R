`getcnt` <-
function(nodehk,cnt.edge,edgesize){
    ## to do: perhaps save time by rewriting this function in C
    ## may speed up exact.network algorithm
    out<-rep(NA,sum(edgesize[nodehk]))
    cnt<-1
    for (i in 1:length(nodehk)){
        cnt2<-cnt+edgesize[nodehk[i]]
        out[cnt:(cnt2-1)]<-cnt.edge[nodehk[i]]:(cnt.edge[nodehk[i]+1]-1)
        cnt<-cnt2
    }
    out
}

