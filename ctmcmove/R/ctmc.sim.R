ctmc.sim <- function(Q,start.state=1,T=1,final.state=NA){
    ##
    ## Code to simulate a continuous-time Markov chain on a finite sample space
    ##
    ## Q is an N-by-N  "rate matrix" containing infinitessimal rates {q_ij}
    ## start.state is the starting state (a number between 1 and N)
    ## T is the total time to run the simulation
    ##
    ## Output is:
    ##  ec=embedded chain (sequence of states)
    ##  t=residence time at each state
    
    if(diag(Q)[1]>0){
        Q=-Q
        diag(Q)=0
    }
    if(diag(Q)[1]<0){
        diag(Q)=0
    }
    
    N=nrow(Q)
    v=apply(Q,1,sum)
    P=Q/v
    X=start.state
    current.state=start.state
    t=0
    if(is.na(final.state)){
        while(max(t)<T){
            wait.time=rexp(1,v[current.state])
            current.state=sample(N,1,prob=P[current.state,])
            X=c(X,current.state)
            t=c(t,max(t)+wait.time)
        }
    }else{
        while(current.state!=final.state){
            wait.time=rexp(1,v[current.state])
            current.state=sample(N,1,prob=P[current.state,])
            X=c(X,current.state)
            t=c(t,max(t)+wait.time)
        }
    }
    list(rt=t,ec=X)
}

