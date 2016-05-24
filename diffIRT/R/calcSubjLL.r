 calcSubjLL=function(pars,rt,x,N,nit,model,A,W,nq,eps,delta,start,constrain){
    par=rep(999,3*nit+2)
    par[constrain!=0]=pars[constrain]
    par[constrain==0]=start[constrain==0]
    ai=par[1:nit]
    vi=par[(nit+1):(2*nit)]
    ter=par[(2*nit+1):(3*nit)]
    sd2A=par[3*nit+1]
    sd2V=par[3*nit+2]
    result=.C("LLdiff",as.double(rt),as.double(x),as.integer(N),as.integer(nit),as.double(ai),
                    as.double(vi),as.double(sd2A),as.double(sd2V),as.double(ter),as.integer(model),
                    as.double(A),as.double(W),as.integer(nq),as.double(eps),
                    as.double(-1),as.double(matrix(0,N)))
    return(result[[16]])
  }
  
  