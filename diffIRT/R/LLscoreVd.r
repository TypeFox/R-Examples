  LLscoreVd=function(pars,rt,x,N,nit,ai,vi,ter,sd2A,sd2V,model,A,W,nq,eps,delta){
    Vp=(pars)
    result=.C("LLfacV",as.double(rt),as.double(x),as.integer(N),as.integer(nit),as.double(ai),
                    as.double(vi),as.double(sd2A),as.double(sd2V),as.double(ter),as.double(Vp),as.integer(model),
                    as.double(A),as.double(W),as.integer(nq),as.double(eps),
                    as.double(-1),as.double(matrix(0,N)))
    return(result[[16]])
  }