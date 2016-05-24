
gflasso<-function(A,y,tp,s1,s2,RmaxIter=100,RaMaxIter=1000,Rrho=5,Rtau=0.15
                  ,Rwt=rep(1,length(tp)),Rtol=1e-3,RaTol=1e-3,x0=rep(0,ncol(A))){
  opts_n<-nrow(A)
  opts_p<-ncol(A)
  opts_g<-length(tp)
  x<-x0
  A<-unlist(A)
  y<-unlist(y)
  tp<-unlist(tp)
  result<-.C("do_gflasso",as.double(x),as.double(A),as.double(y),as.integer(tp),as.double(Rwt),
             as.double(s1),as.double(s2),as.integer(opts_n),as.integer(opts_p),as.integer(opts_g),
             as.integer(RmaxIter),as.integer(RaMaxIter),as.integer(Rrho),as.double(Rtau),
             as.double(Rtol),as.double(RaTol),PACKAGE="FGSG")[[1]]
  return(list(weight=result))
}


goscar<-function(A,y,tp,s1,s2,RmaxIter=100,RaMaxIter=1000,Rrho=5,Rtau=0.15
                 ,Rwt=rep(1,length(tp)),Rtol=1e-3,RaTol=1e-3,x0=rep(0,ncol(A))){
  opts_n<-nrow(A)
  opts_p<-ncol(A)
  opts_g<-length(tp)
  x<-x0
  A<-unlist(A)
  y<-unlist(y)
  tp<-unlist(tp)
  result<-.C("do_goscar",as.double(x),as.double(A),as.double(y),as.integer(tp),as.double(Rwt),
             as.double(s1),as.double(s2),as.integer(opts_n),as.integer(opts_p),as.integer(opts_g),
             as.integer(RmaxIter),as.integer(RaMaxIter),as.integer(Rrho),as.double(Rtau),
             as.double(Rtol),as.double(RaTol),PACKAGE="FGSG")[[1]]
  return(list(weight=result))
}


ncFGS<-function(A,y,tp,s1,s2,RmaxIter=100,RaMaxIter=1000,Rrho=5,Rtau=0.15
                ,Rwt=rep(1,length(tp)),Rtol=1e-3,RaTol=1e-3,x0=rep(0,ncol(A))){
  opts_n<-nrow(A)
  opts_p<-ncol(A)
  opts_g<-length(tp)
  x<-x0
  A<-unlist(A)
  y<-unlist(y)
  tp<-unlist(tp)
  result<-.C("do_ncFGS",as.double(x),as.double(A),as.double(y),as.integer(tp),as.double(Rwt),
             as.double(s1),as.double(s2),as.integer(opts_n),as.integer(opts_p),as.integer(opts_g),
             as.integer(RmaxIter),as.integer(RaMaxIter),as.integer(Rrho),as.double(Rtau),
             as.double(Rtol),as.double(RaTol),PACKAGE="FGSG")[[1]]
  return(list(weight=result))
}



ncTFGS<-function(A,y,tp,s1,s2,RmaxIter=100,RaMaxIter=1000,Rrho=5,Rtau=0.15
                 ,Rwt=rep(1,length(tp)),Rtol=1e-3,RaTol=1e-3,x0=rep(0,ncol(A))){
  opts_n<-nrow(A)
  opts_p<-ncol(A)
  opts_g<-length(tp)
  x<-x0
  A<-unlist(A)
  y<-unlist(y)
  tp<-unlist(tp)
  result<-.C("do_ncTFGS",as.double(x),as.double(A),as.double(y),as.integer(tp),as.double(Rwt),
             as.double(s1),as.double(s2),as.integer(opts_n),as.integer(opts_p),as.integer(opts_g),
             as.integer(RmaxIter),as.integer(RaMaxIter),as.integer(Rrho),as.double(Rtau),
             as.double(Rtol),as.double(RaTol),PACKAGE="FGSG")[[1]]
  return(list(weight=result))
}


ncTL<-function(A,y,tp,s1,s2,RmaxIter=100,RaMaxIter=1000,Rrho=5,Rtau=0.15
               ,Rwt=rep(1,length(tp)),Rtol=1e-3,RaTol=1e-3,x0=rep(0,ncol(A))){
  opts_n<-nrow(A)
  opts_p<-ncol(A)
  opts_g<-length(tp)
  x<-x0
  A<-unlist(A)
  y<-unlist(y)
  tp<-unlist(tp)
  result<-.C("do_ncTL",as.double(x),as.double(A),as.double(y),as.integer(tp),as.double(Rwt),
             as.double(s1),as.double(s2),as.integer(opts_n),as.integer(opts_p),as.integer(opts_g),
             as.integer(RmaxIter),as.integer(RaMaxIter),as.integer(Rrho),as.double(Rtau),
             as.double(Rtol),as.double(RaTol),PACKAGE="FGSG")[[1]]
  return(list(weight=result))
}


ncTF<-function(A,y,tp,s1,s2,RmaxIter=100,RaMaxIter=1000,Rrho=5,Rtau=0.15
               ,Rwt=rep(1,length(tp)),Rtol=1e-3,RaTol=1e-3,x0=rep(0,ncol(A))){
  opts_n<-nrow(A)
  opts_p<-ncol(A)
  opts_g<-length(tp)
  x<-x0
  A<-unlist(A)
  y<-unlist(y)
  tp<-unlist(tp)
  result<-.C("do_ncTF",as.double(x),as.double(A),as.double(y),as.integer(tp),as.double(Rwt),
             as.double(s1),as.double(s2),as.integer(opts_n),as.integer(opts_p),as.integer(opts_g),
             as.integer(RmaxIter),as.integer(RaMaxIter),as.integer(Rrho),as.double(Rtau),
             as.double(Rtol),as.double(RaTol),PACKAGE="FGSG")[[1]]
  return(list(weight=result))
}



ncTLF<-function(A,y,tp,s1,s2,RmaxIter=100,RaMaxIter=1000,Rrho=5,Rtau=0.15
                ,Rwt=rep(1,length(tp)),Rtol=1e-3,RaTol=1e-3,x0=rep(0,ncol(A))){
  opts_n<-nrow(A)
  opts_p<-ncol(A)
  opts_g<-length(tp)
  x<-x0
  A<-unlist(A)
  y<-unlist(y)
  tp<-unlist(tp)
  result<-.C("do_ncTLF",as.double(x),as.double(A),as.double(y),as.integer(tp),as.double(Rwt),
             as.double(s1),as.double(s2),as.integer(opts_n),as.integer(opts_p),as.integer(opts_g),
             as.integer(RmaxIter),as.integer(RaMaxIter),as.integer(Rrho),as.double(Rtau),
             as.double(Rtol),as.double(RaTol),PACKAGE="FGSG")[[1]]
  return(list(weight=result))
}
