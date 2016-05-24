getclf<-function(data, freq)
  {
    # Takes a sorted data frame and returns a likelihood function to be minimized
    nvars<-ncol(data)
    pars<-double(nvars+nvars*(nvars+1)/2)
    testdata<-data[cumsum(freq),]
    presabs<-ifelse(is.na(testdata),0,1)

    data<-t(data)   # convert data to vectors that can be passed to C
    presabs<-t(presabs)
    dim(presabs)<-NULL
    dim(data)<-NULL
    data<-data[!is.na(data)]
    
#    if (!is.loaded(symbol.C("evallf"))) {
#        cat("loading object code...\n")
#        dyn.load("st771/libs/st771.so")
#    }

    function(pars){
      .C("evallf",as.double(data),as.integer(nvars),as.integer(freq),
         as.integer(x=length(freq)),as.integer(presabs),as.double(pars),val=double(1),PACKAGE="mvnmle")$val;
    }
  }
