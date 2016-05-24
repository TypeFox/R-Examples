coxseisim <-
function(parreg,parg,
                      lmd0=function(tt)(1+0.5*cos(2*pi*tt)),
                      g=function(x,parg){
                        ifelse(x<=0,0,parg[1]*parg[2]*exp(-parg[2]*x))
                      },
                      censor=1, m=2, trace=TRUE,
                      Z=function(x)matrix(0,length(x),length(parreg))
                      ){
  ##   browser()
  ets <- numeric(0);#event times
  et <- 0; # initialize last event time to 0;
  dur <- 0; # duration initialized to 0;
  ne <- 0; # no. of events
  Zmat <- matrix(numeric(0),ncol=length(parreg))
  while(et < censor){
    int.at.lasttm <- function(x){
      lmd0(x +et)*exp(Z(x+et)%*%parreg + if(ne>0){
        ## exp(apply(outer(x,ets[,tail(seq_along(ets),m)],"-"),1,g,parg=parg))
        rowSums(g(outer(x+et,ets[tail(1:ne,m)],"-"),parg=parg))
      }else rep(0,length(x)) )
    }
    dur <- RND(1,int.at.lasttm);
    et <- et + dur;
    ets <- c(ets,et);
    Zmat <- rbind(Zmat,Z(et));
    ne <- ne+1; if(trace)cat("ne=",ne,"et=",et,"\n")
  }
  ## ne <- length(ets);
  if(ets[ne]>censor){
    ets[ne] <- censor;
    Zmat[ne,] <- Z(censor)
  }else{
    ets <- c(ets, censor)
    Zmat <- rbind(Zmat,Z(censor))
    ne <- ne+1; 
  }
  data.frame(Y=ets,delta=rep.int(1:0,c(ne-1,1)),Z=Zmat)
}

