Bootstrap.Youden.Cutoff.1run <-
function(x,y,z,method="Normal",seed0=12345,randomStart.N=1,t.minus.start=NULL,t.plus.start=NULL,...)
  {
    ##draw new samples for each diagnosis group
    new.x <- bootSample(x,seed0)
    new.y <- bootSample(y,seed0+1)
    new.z <- bootSample(z,seed0+2)
    ##calculate group means
    x.mean <- mean(new.x,na.rm=TRUE)
    y.mean <- mean(new.y,na.rm=TRUE)
    z.mean <- mean(new.z,na.rm=TRUE)
    
    if(x.mean<=y.mean && y.mean<=z.mean)
      {###if mean are in increasing order from D- to D0 to D+, estimate Youden and cut-point
        res <- Youden3Grp.PointEst(x=new.x,y=new.y,z=new.z, method=method,randomStart.N=randomStart.N,t.minus.start=t.minus.start,t.plus.start=t.plus.start,...)$est
      }else{
        warning("##the group means from D- to D0 to D+ must be in increasing order!")
        res <- rep(NA,17)##the length of the data frame component from Youden3Grp.PointEst.
    }    
    return(res)    
  }

