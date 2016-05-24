synth.peak.error <- function(base=0.07,base.time=6, rise.time=5,
rise.factor, rise.factor2, recession.const=0.2, length.out=240,
rez.time=length.out-base.time-rise.time, err1.factor=c(1.2,1.4,1.6),
err2.factor = c(0.01,0.02,0.04), err3.factor=c(2,4,8), err4.factor =
c(9,18,27), err5.factor = c(0.1,0.2,0.4),err6.factor =c(1.5,2,3),  err9.factor=c(2,3,4.5)){

    n.errors <- 9
    n.levels <- 6
    peaks<-array(dim=c(2,n.errors,n.levels, length.out))

    ref.peak <- synth.peak(base=base, base.time=base.time, rise.time=rise.time, rise.factor=rise.factor, recession.const=recession.const, length.out=length.out, rez.time=rez.time)


    #peak flow too high, recession equally fast
    j=1
    for(factor in err1.factor){
        peaks[1,1,4-j,]<-synth.peak(rise.factor=factor*rise.factor, base=base, base.time=base.time, rise.time=rise.time, recession.const=recession.const, length.out=length.out, rez.time=rez.time)
        #peak flow too low, recession equally fast
        if(factor>=rise.factor){
           warning("In error type 1, factor larger than rise.factor.  Adjusting for underestimation")
           factor=0.95*rise.factor
        }
        peaks[1,1,3+j,]<-synth.peak(rise.factor=1/factor*rise.factor, base=base, base.time=base.time, rise.time=rise.time, recession.const=recession.const, length.out=length.out, rez.time=rez.time)
        peaks[2,1,4-j,]<- ref.peak
        peaks[2,1,3+j,]<- ref.peak
        j=j+1
    }

    #too high values in general
    j=1
    for(factor in err2.factor){
        if(factor >= base){
            warning("In peak type 2. err2.factor larger than base flow. Adjusting factor to avoid unexpected results")
            factor <- 0.95 * base
        }
        peak.hight.from.base <- base*(rise.factor - 1)
        rise.factor.adj <- 1+peak.hight.from.base / (base+factor)
        peaks[1,2,4-j,] <- synth.peak(base=base+factor, base.time=base.time, rise.time=rise.time, rise.factor=rise.factor.adj, recession.const=recession.const, length.out=length.out, rez.time=rez.time)
        rise.factor.adj <- 1+peak.hight.from.base / (base-factor)
        peaks[1,2,3+j,]<-synth.peak(base=base-factor, base.time=base.time, rise.time=rise.time, rise.factor=rise.factor.adj, recession.const=recession.const, length.out=length.out, rez.time=rez.time)
        peaks[2,2,4-j,]<- ref.peak
        peaks[2,2,3+j,]<- ref.peak
        j=j+1
    }

    #faster recession
    j=1
    for(factor in err3.factor){
        peaks[1,3,4-j,]<-synth.peak(base=base, base.time=base.time, rise.time=rise.time, rise.factor=rise.factor, recession.const=recession.const/factor, length.out=length.out, rez.time=rez.time)
        peaks[1,3,3+j,]<-synth.peak(base=base, base.time=base.time, rise.time=rise.time, rise.factor=rise.factor, recession.const=recession.const*factor, length.out=length.out, rez.time=rez.time)
        if(diff(peaks[1,3,3+j,c(length.out-1, length.out)])==0){
            warning("In peak type 3: very fast recession. Consider reducing err3.factor")
        }
        peaks[2,3,4-j,]<- ref.peak
        peaks[2,3,3+j,]<- ref.peak
        j=j+1
    }

    #time lag
    j=1
    for(factor in err4.factor ){
        peaks[1,4,4-j,]<-synth.peak(base=base, base.time=base.time+factor, rise.time=rise.time, rise.factor=rise.factor, recession.const=recession.const, length.out=length.out)
        peaks[1,4,3+j,]<-synth.peak(base=base, base.time=base.time-factor, rise.time=rise.time, rise.factor=rise.factor, recession.const=recession.const, length.out=length.out)
        peaks[2,4,4-j,]<- ref.peak
        peaks[2,4,3+j,]<- ref.peak
        j=j+1
    }

    #too high/low but correct integral

    plot.max=max(ref.peak)*2
    to_opt<-function(par,factor){
        test.peak <-synth.peak(base=base, base.time=base.time, rise.time=rise.time, rise.factor=factor*rise.factor, recession.const=par,  length.out=length.out, rez.time=rez.time)
        #plot(test.peak,type="l", xlab="time",ylab="specific discharge/mm/h", lty=2, ylim=c(0,plot.max))
        #lines(ref.peak)
        return(abs(sum(test.peak)-sum(ref.peak)))
    }
    j=1
    for(factor in err5.factor){
        r.opt<-optim(c(recession.const), to_opt, method="L-BFGS-B", lower=c(0), factor=(1+3*factor))
        k_rez<-r.opt$par
        peaks[1,5,4-j,]<-synth.peak(base=base, base.time=base.time, rise.time=rise.time, rise.factor=(1+3*factor)*rise.factor, recession.const=k_rez, length.out=length.out, rez.time=rez.time)
        r.opt<-optim(c(0.2), to_opt, method="L-BFGS-B", lower=c(0), factor=1/(1+factor))
        k_rez<-r.opt$par
        peaks[1,5,3+j,]<-synth.peak(base=base, base.time=base.time, rise.time=rise.time, rise.factor=1/(1+factor)*rise.factor, recession.const=k_rez, length.out=length.out, rez.time=rez.time)
        peaks[2,5,4-j,]<- ref.peak
        peaks[2,5,3+j,]<- ref.peak
        j=j+1
    }

    print("Check if peak volumes are correct for 'volume-optimized'
peaks:")
    print(paste("Volume for reference peak: ", sum(ref.peak)))
    print("Volumes for error peaks:")
    print(rowSums(peaks[1,5,,]))


    #too wide/narrow
    j=1
    for(factor in err6.factor){
        time.change=round(rise.time*factor,0)
        peaks[1,6,4-j,]<-synth.peak(base=base, base.time=base.time-time.change, rise.time=rise.time+time.change, rise.factor=rise.factor, recession.const=recession.const/factor, length.out=length.out, rez.time=rez.time)
        time.change=round(rise.time/factor,0)
        peaks[1,6,3+j,]<-synth.peak(base=base, base.time=base.time+time.change, rise.time=rise.time-time.change, rise.factor=rise.factor, recession.const=recession.const*factor, length.out=length.out, rez.time=rez.time)
        peaks[2,6,4-j,]<- ref.peak
        peaks[2,6,3+j,]<- ref.peak
        j=j+1
    }

    #Referenzpeaks für weitere peak Fehler
    recession.const2=0.1*recession.const
    base2=0.77*base
    #ref.peak3 <-rnorm(240, sd=0.003)+synth.peak(base=base2, base.time=-0.5, rise.time=0, rise.factor=rise.factor2, recession.const=recession.const2)
    ref.peak4<- synth.peak(base=base2, base.time=-0.5, rise.time=0, rise.factor=rise.factor2, recession.const=recession.const2, length.out=length.out)

    #modeled peak, shifted up and down around the measured recession, where no peak exists

    for(j in 1:6){
        peaks[1,7,j,]<-peaks[1,2,j,]
        peaks[2,7,j,]<- ref.peak4
    }

    # abklingende Rezession über oder unterschätzt
    j=1
    for(factor in err9.factor){
        peaks[1,9,4-j,]<- synth.peak(base=base2, base.time=-0.5, rise.time=0, rise.factor=rise.factor2+0.1*factor, recession.const=recession.const2/factor, length.out=length.out)
        if(0.1*factor +1 >= rise.factor2){
             warning("In peak type 9. (0.1*factor) >= rise.factor2 - 1.  Adjusting factor to avoid unexpected results for underestimation")
             factor = 9.5 * (rise.factor2-1)
        }
        peaks[1,9,3+j,]<- synth.peak(base=base2, base.time=-0.5, rise.time=0, rise.factor=rise.factor2-0.1*factor, recession.const=recession.const2*factor, length.out=length.out)
        peaks[2,9,4-j,]<- ref.peak4
        peaks[2,9,3+j,]<- ref.peak4
        j=j+1
    }

    for(j in 1:6){
        peaks[1,8,j,]<- peaks[1,9,j,]
        peaks[2,8,j,]<- ref.peak
    }
    return(peaks)
}
