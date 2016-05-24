 
## The function is currently defined as
mat.plot.mc<-function (mcObj) 
{
#> names.f(all.sawada)
#     [,1] [,2]           
#[1,] "1"  "sqdist"       
#[2,] "2"  "cumcurve"     
#[3,] "3"  "cutoffs"      
#[4,] "4"  "method"       
#[5,] "5"  "samplesize"   
#[6,] "6"  "replacement"  
#[7,] "7"  "probabilities"
#[8,] "8"  "wascounts"   

    windows(width = 8.5, height = 11, pointsize = 10)
    par(mfrow = c(3, 1), pty = "s")
    theTitle=paste("Method =",mcObj$method,"Sample size = ",mcObj$samplesize)
    hist(mcObj$sqdist, xlim = c(0, 2), xlab = "Squared-chord distance (SCD)", 
        ylab = "Frequency", main=theTitle,col = 16)
    par(usr = c(0, 1, 0, 1))
    vec2 <- mcObj$sqdist
    text(0.02, .87, paste("\nMin = ", round(min(vec2, na.rm = T), 
        3), "\n1st Quartile= ", round(quantile(vec2, 0.25), 
        3), "\nMedian = ", round(median(vec2, na.rm = T), 3), 
        "\nMean = ", round(mean(vec2, na.rm = T), 3), "\n3rd Quartile = ", 
        round(quantile(vec2, 0.75), 3), "\nMax = ", round(max(vec2, 
            na.rm = T), 3), "\nSD = ", round(sqrt(var(vec2, 
            na.rm = T)), 3), sep = ""), cex = 1, adj = 0, 
        col = 1)

        sigx=mcObj$cumcurve
        plot(sigx[, 1], sigx[,2], type = "l", xlab = "Squared-chord distance (SCD)", ylab = "Monte-Carlo probability Pr(SCD < p)")
        

        #PLOT THE BLOW-UP OF THE SIGNIFICANCE REGION

        tlen=length(mcObj$cutoffs$x)
        x=mcObj$cutoffs$y
        y=mcObj$cutoffs$x

        plotsub=mcObj$cumcurve[,1]<=(max(x)+0.1)
        newset=mcObj$cumcurve[plotsub,]

        sigx=newset
        lines(sigx[, 1], sigx[,2],lwd=2,col="blue")

        plot(sigx[, 1], sigx[,2], type = "n", xlab = "Squared-chord distance (SCD)", ylab = "Monte-Carlo probability Pr(SCD < p)")
        lines(sigx[, 1], sigx[,2],lwd=2,col="blue")

        for(i in 1:tlen){
                segments(x[i],0,x[i],y[i],col="red")
                segments(0,y[i],x[i],y[i],col="red")
                xt=round(x,3)
                text(x[i],0,xt[i],adj=c(.1,.1),srt=90)
                text(0,y[i],y[i],adj=c(0,0))

        }


  }



