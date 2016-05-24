lagIngray <- structure(
    function#First-local lag
    ###This function can compute the lag of the first local on the
    ###auto-correlation function (acf) of smoothed grays.
    (
        image,##<<character or matrix. Either path of an image section
        ##or an array representing a gray matrix.
        acf=FALSE, ##<<logical. If TRUE the output is extended with the
        ##acf.
        ...##<< arguments to be passed to \code{\link{imageTogray}}.
    )
    {

        gray <- image
        if(is.character(gray))
            gray <- imageTogray(gray,...)
        
        average.ingray <- function(gray){
        e <- apply(gray, 2, function(x) exp(mean(log(x))))
        data2 <- ts(e,end=ncol(gray))
        return(data2)}
        ser1 <- average.ingray(gray)
        acfs <- acf(ser1,lag.max=dim(gray)[2]/3,plot=FALSE)
        if(acf)acf1 <- acf(ser1,lag.max=dim(gray)[2],plot=FALSE)
        acfs1 <- data.frame(acf=acfs$acf)
        acfs1[,'dif'] <- c(NA,diff(acfs1$acf))                   
        f.lab <- function(i,data){
            ifelse(data[i,'dif']>0,1,0)}
        acfs1[,'flags'] <- sapply(1:nrow(acfs1),f.lab,acfs1)
        acfs1[,'cumflags'] <- c(NA,cumsum(acfs1$flags[2:nrow(acfs1)]))
        cumacf <- cumsum(acfs1$flags[2:nrow(acfs1)])
        Mode <- function(x) {
            ux <- unique(x)
            ux[which.max(tabulate(match(x,ux)))]}
        lim1 <- Mode(cumacf[cumacf>0])
        acfs1[,'window2'] <- with(acfs1,ifelse(cumflags>0&cumflags<=lim1,2,NA))
        window2 <- na.omit(acfs1)                   
        maxacf <- with(window2,max(acf))                   
        local1 <- as.numeric(rownames(window2[window2$acf==maxacf,]))
        if(acf){return(list(local = as.numeric(local1),acf = data.frame(acf=acf1$acf)))}
        else{return(local1)}
        ### constant value of the first local on the acf of the
        ### smoothed gray. If acf is TRUE then the computed acf is
        ### added to the output (see \code{\link{linearDetect}}, and
        ### \code{\link{graySmoothed}}).
    }
,
    ex=function(){
        ## (not run) Read one image sample in folder of package measuRing:
        image1 <- system.file("P105_a.tif", package="measuRing")
        ##(not run) First local in the acf of smoothed grays:       
        local1 <- lagIngray(image1,acf = TRUE)        
        ##(not run) Plot of first local over the acf: 
        Flocal <- local1[['local']]
        Clocal <- ts(local1[['acf']][Flocal,],start=Flocal)
        acf <- ts(local1[['acf']],start=1)    
        {plot(acf,type='h',col='gray',xlab='Lag',main='First local lag')
        points(Clocal,pch=19,cex=0.5)}
    }
)
