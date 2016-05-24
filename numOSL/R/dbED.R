#####
dbED<-
function(EDdata, plot=TRUE, typ=c("pdf","hist"),
         from=NULL, to=NULL, step=NULL, nbin=15,
         pcolor="grey", psize=1.5, outfile=NULL) {
    UseMethod("dbED")
} #
### 2014.10.01; revised in 2016.01.20.
dbED.default<-
function(EDdata, plot=TRUE, typ=c("pdf","hist"),
         from=NULL, to=NULL, step=NULL, nbin=15,
         pcolor="grey", psize=1.5, outfile=NULL) {
    ### Stop if not.
    stopifnot(ncol(EDdata)==2L, nrow(EDdata)>=5L,
              all(EDdata[,2L,drop=TRUE]>0),
              length(plot)==1L, is.logical(plot),
              is.character(typ),
              all(typ %in% c("pdf","hist")),
              is.null(from) || is.numeric(from),
              is.null(to) || is.numeric(to),
              is.null(step) || is.numeric(step), 
              length(nbin)==1L, is.numeric(nbin),
              length(pcolor)==1L, is.character(pcolor),
              length(psize)==1L, is.numeric(psize),
              is.null(outfile) || is.character(outfile))
    ###
    if (!is.null(step)) {
        if (length(step)!=1L) stop("Error: step should be an one-element vector!")
        if (step<=0)  stop("Error: step should exceed zero!")
    } # end if.
    ### 
    if (!is.null(outfile)) {
        if (length(outfile)!=1L)  stop("Error: outfile should be an one-element vector!")
    } # end if.
    ###
    nED<-nrow(EDdata)
    EDdata<-EDdata[order(EDdata[,1L,drop=TRUE],
                   decreasing=FALSE),,drop=FALSE]
    ed1<-as.numeric(EDdata[,1L,drop=TRUE])
    sed1<-as.numeric(EDdata[,2L,drop=TRUE])
    ###
    if (is.null(from)) {
        from<-ifelse(min(ed1)>=0,
                     min(ed1)*0.7,
                     min(ed1)*1.3)
    } # end if
    if (is.null(to)) {
        to<-ifelse(max(ed1)>=0,
                   max(ed1)*1.3,
                   max(ed1)*0.7)
    } # end if
    if (from>=to) {
        stop("Error: from must not exceed to!")
    } # end if
    ###
    ### Calculate weighted skewness and kurtosis.
    weight<-ed1/sed1
    meanED<-mean(ed1)
    sdED<-sd(ed1)
    medianED<-median(ed1)
    skewness<-sum(weight*((ed1-meanED)/sdED)^3L)/sum(weight)     
    Std.skewness<-sqrt(6L/nED)   
    kurtosis<-nED*(nED+1L)/(nED-1L)/(nED-2L)/(nED-3L)*
              sum(((ed1-meanED)/sdED)^4L)-
              3L*(nED-1L)^2L/(nED-2L)/(nED-3L)
    Std.kurtosis<-sqrt(24L/nED)
    ###
    ### Calculate weighted ED.
    weightED<-sum(ed1/sed1^2L)/sum(1.0/sed1^2L)
    sweightED <-1.0/sqrt(sum(1.0/sed1^2L))
    weightED<-c(weightED,sweightED)
    ###
    quantileED<-quantile(ed1, probs=c(0.05,0.1,0.15,
                         0.5,0.85,0.9,0.95))
    ###
    output<-list("weight.ED"=round(weightED, 3L),
                 "skewness"=round(c(skewness, Std.skewness),3L),
                 "kurtosis"=round(c(kurtosis, Std.kurtosis),3L),
                 "quantile.ED"=round(quantileED,3L))
    ###
    if (plot==TRUE) {
        if (typ[1L]=="pdf") {
            if (is.null(step)) { 
                step<-max(diff(sort(ed1)))/10.0
                cat(paste("Default step size:", round(step,3L), "\n\n"))
            } # end if.
            ###
            spreadED<-seq(from=from, to=to, by=step)
            pdfMat<-matrix(nrow=length(spreadED), ncol=nED)
            ###
            for(i in seq(nED)) {
                pdfMat[,i]<-dnorm(x=spreadED, mean=ed1[i], 
                                  sd=sed1[i], log=FALSE)
            } # end if
            pdfED<-rowSums(pdfMat)/sum(rowSums(pdfMat))
            ###
            if (!is.null(outfile))  {
                write.csv(cbind(spreadED, pdfED), 
                  file=paste(outfile,".csv",sep=""))
            } # end if.        
            ### 
            plot(spreadED, pdfED, 
                 main="Equivalent Dose Distribution", 
                 xlab="Equivalent Dose (Gy)", ylab="Density",
                 xlim=c(from,to), type="l", lwd=2)
            xTicks<-axTicks(side=1L)
            maxYx<-spreadED[which.max(pdfED)]
            box(lwd=1)
        } else if (typ[1L]=="hist") {
            breaks<-pretty(ed1, n=nbin)
            HIST<-hist(ed1, breaks=breaks, 
                       main="Equivalent Dose Distribution", 
                       xlab="Equivalent Dose (Gy)", 
                       xlim=c(from,to))
            xTicks<-axTicks(side=1L)
            maxYx<-HIST$mids[which.max(HIST$counts)]
            box(lwd=1)
        } # end if
        ###
        ###
        par("new"=TRUE)
        plot(ed1, seq(nED), xlab="", ylab="", 
             xlim=c(from,to), type="p", 
             pch=23, cex=psize, bg=pcolor, 
             xaxt="n", yaxt="n")
        par("new"=FALSE)
        ###
        options("warn"=-1)
        arrows(ed1-sed1/2L, seq(nED),
               ed1+sed1/2L, seq(nED),
               code=3, lwd=1.5, angle=90, 
               length=0.05, col="black")
        options("warn"=0)
        ###
        legend(ifelse(maxYx>median(xTicks),"left","right"), 
               legend=c(paste("N=",nED,sep=""),
                        paste("mean=",round(meanED,2L)," (Gy)",sep=""),
                        paste("median=",round(medianED,2L)," (Gy)",sep=""),
                        paste("sd=",round(sdED,2L),sep="")), cex=1, bty="n")
        ###
    } # end if
    ###
    return(output)
} # end function dbED.
###           
