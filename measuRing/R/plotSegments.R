plotSegments <- structure(
    function #Image segments
###One or several plots of consecutive segments of gray matrix and
###smoothed grays.
    (
         image,##<<character or matrix. Either path of an image
                ##section or an array representing a gray matrix.
        ratio = NULL,##<<NULL or vector with two values representing the
                     ##aspect of the plots (height, and width). If
                     ##NULL the default aspect in \code{par()} is
                     ##used.
        marker = NULL, 	##<<NULL or a number from 1 to 10 as explained
                        ##in \code{\link{colNarrow}}. If numeric then
                        ##three kind of markers are indicated: those
                        ##narrow rings with averages major than
                        ##\code{marker}, chronological markers
                        ##(decades, centuries, and millenia), and the
                        ##column numbers in gray matrix of the ring
                        ##borders. The markers are highlighted with
                        ##the color in \code{col.marker}. If NULL no
                        ##markers are highlighted.
        col.marker = 'red', ##<< color of the markers.
        tit = TRUE, ##<<logical or character. A title for the plots. If
                   ##TRUE the main title is the image name. For more
                   ##than 1 segment the main title ends with the
                   ##segment number.
        plot = TRUE,##<<logical. If TRUE the image segments are
                    ##plotted.
        ...##<< arguments to be passed to \code{\link{dataSegments}}.


        
    )
    {

        ffig <- function(gray,pixtypes,TRW,origin,
                         segs,ratio,marker,col.marker,tit){

            f.rown <- function(x)as.numeric(rownames(x))            
            range. <- function(n)c(0,n)
            gray <- t(gray)
            mono.chr <- c('white','gray70','gray50','gray30','black')
            poly.chr <- c('orangered')
            poly.chr <- col.marker
            cex. <- par('cex')*0.7
            asp. <- 0
            rowpix <- f.rown(pixtypes)
            included <- pixtypes[pixtypes[,"borders"],]
            rowinc <- f.rown(included)
            xrange <- range(rowpix)
            yrange <- range.(ncol(gray))
            rowt <- nrow(TRW)
            
            f.deca <- function(TRW,decades){
               span <- decades*10
               rangeyr <- range(TRW[,1])
               move <- rangeyr[1]%%span
               decs <- seq(rangeyr[1]-move,rangeyr[2],span)
               decs <- decs[decs%in%rangeyr[1]:rangeyr[2]]
               pixn <- rownames(TRW[TRW[,1]%in%decs,])
               rowns <- f.rown(TRW)
               deca <- rowns + c(diff(rowns)/2,NA)
               marking <- data.frame(deca)
               rownames(marking) <- rownames(TRW)
               decas <- marking[pixn,]
               return(decas)}
            
            f.series <- function(pixtypes,xlab.,ylab.){
               with(pixtypes,
               plot(rowpix,pixtypes[,1],
                 xlim=xrange,axes=FALSE,type='n',xlab=xlab.,ylab=ylab.,
                 cex.lab = 0.8,col.lab=mono.chr[4]))}
            
            f.image <- function(pixtypes,xlab.,ylab.,...){
                with(pixtypes,
                plot(rowpix,rep(yrange[2]/2,length(rowpix)),
                  xlim=xrange,ylim=yrange,
                     axes=FALSE,type='n',
                  xlab=xlab.,ylab=ylab.,asp = asp.,
                     cex.lab=0.8,col.lab=mono.chr[4]))}
            
            f.nodes <- function(pixtypes){
               par(mfg=c(1,1));f.series(pixtypes,'','')
               par(mfg=c(2,1));f.image(pixtypes,'','')
               lim2 <- grconvertY(yrange[2]/2,'user','ndc')
               par(mfg=c(1,1));f.series(pixtypes,'','')
               lim1 <- grconvertY(lim2,'ndc','user')
               with(included,segments(rowinc,included[,1],rowinc,lim1,
               col=mono.chr[2],lty=1,lwd=0.5))}
            
            if(rowt!=0)
            f.nodes(pixtypes)
            par(mfg=c(1,1));f.series(pixtypes,'','Smoothed gray')            
            with(pixtypes,
            lines(rowpix,pixtypes[,1],xlim=xrange,
            lwd=1,type='l',col=mono.chr[3]))
            if(rowt!=0&!is.null(marker))
             lines(xrange,rep(origin,2),lwd=1,lty=3,col=poly.chr[1])
            
            f.labels1 <- function(pixtypes,marker){#
               with(included,
               points(rowinc,included[,1],cex=0.3,pch=19,col=mono.chr[5]))
               
            frect <- function(x1,y1,label,cex.,exp.){
               cha <- paste(" ", label, " ", sep = "")
               ## cha <- label
               xusr <- par("usr")
               xh <- strwidth(cha, cex = cex.)
               yh <- strheight(cha, cex = cex.)#*4/3
               tmp <- xh
               xh <- yh/(xusr[4]-xusr[3])*par("pin")[2]
               xh <- xh/ par("pin")[1] * (xusr[2]-xusr[1])
               yh <- tmp/(xusr[2]-xusr[1])* par("pin")[1]
               yh <- yh/ par("pin")[2] * (xusr[4]-xusr[3])
               plus. <- nchar(cha)*strheight(x1,units='user',cex=cex.)*exp.
               y1 <- y1 + plus.
               rect(x1 - xh/2, y1 - yh/2, x1 + xh/2,y1 + yh/2,
               border=NA, col = "white")
               text(x1, y1, cha, cex = cex.,adj=c(0.5,0.5),
                col=poly.chr[1],srt = 90)}
                if (!is.null(marker))
                with(included,frect(rowinc,included[,1],rowinc,cex.,0.6))}
            
            if(rowt!=0){
               f.labels1(pixtypes,marker)}
            axis(2,col= mono.chr[3],cex.axis = 0.7,cex.main = 0.7)
            par(mfg=c(2,1));f.image(pixtypes,'Gray-matrix column','')#,col.lab =colflag2)
            rasterImage(t(gray),
              xrange[1]*par('cex'),yrange[1]*par('cex'),
              xrange[2]*par('cex'),yrange[2]*par('cex'),interpolate=FALSE)

;            axis(1,col= mono.chr[3],cex.axis = 0.7,cex.main = 0.7)
            if(rowt==0){
                with(pixtypes,
            lines(xrange,rep(yrange[2]/2,2),lwd=1,lty=2,col=poly.chr[1]))}

            f.labels2 <- function(included,marker){
               with(included,
               segments(rowinc,rep(ncol(gray)/2,nrow(included)),
               rowinc,rep(ncol(gray), nrow(included)),
               col='honeydew2',lty=1,lwd=0.5))
               with(included,
               points(rowinc,rep(yrange[2]/2,nrow(included)),
               cex=0.3,pch=19,col=mono.chr[5]))}
            
           if(rowt!=0){
                
            f.labels2(included,marker)                
            if(!is.null(marker)){
              thou. <- f.deca(TRW,100)
              hund. <- f.deca(TRW,10)
              fifth. <- f.deca(TRW,5)
              deca. <- f.deca(TRW,1)
            if(!identical(thou.,numeric(0)))
              hund. <- setdiff(hund.,thou.)
              
             if(!identical(hund.,numeric(0)))
                fifth. <- setdiff(fifth.,hund.)
             if(!identical(fifth.,numeric(0)))
                deca. <- setdiff(deca.,fifth.)             
                points(deca.,rep(3*(yrange[2]/4),length(deca.)),
                       cex=0.30,pch=19,col=poly.chr[1])
                hundfif <- c(hund.,fifth.)
                hundfif <- setdiff(hundfif,thou.)
              
              fpmarker <- function(decs,frgray){
                points(decs,rep(yrange[2]*frgray,length(decs)),
                cex=0.3,pch=19,col=poly.chr[1])}
              
              ## thoum <- par('cex')*c(11:14)/16
              thoum <- par('cex')*c(23:26)/32

              for(i in 1:length(thoum)){fpmarker(thou.,thoum[i])}
              ## hundum <- par('cex')*c(11,13)/16
              hundum <- par('cex')*c(23,25)/32

              for(i in 1:length(hundum)){fpmarker(hundfif,hundum[i])}
              included. <- as.numeric(colNarrow(TRW,marker))
              incum <- c(-5,5)#*par('cex')
              for(j in 1:length(incum))
                  fpmarker((included. + incum[j]),15/16)}}        

            par(mfg=c(1,1));f.series(pixtypes,'','')

            fcyear <- function(pixtypes){
              par(mfg=c(2,1));f.image(pixtypes,'','')
              lim2 <- grconvertY(yrange[2]/2,'user','ndc')
              par(mfg=c(1,1));f.series(pixtypes,'','')
              lim1 <- grconvertY(lim2,'ndc','user')
              range. <- lim1/2
              return(range.)}

            c.year <- fcyear(pixtypes)
            trwh <- paste(' ',as.character(TRW[1,'year']),' ',sep='')
            trwh. <- nchar(as.character(TRW[1,'year']))
            plus.. <- trwh.*strheight(trwh,units='user',cex=0.6)

            if(rowt!=0){
               segments(f.rown(TRW),rep(c.year - 0.5*plus..,nrow(TRW)),
               f.rown(TRW),rep(c.year + 0.5*plus..,nrow(TRW)),
               col='white',lty=1,lwd=1,lend='square')# 
               with(TRW,
               text(f.rown(TRW),rep(c.year ,nrow(TRW)),
               year,cex=cex.,adj=c(0.5,0.5),col=mono.chr[4],srt=90))}
               
            par(mfg=c(2,1));f.image(pixtypes,'','')
            
            par(oma=c(1.5,1,0,0))
            libname <- 'measuRing V0.3; Printing time:'
            mtext(paste(libname,
              " ",format(Sys.time(), "%Y-%m-%d %H:%M")),
               cex=0.5, line=0, side=SOUTH<-1,
               adj=0, col=mono.chr[2],outer=TRUE)
            par(oma=c(0,0,3,0))

            if(segs>1&is.character(tit))
                {tit <- paste(tit,': ',k)}
            title(tit,outer=TRUE,cex.main=1,col.main=mono.chr[3])
        }
                
        options(warn = -1)
        segments. <- dataSegments(image,...)

        f.rown <- function(x)as.numeric((rownames(x)))

        f.tit <- function(image){
            p <- '.tif'
            if(any(grepl('.png',image)))p <- '.png'
            bn <- basename(image)
            gsub(p,'',bn)}

        opar <- list(mfrow=c(2,1),mar=c(3,2,3,2),oma=c(2,3,0,0),xpd=NA)
        if(plot){
        graphics.off()
        if(is.character(tit)){tit <- tit}
        else{if(tit)                 
                 {tit <- f.tit(attributes(segments.)[['image']])}
             else{tit <- NULL}}
        origin <- attributes(segments.)[['origin']]
        segs <- length(segments.[[2]])
        for(k in segs:1){
            if(is.null(ratio)){dev.new()}
            else{dev.new(width = ratio[1],height = ratio[2])}
            par(opar)
         ffig(segments.[['imageTogray']][[k]],
              segments.[['ringBorders']][[k]],
              segments.[['ringWidths']][[k]],
              origin,segs,ratio,marker,col.marker,tit)}}
        attributes(segments.)[['opar']] <- opar        
        return(segments.)
        ###the image segments and a list such as that produced by
        ###\code{\link{dataSegments}}.
    }
   ,
    ex=function(){
        ## (not run) Read one image sample in folder of package measuRing:
        image1 <- system.file("P105_a.tif", package="measuRing")        
        ## column numbers to be included/avoided:
        Toinc <- c(196,202,387,1564) 
        Toexc <- c(21,130,197,207,1444,1484)        
        ##(not run) Plotting of five image segments:
        plots <- plotSegments(image1,rgb=c(0.5,0,0.5),last.yr=2011,
            marker=8,segs=3,inclu = Toinc,exclu = Toexc)
        ## plots <- plotSegments(rwidths,segs = 4,marker=8)
        ## (not run) kill all the image segments:
        graphics.off()
    
    }
)
