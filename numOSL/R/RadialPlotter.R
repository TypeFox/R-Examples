#####
RadialPlotter<-
function(EDdata, ncomp=0, addsigma=0,
         maxcomp=6, algorithm=c("port","lbfgsb"),
         plot=TRUE, pcolor="blue", psize=1.5, 
         kratio=0.3, zscale=NULL) {
    UseMethod("RadialPlotter")
} #
### 2014.10.01.
RadialPlotter.default<-
function(EDdata, ncomp=0, addsigma=0,
         maxcomp=6, algorithm=c("port","lbfgsb"),
         plot=TRUE, pcolor="blue", psize=1.5, 
         kratio=0.3, zscale=NULL) {
    ### Stopifnot.
    stopifnot(ncol(EDdata)==2L, nrow(EDdata)>=5L,
              all(EDdata[,1L,drop=TRUE]>0),
              all(EDdata[,2L,drop=TRUE]>0),
              length(ncomp)==1L, is.numeric(ncomp),
              length(addsigma)==1L, addsigma>=0,
              length(maxcomp)==1L, maxcomp>=2L, maxcomp<=9L,
              ncomp %in% (-2L:maxcomp),
              is.character(algorithm),
              all(algorithm %in% c("lbfgsb","port")),
              length(plot)==1L, is.logical(plot),
              length(pcolor)==1L, is.character(pcolor),
              length(psize)==1L, is.numeric(psize),
              length(kratio)==1L, is.numeric(kratio),
              is.null(zscale) || is.numeric(zscale))
    ###
    ### R function for plot a radial plot.
    RadialPlot<-function(Data, Pars, addsigma, zscale,
                         kratio, pcolor, psize)  {
        z.i<-log(Data[,1L,drop=TRUE])
        se.i<-sqrt((Data[,2L,drop=TRUE]/
                    Data[,1L,drop=TRUE])^2L + addsigma^2L)
        ### If Pars=NULL, lines will not be plotted out.
        if(!is.null(Pars))  {
            Pars<-log(Pars)
        } # end if
        plot.area_ratio<-kratio
        z.0<-sum(z.i/se.i^2L)/sum(1.0/se.i^2L)
        x.i<-1.0/se.i
        y.i<-(z.i-z.0)/se.i
        h<-plot.area_ratio*(max(x.i)-min(x.i))/
                           (max(y.i)-min(y.i))
        r.0<-max(sqrt(x.i^2L+y.i^2L*h^2L))
        if(is.null(zscale))  {
            mkest<-round(exp(pretty(z.i)), 0L)
            if(sd(z.i)/mean(z.i)<0.05)  {
                mkest<-round(exp(pretty(c(min(z.i)*0.98, 
                       max(z.i)*1.02))), 0L)
            } # end if
        } else  {
            mkest<-zscale
        } # end if
        mkr<-range(mkest)
        circ.x<-(0:200)/200
        circ.x<-mkr[1L]+circ.x*(mkr[2L]-mkr[1L])
        circ.z<-log(circ.x)
        zmkest<-log(mkest) 
        circ.x<-r.0/sqrt(1.0+h^2L*(circ.z-z.0)^2L)
        circ.y<-(circ.z-z.0)*circ.x
        xmk1<-r.0/sqrt(1.0+h^2L*(zmkest-z.0)^2L)
        xmk2<-xmk1*1.01
        xmk3<-1.01*1.04*r.0/
              sqrt(1.0+h^2L*(zmkest-z.0)^2L)
        ymk1<-(zmkest-z.0)*xmk1
        ymk2<-(zmkest-z.0)*xmk2
        ymk3<-(zmkest-z.0)*xmk3 
        may<-max(abs(y.i))
        yaxis.max<-if(may>12)   {
            may  } else if(may<3) {
            8    } else  {
            12   } # end if
        par(mfrow=c(1,1), oma=c(2,0,0,0), 
            mar=c(5,4,4,4.5), xpd=TRUE, las=1, cex=1)
        plot(NA, NA, xlim=c(0,max(x.i)), ylim=c(-yaxis.max,yaxis.max), 
             xlab="", ylab="", xaxt="n", yaxt="n", main=NULL, 
             bty="n", typ="n", xaxs="i", yaxs="i")      
        lenpar<-length(Pars) 
        maxx<-which.max(circ.x)
        ### If Pars!=NULL, drawing lines.
        if(!is.null(Pars))  {
            for(i in 1L:lenpar)  {
                if(Pars[i]-z.0>=0.0)  {
                    Ci<-approx(circ.y[maxx:201L], circ.x[maxx:201L], 
                               ties="ordered", xout=(Pars[i]-z.0)*max(x.i), 
                               rule=2, method="constant", f=0.5)$y
                }  else  {
                    Ci<-approx(circ.y[1L:(maxx-1L)], circ.x[1L:(maxx-1L)], 
                               ties="ordered", xout=(Pars[i]-z.0)*max(x.i),
                               rule=2, method="constant",f=0.5)$y
                } # end if
               lines(c(0,Ci), c(0,(Pars[i]-z.0)*Ci), 
                     lty=1, col="black", lwd=1.5)
            } # end for
        } # end if
        lines(circ.x, circ.y, lwd=2)
        segments(xmk1, ymk1, xmk2, ymk2)
        text(xmk3, ymk3, round(mkest,2L), cex=1)
        text(max(x.i)*1.15,0,"De (Gy)", cex=1, srt=90)
        par(mfrow=c(1,1), oma=c(2,0,0,0.5), mar=c(5,4,4,4.5), 
            xpd=TRUE, las=1, cex=1, new=TRUE)
        plot(NA, NA, xlim=c(0,max(x.i)), ylim=c(-yaxis.max, yaxis.max), 
             yaxt="n", xaxt="n", xaxs="i", yaxs="i", 
             ylab="Standardised Estimate", 
             xlab="", type="n", bty="n")
        points(x.i, y.i, pch=21, col="black", 
               cex=psize, bg=pcolor)
        axis(side=1, line=4, cex.axis=1, lwd=2)
        mtext(side=1, line=6, "Precision", cex=1)
        reticks.labels<-round(1/axTicks(side=1)*100,1L)
        reticks.values<-axTicks(side=1)
        axis(side=1, at=reticks.values[-1L], lwd=2, 
             labels=reticks.labels[-1L], 
             line=4, cex.axis=1, tck=0.02, padj=-4)
        mtext(side=1, line=1.5,"Relative Error (%)", cex=1)
        axis(side=2,at=c(-2,-1,0,1,2), lwd=2, 
             labels=c("-2","","0","","2"), cex.axis=1)
        par(oma=c(0,0,0,0), xpd=FALSE,
            las=0, new=FALSE, mar=c(5,4,4,2)+0.1)
    } ### end function RadialPlot.
    ###
    ###
    ### R function for estimating standard errors of MAM.
    apMamStd<-function(ed,sed,pars) {
        ###
        ndat<-length(ed)
        np<-length(pars)
        stdp<-vector(length=np)
        iflag<-0
        ###
        res<-.Fortran("apmamstd",as.double(ed),as.double(sed),
                      as.integer(ndat),as.double(pars),
                      stdp=as.double(stdp),as.integer(np),
                      iflag=as.integer(iflag),PACKAGE="numOSL")
        if (res$iflag==0) {
            return(res$stdp)
        } else {
            stop("Error!")
        } # end if
        ###
    } # end function apMamStd.
    ###
    ###
    ### R function for MAM optimization.
    Rmam<-function(EDdata,ncomp,addsigma) {
        stopifnot(ncomp %in% c(-1L,-2L))
        x<-sqrt((EDdata[,2L,drop=TRUE]/
                 EDdata[,1L,drop=TRUE])^2L+addsigma^2L)
        y<-log(EDdata[,1L,drop=TRUE])
        ### 
        ### Function to be minimized.
        minfunc<-function(p)  {
            if(ncomp==-1L) {
                ### For MAM3.
                u0<-(p[2L]/(p[3L])^2L+y/x^2L)/
                    (1.0/(p[3L])^2L+1.0/x^2L)
                sigma0<-1.0/(sqrt(1.0/(p[3L])^2L+1.0/x^2L)) 
                prop1<-p[1L]/sqrt(2.0*pi)/x*
                       exp(-(y-p[2L])^2L/2.0/x^2L)
                prop2<-(1.0-p[1L])/sqrt(2.0*pi*((p[3L])^2L+x^2L))* 
                        2.0*(1.0-pnorm((p[2L]-u0)/sigma0))*
                        exp(-(y-p[2L])^2L/2.0/((p[3L])^2L+x^2L))
                return(-sum(log(prop1+prop2)))
            } else if(ncomp==-2L) {
                ### For MAM4.
                u0<-(p[3L]/(p[4L])^2L+y/x^2L)/
                    (1.0/(p[4L])^2L+1.0/x^2L)
                sigma0<-1.0/(sqrt(1.0/(p[4L])^2L+1.0/x^2L)) 
                prop1<-p[1L]/sqrt(2.0*pi)/x*exp(-(y-p[2L])^2L/2.0/x^2L)
                prop2<-(1.0-p[1L])/sqrt(2.0*pi*((p[4L])^2L+x^2L))* 
                       (1.0-pnorm((p[2L]-u0)/sigma0))/
                       (1.0-pnorm((p[2L]-p[3L])/p[4L]))*
                        exp(-(y-p[3L])^2L/2.0/((p[4L])^2L+x^2L))
                return(-sum(log(prop1+prop2)))
            } # end if
        } # end function minfunc.
        ###
        ### Set boundaries.
        if (ncomp==-1L) {
            lower<-c(1e-4, min(y), 1e-3)
            upper<-c(1-1e-4, max(y), 5.0)
        } else if (ncomp==-2L) {
            lower<-c(1e-4, min(y), min(y), 1e-3)
            upper<-c(1-1e-4, max(y), max(y), 5.0)
        } # end if
        ### 
        ### 
        kclus<-kmeans(x=y, centers=3L, 
                      iter.max=50L, nstart=100L)
        kclus<-sort(kclus$centers)
        ###
        ### Set gamma, mu and sdy.
        gama<-kclus[1L]
        mu<-c(kclus[2L], mean(kclus), mean(y))
        sdy<-sd(y)
        ### 
        cmaxlik<-maxlik<-1e20
        errorflag<-1
        bexist<-0
        ### Do optimization with various initial values.
        if(ncomp==-1L) {
            ### For MAM3.
            for(i in seq(3L)) {
                for(j in seq(5L)) {
                    for (k in seq(3L)) {
                        ###
                        ini<-c(0.01*(5.0)^(i-1L), mu[k], sdy*0.4*j)
                        ###
                        res<-if(algorithm[1L]=="lbfgsb") {
                            suppressWarnings(try(optim(par=ini, fn=minfunc, gr=NULL, 
                            method="L-BFGS-B", control=list(maxit=1000), lower=lower,
                            upper=upper, hessian=TRUE), silent=TRUE))
                        } else if(algorithm[1L]=="port") {
                            suppressWarnings(try(nlminb(start=ini, objective=minfunc, 
                            gradient=NULL, hessian=NULL, scale=1, control=list(iter.max=1000), 
                            lower=lower, upper=upper), silent=TRUE))
                        } # end if                     
                        ###
                        if(class(res)!="try-error" && 
                           res$convergence==0) { 
                            ###
                            min.obj<-ifelse(algorithm[1L]=="lbfgsb",
                                            res$value, res$objective)
                            ###
                            stdp<-if(algorithm[1L]=="lbfgsb") {
                                suppressWarnings(try(sqrt(diag(solve(res$hessian,
                                                     tol=1e-10))),silent=TRUE))
                            } else if(algorithm[1L]=="port") {
                                suppressWarnings(try(apMamStd(y,x,res$par),silent=TRUE))
                            } # end if    
                            ###
                            if(class(stdp)!="try-error" &&       
                               all(is.finite(stdp)) && 
                               min.obj<maxlik &&            
                               all(abs(res$par-lower)>=1e-5) &&   
                               all(abs(res$par-upper)>=1e-5) ) {  
                                pars<-res$par
                                error<-stdp
                                maxlik<-min.obj
                                errorflag<-0
                            } # end if
                            ###
                            if(class(stdp)!="try-error" &&  
                               all(is.finite(stdp)) &&               
                               min.obj<cmaxlik &&           
                               errorflag==1) {                    
                                cpars<-res$par
                                cerror<-stdp
                                cmaxlik<-min.obj
                                bexist<-1
                            } # end if
                            ###
                        } # end if
                        ###
                    } # end if
                    ###
                } # end for
                ###
            } # end for
            ###
        } else if(ncomp==-2L) {
            ### MAM4.
            for(i in seq(3L)) {
                for(j in seq(3L)) {
                    for(k in seq(3L)) {
                        ###
                        ini<-c(0.01*(5.0)^(i-1L), gama, mu[k], sdy*0.5*j)
                        ###
                        res<-if (algorithm[1L]=="lbfgsb") {
                            suppressWarnings(try(optim(par=ini, fn=minfunc, gr=NULL, 
                            method="L-BFGS-B", control=list(maxit=1000), lower=lower,
                            upper=upper, hessian=TRUE), silent=TRUE))
                        } else if (algorithm[1L]=="port") {
                            suppressWarnings(try(nlminb(start=ini, objective=minfunc, 
                            gradient=NULL, hessian=NULL, scale=1, control=list(iter.max=1000), 
                            lower=lower, upper=upper), silent=TRUE))
                        } # end if  
                        ###
                        if(class(res)!="try-error" && 
                           res$convergence==0) { 
                            ###
                            min.obj<-ifelse(algorithm[1L]=="lbfgsb",
                                            res$value, res$objective)
                            ###
                            stdp<-if(algorithm[1L]=="lbfgsb") {
                                suppressWarnings(try(sqrt(diag(solve(res$hessian))),silent=TRUE))
                            } else if(algorithm[1L]=="port") {
                                suppressWarnings(try(apMamStd(y,x,res$par),silent=TRUE))
                            } # end if  
                            ###
                            if(class(stdp)!="try-error" &&       
                               all(is.finite(stdp)) &&   
                               min.obj<maxlik &&           
                               all(abs(res$par-lower)>=1e-5) &&   
                               all(abs(res$par-upper)>=1e-5)) {      
                                pars<-res$par
                                error<-stdp
                                maxlik<-min.obj
                                errorflag<-0
                            } # end if
                            ###
                            if(class(stdp)!="try-error" && 
                               all(is.finite(stdp)) &&                
                               min.obj<cmaxlik &&           
                               errorflag==1) {  
                                cpars<-res$par
                                cerror<-stdp
                                cmaxlik<-min.obj
                                bexist<-1
                            } # end if
                            ###
                        } # end if
                        ###
                    } # end for
                    ###
                } # end for
                ###
            } # end for
            ###
        } # end if
        ###
        ### Output the results.
        if (errorflag==0) {
            ### Reset parameters before output.
            if (ncomp==-1L) {
                pars[2L]<-exp(pars[2L])
                error[2L]<-pars[2L]*error[2L]
            } else if (ncomp==-2L) {
                pars[2L:3L]<-exp(pars[2L:3L])
                error[2L:3L]<-pars[2L:3L]*error[2L:3L]
            } # end if
            return(list("pars"=cbind(pars,error),
                        "maxlik"=-maxlik,
                        "bic"=ifelse(ncomp==-1L,
                         2.0*maxlik+3L*log(length(y)),
                         2.0*maxlik+4L*log(length(y)))))
        } else if (bexist==1){
            ### Reset parameters before output.
            if (ncomp==-1L) {
                cpars[2L]<-exp(cpars[2L])
                cerror[2L]<-cpars[2L]*cerror[2L]
            } else if (ncomp==-2L) {
                cpars[2L:3L]<-exp(cpars[2L:3L])
                cerror[2L:3L]<-cpars[2L:3L]*cerror[2L:3L]
            } # end if
            return(list("pars"=cbind(cpars,cerror),
                        "maxlik"=-cmaxlik,
                        "bic"=ifelse(ncomp==-1L,
                         2.0*cmaxlik+3L*log(length(y)),
                         2.0*cmaxlik+4L*log(length(y)))))
        } else {
           return(NULL)
        } # end if
    } # end function Rmam.
    ###
    ###
    ndat<-nrow(EDdata)
    ed1<-as.numeric(EDdata[,1L,drop=TRUE])
    sed1<-as.numeric(EDdata[,2L,drop=TRUE])
    ###
    if (ncomp %in% c(-1L,-2L)) {
        res<-Rmam(EDdata=EDdata, ncomp=ncomp, addsigma=addsigma)
        if (is.null(res)) {
            if (plot==TRUE) {
             RadialPlot(Data=EDdata, Pars=NULL,
                        addsigma=addsigma, zscale=zscale,
                        kratio=kratio, pcolor=pcolor, psize=psize)
            } # end if
            stop("Error: fail in optimization!")
        } # end if
        ###
        maxlik<-res$maxlik
        bic<-res$bic
        ###
        ParsAndErrors<-res$pars
        colnames(ParsAndErrors)<-c("Pars", "Std.Pars")
        if (ncomp==-1L)  {
            rownames(ParsAndErrors)<-c("Proportion", "Minimum.ED", "Sigma")
        } else if (ncomp==-2L) {
            rownames(ParsAndErrors)<-c("Proportion", "Minimum.ED", 
                                       "Central.ED", "Sigma")
        } # end if
        ###
        if (plot==TRUE) {
             RadialPlot(Data=EDdata, Pars=ParsAndErrors[2L,1L],
                        addsigma=addsigma, zscale=zscale,
                        kratio=kratio, pcolor=pcolor, psize=psize)
        } # end if
    } else {
        if (ncomp==0L) {
            gcomp<--99
            message<-0
            res<-.Fortran("goodComp",as.double(ed1),as.double(sed1),
                          as.integer(ndat),as.integer(maxcomp),
                          gcomp=as.integer(gcomp),as.double(addsigma),
                          message=as.integer(message),PACKAGE="numOSL")
            if (res$message!=0) {
                if (plot==TRUE) {
                    RadialPlot(Data=EDdata, Pars=NULL,
                               addsigma=addsigma, zscale=zscale,
                               kratio=kratio, pcolor=pcolor, psize=psize)
                } # end if
                stop("Error: fail in optimization!")
            } # end if
            ncomp<-res$gcomp
        } # end if
        ###
        stdp<-pars<-matrix(0,nrow=2L,ncol=ncomp)
        maxlik<-bic<-0
        message<-0
        ###
        res<-.Fortran("compED",as.double(ed1),as.double(sed1),as.integer(ndat),
                      as.integer(ncomp),as.double(addsigma),pars=as.double(pars),
                      stdp=as.double(stdp),maxlik=as.double(maxlik),bic=as.double(bic),
                      message=as.double(message),PACKAGE="numOSL")
        if (res$message!=0) {
            if (plot==TRUE) {
                RadialPlot(Data=EDdata, Pars=NULL,
                           addsigma=addsigma, zscale=zscale,
                           kratio=kratio, pcolor=pcolor, psize=psize)
            } # end if
            stop("Error: fail in optimization!")
        } # end if
        ###
        maxlik<-res$maxlik
        bic<-res$bic
        ParsAndErrors<-cbind(matrix(res$pars, byrow=TRUE, ncol=2L),
                             matrix(res$stdp, byrow=TRUE, ncol=2L))
        if (ncomp==1L) {
            ParsAndErrors<-matrix(ParsAndErrors[,c(1L,3L,2L,4L)])
            rownames(ParsAndErrors)<-c("Overdispersion", "Std.Overdispersion", 
                                       "CAM.ED", "Std.CAM.ED")
            colnames(ParsAndErrors)<-"CAM"
            ###
            if (plot==TRUE) {
                RadialPlot(Data=EDdata, Pars=ParsAndErrors[3L],
                           addsigma=addsigma, zscale=zscale,
                           kratio=kratio, pcolor=pcolor, psize=psize)
            } # end if
        } else {
            ParsAndErrors<-ParsAndErrors[,c(1L,3L,2L,4L)][
            order(ParsAndErrors[,2L,drop=TRUE]),,drop=FALSE]
            colnames(ParsAndErrors)<-c("P","Std.P","ED","Std.ED") 
            rownames(ParsAndErrors)<-paste(rep("Comp", ncomp), 
                                     seq(ncomp), sep="") 
            if (plot==TRUE) {
                RadialPlot(Data=EDdata, Pars=ParsAndErrors[,3L,drop=TRUE],
                           addsigma=addsigma, zscale=zscale, kratio=kratio,
                           pcolor=pcolor, psize=psize)
            } # end if
        } # end if
    } # end if
    ###
    output<-list("pars"=round(ParsAndErrors,5L), 
                 "bic"=round(bic,5L), 
                 "maxlik"=round(maxlik,5L))
    ### 
    class(output)<-"RadialPlotter"
    ###
    return(output)
} # end function RadialPlotter.default.
