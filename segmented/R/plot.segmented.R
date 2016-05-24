plot.segmented<-function (x, term, add = FALSE, res = FALSE, conf.level = 0, 
    interc=TRUE, link = TRUE, res.col = 1, rev.sgn = FALSE, const = 0, 
    shade=FALSE, rug=TRUE, dens.rug=FALSE, dens.col = grey(0.8),
    show.gap=FALSE, transf=I, ...){
#funzione plot.segmented che consente di disegnare anche i pointwise CI
        f.U<-function(nomiU, term=NULL){
        #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
        #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
            k<-length(nomiU)
            nomiUsenzaU<-strsplit(nomiU, "\\.")
            nomiU.ok<-vector(length=k)
            for(i in 1:k){
                nomi.i<-nomiUsenzaU[[i]][-1]
                if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".")
                nomiU.ok[i]<-nomi.i
                }
          if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
          return(nomiU.ok)
        }
#-------------- 
        enl.range<-function(..., enlarge=TRUE){
        #modifica il min dei valori in ...
          r<-range(...)
          if(enlarge) r[1]<-if(sign(r[1])>0) r[1]*.9 else r[1]*1.1
          r
        }
#--------------
  #se l'oggetto e' segmented.Arima il nome dell'eventuale interc va sostituito..
  if((all(class(x)==c("segmented", "Arima")))) names(x$coef)<-gsub("intercept", "(Intercept)", names(coef(x)))
#--------------
    linkinv <- !link
    if (inherits(x, what = "glm", which = FALSE) && linkinv && !is.null(x$offset) && res) stop("residuals with offset on the response scale?")
    if(conf.level< 0 || conf.level>.9999) stop("meaningless 'conf.level'")
    if ((inherits(x, what = "glm", which = FALSE) && linkinv) || res) {
        if(!(identical(transf, I) || identical(transf, "I"))) {transf<-I; warning("'transf' set to I..")}
        }
    show.gap<-FALSE
    if (missing(term)) {
        if (length(x$nameUV$Z) > 1) {
            stop("please, specify `term'")
        }
        else {
            term <- x$nameUV$Z
        }
    }
    else {
        if (!term %in% x$nameUV$Z)
            stop("invalid `term'")
    }
    opz <- list(...)
    cols <- opz$col
    if (length(cols) <= 0)
        cols <- 1
    lwds <- opz$lwd
    if (length(lwds) <= 0)
        lwds <- 1
    ltys <- opz$lty
    if (length(ltys) <= 0)
        ltys <- 1
    cexs <- opz$cex
    if (length(cexs) <= 0)
        cexs <- 1
    pchs <- opz$pch
    if (length(pchs) <= 0)
        pchs <- 1
    xlabs <- opz$xlab
    if (length(xlabs) <= 0)
        xlabs <- term
    ylabs <- opz$ylab
    if (length(ylabs) <= 0)
        ylabs <- paste("Effect  of ", term, sep = " ")
    a <- intercept(x, term, gap = show.gap)[[1]][, "Est."]
    #Poiche' intercept() restituisce quantita' che includono sempre l'intercetta del modello, questa va eliminata se interc=FALSE
    if(!interc && ("(Intercept)" %in% names(coef(x)))) a<- a-coef(x)["(Intercept)"]
    b <- slope(x, term)[[1]][, "Est."]
    
    #id <- grep(paste("\\.", term, "$", sep = ""), rownames(x$psi), value = FALSE) #confondeva "psi1.x","psi1.neg.x"
    id <- f.U(rownames(x$psi), term)
    
    est.psi <- x$psi[id, "Est."]
    K <- length(est.psi)
    val <- sort(c(est.psi, x$rangeZ[, term]))
    #---------aggiunta per gli IC
    rangeCI<-NULL
    n<-length(x$residuals) #fitted.values - Arima non ha "fitted.values", ma ha "residuals"..
    tipo<- if(inherits(x, what = "glm", which = FALSE) && link) "link" else "response"
    
    vall<-sort(c(seq(min(val), max(val), l=120), est.psi))
    #ciValues<-predict.segmented(x, newdata=vall, se.fit=TRUE, type=tipo, level=conf.level)
    vall.list<-list(vall)
    names(vall.list)<-term
    ciValues<-broken.line(x, vall.list, link=link, interc=interc, se.fit=TRUE)
    
    if(conf.level>0) {
        k.alpha<-if(inherits(x, what = c("glm","Arima"), which = FALSE)) abs(qnorm((1-conf.level)/2)) else abs(qt((1-conf.level)/2, x$df.residual))
        ciValues<-cbind(ciValues$fit, ciValues$fit- k.alpha*ciValues$se.fit, ciValues$fit + k.alpha*ciValues$se.fit)
        #---> transf...
        ciValues<-apply(ciValues, 2, transf)
        rangeCI<-range(ciValues)
        #ciValues  e' una matrice di length(val)x3. Le 3 colonne: stime, inf, sup
        #polygon(c(vall, rev(vall)), c(ciValues[,2],rev(ciValues[,3])), col = "gray", border=NA)
        }

    #---------
    a.ok <- c(a[1], a)
    b.ok <- c(b[1], b)
    y.val <- a.ok + b.ok * val + const
    a.ok1 <- c(a, a[length(a)])
    b.ok1 <- c(b, b[length(b)])
    y.val <- y.val1 <- a.ok1 + b.ok1 * val + const
    s <- 1:(length(val) - 1)
    xvalues <- if(all(class(x)==c("segmented", "Arima"))) x$Z[,1] else  x$model[, term]
    if (rev.sgn) {
        val <- -val
        xvalues <- -xvalues
    }
    m <- cbind(val[s], y.val1[s], val[s + 1], y.val[s + 1])
    #values where to compute predictions (useful only if res=TRUE)
    if(res){
        new.d<-data.frame(ifelse(rep(rev.sgn, length(xvalues)),-xvalues, xvalues))
        names(new.d)<-term
        fit0 <- broken.line(x, new.d, link = link, interc=interc, se.fit=FALSE)$fit
        }
#-------------------------------------------------------------------------------
    if (inherits(x, what = "glm", which = FALSE) && linkinv) { #se GLM con linkinv
        fit <- if (res)
            #predict.segmented(x, ifelse(rep(rev.sgn, length(xvalues)),-xvalues,xvalues), type=tipo) + resid(x, "response") + const
            #broken.line(x, term, gap = show.gap, link = link) + resid(x, "response") + const
              fit0 + resid(x, "response") + const        
                else x$family$linkinv(c(y.val, y.val1))
        xout <- sort(c(seq(val[1], val[length(val)], l = 120), val[-c(1, length(val))]))
        l <- approx(as.vector(m[, c(1, 3)]), as.vector(m[, c(2, 4)]), xout = xout)
        id.group <- cut(l$x, val, FALSE, TRUE)
        yhat <- l$y
        xhat <- l$x
        m[, c(2, 4)] <- x$family$linkinv(m[, c(2, 4)])
        if (!add) {
            plot(as.vector(m[, c(1, 3)]), as.vector(m[, c(2,
                4)]), type = "n", xlab = xlabs, ylab = ylabs,
                main = opz$main, sub = opz$sub, 
                xlim = opz$xlim,
                ylim = if(is.null(opz$ylim)) enl.range(fit, rangeCI, enlarge=dens.rug) else opz$ylim )
        if(dens.rug){
          density <- density( xvalues )
          # the height of the densityity curve
          max.density <- max(density$y)
          # Get the boundaries of the plot to
          # put the density polygon at the x-line
          plot_coordinates <- par("usr")
          # get the "length" and range of the y-axis
          y.scale <- plot_coordinates[4] - plot_coordinates[3]
          # transform the y-coordinates of the density
          # to the lower 10% of the plotting panel
          density$y <- (0.1 * y.scale / max.density) * density$y + plot_coordinates[3]
          ## plot the polygon
          polygon( density$x , density$y , border = F , col = dens.col) 
          box()
          }
          
        if(rug) {
            segments(xvalues, rep(par()$usr[3],length(xvalues)), xvalues,
              rep(par()$usr[3],length(xvalues))+ abs(diff(par()$usr[3:4]))/40)}
            }
       
        if(conf.level>0){
          if(rev.sgn) vall<- -vall
          if(shade) polygon(c(vall, rev(vall)), c(ciValues[,2],rev(ciValues[,3])),
            col = "gray", border=NA) else matlines(vall, ciValues[,-1], type="l", lty=2, col=cols)
            }
        if (res) points(xvalues, fit, cex = cexs, pch = pchs, col = res.col)
        yhat <- x$family$linkinv(yhat)
        if (length(cols) == 1)
            cols <- rep(cols, max(id.group))
        if (length(lwds) == 1)
            lwds <- rep(lwds, max(id.group))
        if (length(ltys) == 1)
            ltys <- rep(ltys, max(id.group))
        for (i in 1:max(id.group)) {
            lines(xhat[id.group == i], yhat[id.group == i], col = cols[i],
                lwd = lwds[i], lty = ltys[i])
        }
#-------------------------------------------------------------------------------
    } else { #se LM o "GLM con link=TRUE (ovvero linkinv=FALSE)"
        ##---> transf!!!
        y.val<- do.call(transf, list(y.val)) 
        y.val1<-do.call(transf, list(y.val1))
        r <- cbind(val, y.val)
        r1 <- cbind(val, y.val1)
        rr <- rbind(r, r1)
        fit <- c(y.val, y.val1)
        if (res) {
            ress <- if (inherits(x, what = "glm", which = FALSE))
                residuals(x, "working") * sqrt(x$weights)
            else resid(x)
            #if(!is.null(x$offset)) ress<- ress - x$offset
            #fit <- broken.line(x, term, gap = show.gap, link = link, interc = TRUE) + ress + const
            #fit <- predict.segmented(x, ifelse(rep(rev.sgn, length(xvalues)),-xvalues,xvalues), type=tipo) + ress + const
            fit <- fit0 + ress + const
        }
        if (!add)
            plot(rr, type = "n", xlab = xlabs, ylab = ylabs,
                main = opz$main, sub = opz$sub, 
                xlim = opz$xlim,
                #ylim = if(is.null(opz$ylim)) enl.range(fit, rangeCI, enlarge=dens.rug) else opz$ylim)
                ylim = if(is.null(opz$ylim)) enl.range(fit, rangeCI, do.call(transf, list(m[, c(2,4)])), enlarge=dens.rug) else opz$ylim)
        if(dens.rug){
          density <- density( xvalues )
          # the height of the densityity curve
          max.density <- max(density$y)
          # Get the boundaries of the plot to
          # put the density polygon at the x-line
          plot_coordinates <- par("usr")
          # get the "length" and range of the y-axis
          y.scale <- plot_coordinates[4] - plot_coordinates[3]
          # transform the y-coordinates of the density
          # to the lower 10% of the plotting panel
          density$y <- (0.1 * y.scale / max.density) * density$y + plot_coordinates[3]
          ## plot the polygon
          polygon( density$x , density$y , border = F , col = dens.col) 
          box()
          }
        if(rug) {segments(xvalues, rep(par()$usr[3],length(xvalues)), xvalues,
            rep(par()$usr[3],length(xvalues))+ abs(diff(par()$usr[3:4]))/40)}


        if(conf.level>0) {
          if(rev.sgn) vall<- -vall
          if(shade) polygon(c(vall, rev(vall)), c(ciValues[,2],rev(ciValues[,3])),
            col = "gray", border=NA) else matlines(vall, ciValues[,-1], type="l", lty=2, col=cols)
            }
        if (res)
            points(xvalues, fit, cex = cexs, pch = pchs, col = res.col)
            segments(m[, 1], do.call(transf, list(m[, 2])), m[, 3], do.call(transf, list(m[, 4])), 
                col = cols, lwd = lwds, lty = ltys)
        }
    invisible(NULL)
}
