plot.gcrq <-
function(x, term, add=FALSE, y=FALSE, legend=FALSE, select.tau, deriv=FALSE, cv=FALSE, 
  transf=NULL, lambda0=FALSE, ...){ #, se=FALSE, intercept=FALSE, resid=TRUE, alpha=0.01, legend=TRUE, ...){
#x: un oggetto restituito da gcrq()
#add: se TRUE aggiunge le linee.....
#y: se TRUE e se l'oggetto x contiene y (i dati) allora li disegna. Se add=TRUE, y viene posto a FALSE
#... argomenti da passare a plot() e a matplot/matlines e text(). Quindi se c'e' col questo viene applicato sia a plot (se i dati devono essere disegnati) 
#   sia a matlines/matplot se lwd solo a matlines() se cex a plot e a text per legend.
#  alla legenda (se legend=TRUE) e sia a plot() 
#select.tau: which quantile curves should be drawn? default (missing) is all
#f.deriv: if TRUE, the first derivatives ofthe growth curves  are plotted
#cv: se TRUE disegna la cross validation versus lambdas
          if(is.null(x$info.smooth)) stop("plot for simple linear fits is not allowed")
          if(is.null(x$BB)) stop(" plot.gcrq() only works with smooth terms")          
          if(cv){
               if(is.null(x$cv)) stop("the object does not include the 'cv' component")
               valoriL<-x$cv[,1]
               valoriCV<-apply(x$cv[,-1],1,mean)
               if(!lambda0) {
                  valoriL<-valoriL[-1]
                  valoriCV<-valoriCV[-1]
                  }
               plot(valoriL, valoriCV ,type="o",ylab="Cross Validation score",
                xlab="lambda values",xaxt="n",...)
               axis(1, at=x$cv[,1], labels=round(x$cv[,1],2), las=2, cex.axis=.7)
               #boxplot(unlist(apply(m$cv[,-1],1,function(x)list(x)), recursive=FALSE))
               return(invisible(NULL))
              }
          if(length(x$BB)>1) {if(missing(term)) stop("please provide 'term'")} else {term=names(x$BB)}
          if(!term%in%names(x$BB)) stop("'term' is not a smooth variable")
          if(length(x$tau)<=1) select.tau<-1
          if(missing(select.tau)) {
              select.tau<-1:ncol(x$coefficients)
              } else {
              if(length(select.tau)>length(x$taus)) stop("`length(select.tau)<=length(taus)' is requested")
              if(all(select.tau<1 & select.tau>0)) select.tau<- match(select.tau,x$taus)
              if(!all(select.tau<=length(x$taus) & select.tau>=1)) stop("'select.tau' is not correctly specified")
              }
          if(deriv) y<-FALSE
          if(add) y<-FALSE
          BB<-x$BB[[term]]
          xvar.n<-attr(BB,"covariate.n")
          xvar.35<-attr(BB,"covariate.35")
          nomi.ok<-paste(term,"ps",1:ncol(BB),sep=".")
          b<-if(length(x$tau)<=1) x$coefficients[nomi.ok] else x$coefficients[nomi.ok,select.tau]
          fit.35<-if(deriv) x$Bderiv%*%b else BB%*%b #matrici
          #se <-. BB%*% V(b) %*%t(BB)
          if("(Intercept)"%in%rownames(as.matrix(x$coefficients))) {
                fit.35<-fit.35 + matrix(as.matrix(x$coefficients)["(Intercept)",], ncol=ncol(fit.35), nrow=nrow(fit.35), byrow=TRUE)
                }
          l<-c(list(x=xvar.35, y=fit.35),list(...))
          cexL<-if(is.null(l$cex)) .6 else l$cex #sara' usato solo se legend=TRUE
          if(!is.null(transf)) l$y <- eval(parse(text=transf), list(y=l$y))
          if(y && is.null(x$y)) warning("'y=TRUE' ignored.. the fit does not include the data", call.=FALSE)
          if(y && !is.null(x$y)) {
              l1<-c(list(x=xvar.n, y=x$y),list(...))
              if(!is.null(transf)) l1$y <- eval(parse(text=transf), list(y=l1$y))              
              if(is.null(l1$xlab)) l1$xlab<-term
              if(is.null(l1$ylab)) l1$ylab<-"Growth variable"
              if(!is.null(l1$col.p)) l1$col<-l1$col.p;l1$col.p<-NULL
              if(!is.null(l1$cex.p)) l1$cex<-l1$cex.p;l1$cex.p<-NULL
              if(!is.null(l1$pch.p)) l1$pch<-l1$pch.p;l1$pch.p<-NULL              
              if(legend) l1$xlim <- c(min(xvar.n),1.1*max(xvar.n))
              do.call(plot, l1)              
              #plot(xvar.n, x$y, xlab=term, ylab="Growth variable")
              if(legend) {
                  #cexL<-if(is.null(l1$cex)) .6 else l1$cex TOGLIERE
                  #text(1.05*max(xvar.n),  l$y[nrow(l$y),], x$taus[select.tau], cex=cexL)
                  text(1.05*max(xvar.n),  l$y[nrow(l$y),], formatC(x$taus[select.tau], digits=2, format="f"),
                      cex=cexL)
                  legend<-FALSE
                  }
              add<-TRUE
              }
          if(legend) {
            x.leg<-l$x[(length(l$x)-8)]
            l$x[(length(l$x)-11):(length(l$x)-5)]<-NA
          }
          if(add){
                do.call(matlines, l)
              } else {
                if(is.null(l$xlab)) l$xlab<-term
                if(is.null(l$ylab)) {l$ylab<-if(deriv) "Growth variable (first derivative)" else "Growth variable"}
                l$type<-"l"
                l$col.p<-NULL
                l$cex.p<-NULL
                l$pch.p<-NULL
                do.call(matplot, l)
                #if(y && !is.null(x$y)) points(xvar.n, x$y)
                #matplot(xvar.35, fit.35, type="l", xlab=term, ylab="", ...)
              }
          if(legend)  {
              #cexL<-if(is.null(l$cex)) .6 else l$cex TOGLIERE
              text(x.leg,  l$y[(length(l$x)-8),], formatC(x$taus[select.tau], digits=2, format="f"), cex=cexL)
            #tau<-x$tau; mtext(bquote(tau == .(tau)),line=-2)
            }
          }
