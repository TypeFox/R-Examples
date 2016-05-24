Graf.Diagnostic<-function(x,k,m,Alpha,Beta,Color,LTY,xlab="x",ylab="Estimated DF's",
 main="Model Diagnostic",OneLegend=TRUE,lgnd1, lgnd2,arw1,arw2,...)
{
    #require(stepfun)
    #require(stats)
    if(is.vector(x)){
        n<-length(x)
        z<-sort(x)
        x<-matrix(x,n,1)
    }
    if(!is.matrix(x)) x<-as.matrix(x)
    if(is.matrix(x)){
        n<-nrow(x)
        z<-sort(x[,1])
    }
    
    p<-NULL
    eF<-NULL; sF<-NULL; eG<-NULL; sG<-NULL
    eF[1]<-0; sF[1]<-0; eG[1]<-0; sG[1]<-0
    p<-1/((m-k)*exp(Alpha+x%*%Beta)+(n-m+k))
    for(i in 1:n){
        if(m<n) xx<-c(x[1:k,1],x[(m+1):n,1])
        else xx<-x[1:k,1]
        eF[i+1]<-sum(xx<=z[i])/(n-m+k)
        eG[i+1]<-sum(x[(k+1):m,1]<=z[i])/(m-k)
        sF[i+1]<-sum(p*(x[,1]<=z[i]))
        sG[i+1]<-sum(exp(Alpha+x%*%Beta)*p*(x[,1]<=z[i]))
    }
    sfunc.EF<-stepfun(z,eF,f = 0, ties = "ordered")
    sfunc.SF<-stepfun(z,sF,f = 0, ties = "ordered")
    sfunc.EG<-stepfun(z,eG,f = 0, ties = "ordered")
    sfunc.SG<-stepfun(z,sG,f = 0, ties = "ordered")

    plot(sfunc.EF,do.points = FALSE,col.h=Color[1],col.v=Color[1],lty=LTY[1],xlab=xlab, ylab=ylab, main = NULL,...)
    plot(sfunc.SF,add=TRUE,do.points = FALSE,col.h=Color[2],col.v=Color[2],lty=LTY[2])
    plot(sfunc.EG,add=TRUE,do.points = FALSE,col.h=Color[3],col.v=Color[3],lty=LTY[3])
    plot(sfunc.SG,add=TRUE,do.points = FALSE,col.h=Color[4],col.v=Color[4],lty=LTY[4])
    #title(main=expression(paste("Plots of (",eF[k]," , ", sF[k]," ) and (",eG[k]," , ", sG[k]," )")))
    title(main=main)
    zp<-quantile(z,c(0, 0.85))
    if(OneLegend==TRUE){
        if(m==n) legend(lgnd1[1],lgnd1[2],lty=LTY,expression(hat(F)[hat(k)],tilde(F)[hat(k)],hat(F)[hat(k)],tilde(G)[hat(k)]),col=Color)
        else legend(lgnd1[1],lgnd1[2],lty=LTY,expression(hat(F)[hat(k)][hat(m)],tilde(F)[hat(k)][hat(m)],hat(F)[hat(k)][hat(m)],tilde(G)[hat(k)][hat(m)]),col=Color)
    }
    else{
        if(m==n){ 
            legend(lgnd1[1],lgnd1[2],lty=LTY[1:2],expression(paste(hat(F)[hat(k)],phantom(0)),tilde(F)[hat(k)]),col=Color[1:2])
            legend(lgnd2[1],lgnd2[2],lty=LTY[3:4],expression(paste(hat(G)[hat(k)],phantom(0)),tilde(G)[hat(k)]),col=Color[3:4])
        }
        else{
            legend(lgnd1[1],lgnd1[2],lty=LTY[1:2],expression(paste(hat(F)[hat(k)][hat(m)],phantom(0)),tilde(F)[hat(k)][hat(m)]),col=Color[1:2])
            legend(lgnd2[1],lgnd2[2],lty=LTY[3:4],expression(paste(hat(G)[hat(k)][hat(m)],phantom(0)),tilde(G)[hat(k)][hat(m)]),col=Color[3:4])
        }
    }
    arrows(arw1[1], arw1[2],arw1[3], arw1[4], length = 0.07, angle = 30, code = 2,
       col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    arrows(arw2[1], arw2[2], arw2[3],arw2[4], length = 0.07, angle = 30, code = 2,
      col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
}


Plot.ll<-function(x, ll, col, xaxis.lab = NULL, xlab="k", ylab="Loglikelihood", main="Plot of Loglikelihood",...)
{
    if(is.vector(x)){
        n<-length(x)
        z<-sort(x)
        x<-matrix(x,n,1)
    }
    if(!is.matrix(x)) x<-as.matrix(x)
    if(is.matrix(x)){
        n<-nrow(x)
        z<-sort(x[,1])
    }
    if(n<3) stop("too few observations")
    if(length(ll)!=(n+1)) cat("Error: Length of Log-likelihood is not the same as samplesize+1!\n")
    chpt<-order(ll)[n+1]-1
    cat("change-point = ",chpt,"\n")
    plot(0:n,ll,xlab=xlab,ylab=ylab, main=main,col=col, axes=FALSE, frame.plot=T,...)#,pch=".")
    segments(chpt,1.5*ll[n+1],chpt,ll[chpt+1],lty=2)
    axis(1, chpt-n*2.5/100,  tick = F, lab = expression(hat(k)), cex.axis = 1, col.axis = "red")
    axis(1, chpt,  tick = TRUE, "=", cex.axis = 1, col.axis = "black")
    sq<-3:10
    at.x<-seq(0, n, length = 1+min(sq[n%%sq==0]))
    at.x<-setdiff(at.x,(floor(chpt-n*6/100):ceiling(chpt+n*6/100)))
    x.lab<-NULL
    if(!is.null(xaxis.lab)){
        len<-length(xaxis.lab)
#        xaxis.lab<-as.vector(xaxis.lab)
        if(!is.character(xaxis.lab)){
            if(any(xaxis.lab!=sort(xaxis.lab))) stop("numeric xaxis.lab must be increasing")
            if(len>(n+1) || len<n) stop("wrong length of xaxis.lab")
            else{
                if(all((xaxis.lab[-1]-xaxis.lab[-len]) == (xaxis.lab[2]-xaxis.lab[1]))){
                    if(len==n){
                        x.lab[1]<-2*xaxis.lab[1]-xaxis.lab[2]
                        x.lab[2:(n+1)]<-xaxis.lab[1:n]
                    }
                    else x.lab<-xaxis.lab
                }
                else{
                    x.lab<-0:n
                    warning("xaxis.lab is ignored because it is not equal-paced sequence")
                }
            }
        }
        else{
            if(len>(n+1) || len<n) stop("wrong length of xaxis.lab")
            else{
                if(len==n){
                    x.lab[1]<-"start"
                    x.lab[2:(n+1)]<-xaxis.lab[1:n]
                }
                else x.lab<-xaxis.lab
            }
        }
        tempstr<-as.character(x.lab[chpt+1])
        axis(1, chpt+n*(nchar(tempstr)+.5)/100,  tick = F, tempstr, cex.axis = 1, col.axis = "red")
        axis(1, at = at.x,  tick = T, labels = formatC(x.lab[at.x+1], format = "fg"), cex.axis = 1, col.axis = "black")
    }
    else{ axis(1, at = at.x,  tick = T, labels = formatC(at.x, format = "fg"), cex.axis = 1, col.axis = "black")
        tempstr<-as.character(chpt)
        axis(1, chpt+n*(nchar(tempstr)+.5)/100,  tick = F, tempstr, cex.axis = 1, col.axis = "black")
    }
#    axis(1, chpt,  tick = TRUE, paste("k=",as.character(chpt),sep = "\0"), cex.axis = 1, col.axis = "black")
    axis(2)
}
