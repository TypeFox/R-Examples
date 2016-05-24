'Frechet.bounds.cat' <-
function(tab.x, tab.xy, tab.xz, print.f="tables", tol= 0.001)
{
    
############# start fb.yz function
    fb.yz <-  function(y, z, prn="tables"){
        lab.y <- names(y)
        if(is.null(lab.y)) lab.y <- paste("y",1:length(y), sep="")
        lab.z <- names(z)
        if(is.null(lab.z)) lab.z <- paste("z",1:length(z), sep="")
  
        p.y <- prop.table(y)
        p.z <- prop.table(z)
  
        ll <- outer(p.y, p.z, FUN="+") - 1
        m0 <- matrix(0, nrow(ll), ncol(ll))
        low <- pmax(m0, ll)
        upper <- outer(p.y, p.z, FUN="pmin")
  
        ind <-  outer(p.y, p.z, FUN="*")
        
        dimnames(low) <- dimnames(upper) <- dimnames(ind) <- list(lab.y, lab.z)
        class(low) <- class(upper) <- class(ind)  <- "table"
        
        H.y <- sum(-1*p.y*log(p.y), na.rm=TRUE)
        H.z <- sum(-1*p.z*log(p.z), na.rm=TRUE)
        
        res.0 <- list(low.u=low, up.u=upper, IA=ind, H=c(H.y, H.z), uncertainty=mean(upper-low))
  
        if(prn=="tables"){
            out <- res.0
        }
        else if(prn=="data.frame"){
            df <- data.frame(low)
            colnames(df) <- c("Y", "Z", "low.u")
            df$IA <- c(ind)
            df$up.u <- c(upper)
            out <- list(bounds=df, H=c(H.y, H.z), uncertainty=mean(upper-low))
        }  
        out
    }

### end  function fb.yz
########################

    if(is.null(tab.x)) out <- fb.yz(y=tab.xy, z=tab.xz, prn=print.f)    
    else{
        lab.x <- names(dimnames(tab.x))
        if(all(nchar(lab.x)==0)) lab.x <- paste("x",1:length(lab.x), sep="")
        names(attr(tab.x, "dimnames")) <- lab.x

        lab.xy <- names(dimnames(tab.xy))
        if(all(nchar(lab.xy)==0)) lab.xy <- c(lab.x, "y")
        names(attr(tab.xy, "dimnames")) <- lab.xy

        lab.y <- setdiff(lab.xy, lab.x)
        pos.y <- match(lab.y, lab.xy)
       
        lab.xz <- names(dimnames(tab.xz))
        if(all(nchar(lab.xz)==0)) lab.xz <- c(lab.x, "z")
        names(attr(tab.xz, "dimnames")) <- lab.xz

        lab.z <- setdiff(lab.xz, lab.x)
        pos.z <- match(lab.z, lab.xz)

        p.x <- prop.table(tab.x)
        p.xy <- prop.table(tab.xy)
        p.y <- margin.table(p.xy, pos.y)
        
        p.xz <- prop.table(tab.xz)
        p.z <- margin.table(p.xz, pos.z)

# check marginal distribution of the X variables

        d1.x <- 1:length(dim(p.xy))
        m1.x <- margin.table(p.xy, d1.x[-pos.y])
        if(any(abs(m1.x-p.x)>tol) )
            warning("The marginal distr. of the X variables \n in tab.xy is not equal to tab.x")

        d2.x <- 1:length(dim(p.xz))
        m2.x <- margin.table(p.xz, d2.x[-pos.z])
        if(any(abs(m2.x-p.x)>tol) )
            warning("The marginal distr. of the X variables \n in tab.xz is not equal to tab.x")

        if(any(abs(m1.x-m2.x)>tol) )
            warning("The marginal distr. of the X variables \n in tab.xy and in tab.xz are not equal")

########################################################
# computes Frechet bounds _without_ using X variables

        ll <- outer(p.y, p.z, FUN="+") - 1
        low <- pmax(matrix(0, nrow(ll), ncol(ll)), ll)
        upper <- outer(p.y, p.z, FUN="pmin")

        dimnames(low) <- dimnames(upper) <- list(names(p.y), names(p.z))
        class(low) <- class(upper)  <- "table"
        res.0 <- list(low.u=low, up.u=upper)

#############################################
# computes Frechet bounds using X variables

        if(length(dim(p.x))==1) xx.0 <- cbind(p.x)
        else xx.0 <- ftable(p.x, row.vars=1:length(dim(p.x)))
        xx.y <- ftable(p.xy, col.vars = length(dim(p.xy)))
        ygxx <- prop.table(xx.y, margin=1)
        xx.z <- ftable(p.xz, col.vars = length(dim(p.xz)))
        zgxx <- prop.table(xx.z, margin=1)
        
        H <- nrow(xx.0)
        out.CIA <- out.low <- out.up <- array(0, dim=c(ncol(xx.y), ncol(xx.z), H))
#       out.D <- array(0, dim=c(ncol(xx.y), ncol(xx.z), H))
        for(h in 1:H){
            out.CIA[ , ,h] <- outer(ygxx[h,], zgxx[h,], FUN="*") * xx.0[h,]
            thetas <- outer(ygxx[h,], zgxx[h,], FUN="+")-1
            ll <- pmax(thetas, matrix(0, nrow=nrow(thetas), ncol=ncol(thetas)))
            uu <- outer(ygxx[h, ], zgxx[h, ], FUN="pmin")
            out.low[ , , h] <- ll*xx.0[h, ]
            out.up[  , , h] <- uu*xx.0[h, ]
#            jgi <- outer(ygxx[h, ], rep(1,length(zgxx[h, ])), FUN="*")
#            kgi <- outer(rep(1,length(ygxx[h, ])), zgxx[h, ], FUN="*")
#            out.D[ , ,h] <- (uu - ll)*jgi*kgi*xx.0[h, ]
        }

        fine.CIA <- apply(out.CIA, c(1,2), sum, na.rm=TRUE)
        fine.low <- apply(out.low, c(1,2), sum, na.rm=TRUE)
        fine.up <- apply(out.up, c(1,2), sum, na.rm=TRUE)

        l.y <- attr(ygxx, "col.vars")[[1]]
        l.z <- attr(zgxx, "col.vars")[[1]]

        class(fine.CIA) <- class(fine.low) <- class(fine.up) <- "table"
        dimnames(fine.CIA) <- dimnames(fine.low) <-  dimnames(fine.up) <- list(l.y, l.z)

# uncertainty        
#        H.x <- sum(-1 * p.x * log(p.x), na.rm=TRUE)
#        H.y <- sum(-1 * p.y * log(p.y), na.rm=TRUE)
#        H.z <- sum(-1 * p.z * log(p.z), na.rm=TRUE)
#        H.ygx <- sum(-1 * xx.y * log(ygxx), na.rm=TRUE)
#        H.zgx <- sum(-1 * xx.z * log(zgxx), na.rm=TRUE)
#        vet.H <- c(H.x=H.x, H.y=H.y, H.ygx=H.ygx, H.z=H.z, H.zgx=H.zgx)
#        vet.U <- c(U.ygx = 1-H.ygx/H.y, U.zgx = 1-H.zgx/H.z)
#
#        vet.unc <- c(av.u=mean(c(upper-low)), av.cx=mean(c(fine.up-fine.low)), delta.CMS=sum(c(out.D), na.rm=TRUE))
        vet.unc <- c(av.u=mean(c(upper-low)), av.cx=mean(c(fine.up-fine.low)))
        res.1 <- list(CIA=fine.CIA, low.cx=fine.low, up.cx=fine.up) 
#        res.2 <- list(H=vet.H, U=vet.U, uncertainty=vet.unc)
        res.2 <- list(uncertainty=vet.unc)
        
        if(print.f=="tables"){
            out <- c(res.0, res.1, res.2)
        }
        else if(print.f=="data.frame"){
            dataf <- data.frame(res.0$low.u)
            labdf <- c(lab.y, lab.z, "low.u")
            colnames(dataf) <- labdf
            dataf$low.cx <- c(res.1$low.cx)
            dataf$CIA <- c(res.1$CIA)
            dataf$up.cx <- c(res.1$up.cx)
            dataf$up.u <- c(res.0$up.u)
            out <- c(list(bounds=dataf), res.2)
        }
    }
    out
}
