varsup <- function(resmca,var) {
        n <- sum(resmca$call$row.w)
        #n <- length(resmca$call$row.w)
        FK <- colSums(resmca$call$row.w*(dichotom(as.data.frame(factor(var)),out='numeric')))/n
        classe <- class(resmca)[1] # new
        if(classe=='stMCA') classe=resmca$call$input.mca # new
        if(classe %in% c('MCA','speMCA','multiMCA')) { # new
           v <- factor(var)
           wt <- resmca$call$row.w
           ind <- resmca$ind$coord
           coord <- aggregate(wt*ind,list(v),sum)[,-1]/n/FK
           vrc <- aggregate(wt*ind*ind,list(v),sum)[,-1]/n/FK-coord*coord
           for(i in 1:resmca$call$ncp) coord[,i] <- coord[,i]/resmca$svd$vs[i]
           cos2 <- coord*coord/((1/FK)-1)
           weight=n*FK
           }
        if(classe == 'csMCA') { # new
           v <- factor(var[resmca$call$subcloud])
           wt <- resmca$call$row.w[resmca$call$subcloud]
           #ind <- resmca$ind$coord[resmca$call$var %in% resmca$call$lev,]
           ind <- resmca$ind$coord
           fK <- colSums(wt*(dichotom(as.data.frame(v),out='numeric')))/length(wt)
           n.w  <- sum(wt)
           #coord <- aggregate(wt*ind,list(v),sum)[-1]/n.w/FK
           #vrc <- aggregate(wt*ind*ind,list(v),sum)[,-1]/n.w/FK-coord*coord
           coord <- aggregate(wt*ind,list(v),sum)[-1]/n.w/fK
           vrc <- aggregate(wt*ind*ind,list(v),sum)[,-1]/n.w/fK-coord*coord
           for(i in 1:resmca$call$ncp) coord[,i] <- coord[,i]/resmca$svd$vs[i]
           cos2 <- coord*coord*FK*FK/fK/(1-fK)
           weight <- length(wt)*fK
           }
        names(weight) <- levels(v)# v au lieu de var
        rownames(coord) <- levels(v)#[as.numeric(table(v))>0]
        rownames(cos2) <- levels(v)#[as.numeric(table(v))>0]
        wi <- apply(vrc,2,weighted.mean,w=weight)
        be <- resmca$eig[[1]][1:resmca$call$ncp]-wi
        eta2 <- be/resmca$eig[[1]][1:resmca$call$ncp]
	vrc <- rbind(vrc,wi,be,resmca$eig[[1]][1:resmca$call$ncp],eta2)
        vrc <- round(vrc,6)
	rownames(vrc) <- c(levels(v),'within','between','total','eta2')
        coord <- round(coord,6)
        v.test <- sqrt(cos2)*sqrt(length(v)-1)
        v.test <- (((abs(coord)+coord)/coord)-1)*v.test
        list(weight=round(weight,1),coord=coord,cos2=round(cos2,6),var=round(vrc,6),v.test=round(v.test,6))
        }
