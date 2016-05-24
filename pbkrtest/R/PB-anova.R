

### ###########################################################
###
### PBanova
###
### ###########################################################


.PBanova <- function(largeModel, smallModel=NULL, nsim=200, cl=NULL){
    if (is.null(smallModel)){
        fixef.name <- rev(attr(terms(largeModel),"term.labels"))
        ans1 <- list()
        ans2 <- list()

        for (kk in seq_along(fixef.name)){
            dropped <- fixef.name[kk]
            newf <- as.formula(paste(".~.-", dropped))
            smallModel <- update(largeModel, newf)
            #cat(sprintf("dropped: %s\n", dropped))
            rr <- PBrefdist(largeModel, smallModel, nsim=nsim, cl=cl)
            ans1[[kk]] <- PBmodcomp(largeModel, smallModel, ref=rr)
            #ans2[[kk]] <- .FFmodcomp(largeModel, smallModel, ref=rr)
            largeModel <- smallModel
        }

        ans12 <- lapply(ans1, as.data.frame)
        ans22 <- lapply(ans2, as.data.frame)

        ans3 <- list()
        for (kk in seq_along(fixef.name)){
            dropped <- fixef.name[kk]
            ans3[[kk]] <-
		rbind(
                      cbind(term=dropped, test=rownames(ans12[[kk]]), ans12[[kk]],df2=NA),
                      cbind(term=dropped, test=rownames(ans22[[kk]][2,,drop=FALSE]), ans22[[kk]][2,,drop=FALSE]))

	}

        ans3 <- rev(ans3)
        ans4 <- do.call(rbind, ans3)
        rownames(ans4) <- NULL
        ans4$p<-round(ans4$p,options("digits")$digits)
        ans4$tobs<-round(ans4$tobs,options("digits")$digits)

        ans4
    } else {


    }
}




