"multilap" <-
function(df=DAAG::nsw74psid1, maxf=20, colnames=c("educ", "age", "re74", "re75", "re78")){
    if(maxf==Inf) return(rep(TRUE, dim(df)[1]))
    if (length(maxf)==1) maxf <- c(1/maxf, maxf)
    trt <- df$trt
    common <- rep(TRUE, length(trt))
    for (vname in colnames){
        y0 <- df[trt==0, vname]
        y1 <- df[trt==1, vname]
        xchop <- overlapDensity(y0, y1, ratio=maxf,
                           ratio.number=TRUE, plotvalues="")
        common <- common & df[,vname]>=xchop[1] & df[,vname] <= xchop[2]
    }
    invisible(common)
}

