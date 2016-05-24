NRadjClass <- function(ID, NRclass, resp, preds=NULL, wts=NULL, type){
    if (!(type %in% 1:5))
        stop("type must be in 1:5.\n")

    if (type == 1){
        xx <- by(data = preds, INDICES = NRclass, FUN = mean)
    }
    if (type == 2){
        xx <- by(data = data.frame(preds = preds, wt = wts), INDICES = NRclass,
           function(x) {weighted.mean(x$preds, x$wt)})
    }
    if (type == 3){
        xx <- by(data = resp, INDICES = NRclass, FUN = mean)
    }
    if (type == 4){
        xx <- by(data = data.frame(resp = resp, wt = wts), INDICES = NRclass,
           function(x) {weighted.mean(x$resp, x$wt)})
    }
    if (type == 5){
        xx <- by(data = preds, INDICES = NRclass, FUN = median)
    }
    xx <- data.frame(cbind(NRcl.no = 1:length(xx), RR = xx))
    mdata <- data.frame(cbind(ID=ID, NRcl.no = as.numeric(NRclass), resp=resp))
    tmp <- merge(mdata, xx, by="NRcl.no")
    tmp <- tmp[tmp$resp == 1,]
    data.frame(tmp)
}
