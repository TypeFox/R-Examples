seeBarcode <- 
function(seq, col=c("green", "blue", "black", "red")){
                if(!is.null(dim(seq))) if(dim(seq)[1] > 1) stop("Single sequences only please!")
                bases <- c(136, 40, 72, 24)
                pos <- 1:length(seq)
                ind <- match(as.numeric(seq), bases)
                plot(1, 1, xlim=c(0,max(pos)), ylim=c(0,1), xaxt="n", yaxt="n", xlab=NA, ylab=NA,
                                bty="n", type="n", main=dimnames(seq)[[1]])
                abline(v=pos, col=col[ind])
}