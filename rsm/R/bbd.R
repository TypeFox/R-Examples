### Generate a Box-Behnken design (up thru k=7)

bbd = function(k, n0=4, block = (k==4|k==5), randomize=TRUE, coding) 
{
    reftbl = list(NULL, NULL,
        list(c(1,2),c(1,3),c(2,3)),                        # k=3
        list(c(1,2),c(3,4), c(1,4),c(2,3), c(1,3),c(2,4)), # k=4
        list(c(1,2),c(1,3),c(3,4),c(4,5),c(2,5),
             c(1,4),c(1,5),c(2,3),c(2,4),c(3,5)),          # k=5
        list(c(1,2,4),c(2,3,5),c(3,4,6),
             c(1,4,5),c(2,5,6),c(1,3,6)),                  # k=6
        list(c(4,5,6),c(1,6,7),c(2,5,7),c(1,2,4),
             c(3,4,7),c(1,3,5),c(2,3,6))                   # k=7
    )
    
    CALL = match.call()
    yvars = NULL
    if (inherits(k, "formula")) {
        names = all.vars (k[[length(k)]])
        if (length(k) > 2) yvars = all.vars(k[[2]])
        k = length(names)
    }
    else
        names = paste("x", 1:k, sep="")
    
    if ((k<3) | (k>7))
        stop("Box-Behnken designs are available only for k=3:7")

    clist = reftbl[[k]]
    if (length(clist[[1]])==2)
        tbl = expand.grid(c(-1,1),c(-1,1))
    else
        tbl = expand.grid(c(-1,1),c(-1,1), c(-1,1))
    n = nrow(tbl)
    des = as.data.frame(matrix(0, nrow=n*length(clist), ncol=k))
    idx = 1:n - n
    for (i in 1:length(clist))
        des[idx + i*n, clist[[i]]] = tbl
    
    if (is.character(block)) {
        blkname = block
        block = TRUE
    }
    else
        blkname = "Block"
        
    blk = 0
    if (block) {
        if (k==4) 
            blk = c(rep(1:3, rep(2*n, 3)), rep(1:3, n0))
        else if (k==5) 
            blk = c(rep(1:2, rep(5*n, 2)), rep(1:2, n0))
        else
            stop("Can only block when k=4 or k=5")
        nblk = ifelse(k==4, 3, 2)
    }
    else 
        nblk = 1
    des = rbind(des, matrix(0, nrow=n0*nblk, ncol=k))
    names(des) = names
    if (block) {
        des = cbind(factor(blk), des)
        names(des)[1] = blkname
        des = des[order(blk), ]
    }
    row.names(des) = 1:nrow(des)
    
    if (!is.null(yvars))
        for (v in yvars)  des[[v]] = NA
    
    if (missing(coding))
        coding = sapply(names, function(v) as.formula(paste(v,"~",v,".as.is", sep="")))
    des = as.coded.data (des, formulas=coding)

    # create design info as if each block is a CCD in  one block
    N = nrow(des) / nblk
    rsd = list(
        primary = names,
        call = CALL
#         n0 = c(n0,0),
#         non0 = c(N-n0,0),
#         alpha = 1
    )
     if (block) {
#         rsd$blk.info = rep(list(rsd), nblk)
         rsd$block = blkname
     }
#     rsd$call = CALL
    attr(des, "rsdes") = rsd

    des = .randomize(des, randomize=randomize)
    
    des
}
