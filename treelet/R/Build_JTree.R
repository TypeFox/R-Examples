Build_JTree <-
function(C, cc, maxlev,whichsave){
	myCs=list()
    dim_C = dim(C)[1]
    J = maxlev
    Z = matrix(rep(0, J * 2), ncol = 2)
    T = list()
    theta = rep(0, J)
    PCidx = matrix(rep(0, J * 2), ncol = 2)
    L = 1
    maskno = matrix()
    nodes = seq(1, dim_C, by = 1)
    dlabels = rep(0, dim_C)
    PC_ratio = rep(0, dim_C - 1)
    Zpos = matrix(rep(0, J * 2), ncol = 2)
    all_d = matrix(rep(0, J * dim_C), ncol = dim_C)
    all_nodes = matrix(rep(0, J * dim_C), ncol = dim_C)
    cc.out = rep(NA,maxlev)
    for (lev in 1:J) {
        mask_C = upper.tri(cc) * cc
        k = (mask_C == 0)
        mask_C[k] = -1
        mask_C[maskno, ] = -1
        mask_C[, maskno] = -1
        compno = which(mask_C == max(mask_C), arr.ind = TRUE)[1,]
        Cred = C[compno, compno]
	cc.out[lev] = cc[compno,compno][1,2]
        if (Cred[1, 2] == 0) {
            Cnew = C
            ccnew = cc
            R = diag(c(1, 1))
            theta = 0
            idx = c(1, 2)
        }
        else {
            C11 = Cred[1, 1]
            C22 = Cred[2, 2]
            C12 = Cred[1, 2]
            th = 1/2 * atan(2 * C12/(C11 - C22))
            cs = cos(th)
            sn = sin(th)
            R = rbind(c(cs, -sn), c(sn, cs))
            M = C
            M[compno, ] = t(R) %*% C[compno, ]
            C = M
            C[, compno] = M[, compno] %*% R
            Cred = C[compno, compno]
            idx = c(Cred[1, 1], Cred[2, 2])
            idx = sort.list(idx, decreasing = TRUE)
            dnew = diag(C)
            temp = sqrt(matrix(dnew[compno], ncol = 1) %*% dnew)
            temp = C[compno, ]/temp
            cc[compno, ] = temp
            cc[, compno] = t(temp)
        }
        PCidx[lev, ] = idx
        theta[lev] = th
        T[[lev]] = R
        Z[lev, ] = nodes[compno]
        pind = compno[idx]
        p1 = pind[1]
        p2 = pind[2]
        nodes[pind] = cbind(dim_C + lev, 0)
        dlabels[p2] = lev
        if (lev == 1) {
            maskno = p2
            maskno = as.matrix(maskno)
        }
        else {
            maskno = cbind(maskno, p2)
        }
        PC_ratio[lev] = C[p2, p2]/C[p1, p1]
        Zpos[lev, ] = compno
        all_d[lev, ] = t(dlabels)
        all_nodes[lev, ] = nodes
		if (lev %in% whichsave) {myCs[[lev]]=C}
		else myCs[[lev]]=NULL
if(lev%%100==0){print(lev);flush.console()}
    }
    return(list(Zpos = Zpos, T = T, PCidx = PCidx, all_nodes = all_nodes, 
        TreeCovs=myCs))

}
