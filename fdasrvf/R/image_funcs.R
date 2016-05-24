apply_gam_to_gam <- function(gamnew, gam){
    # gam \circ gam0
    m = dim(gam)[1]; n = dim(gam)[2]; D = dim(gam)[3]
    md = 8
    mt = md*m
    nt = md*n
    U = seq(0, 1, length.out=m)
    V = seq(0, 1, length.out=n)
    xlim = c(0,1);
    ylim = c(0,1);
    dx1 = (1-0)/(m-1)
    dy1 = (1-0)/(n-1)
    Ut = seq(0, 1, length.out=mt)
    Vt = seq(0, 1, length.out=nt)
    dx = (1-0)/(mt-1)
    dy = (1-0)/(nt-1)
    gam_tmp = array(0, dim=c(mt,nt,D))
    gam_new_tmp = array(0, dim=c(mt,nt,D))
    for (i in 1:D) {
        if (requireNamespace("akima", quietly = TRUE)) {
            gam_tmp[,,i] = akima::bicubic.grid(U,V,gam[,,i],xlim,ylim,dx,dy)$z
            gam_new_tmp[,,i] = akima::bicubic.grid(U,V,gamnew[,,i],xlim,ylim,dx,dy)$z
        } else {
            grid.list<- list(x=Ut, y=Vt)
            obj<-list(x=U, y=V, z=gam[,,i])
            gam_tmp[,,i] = interp.surface.grid(obj, grid.list)$z
            obj<-list(x=U, y=V, z=gamnew[,,i])
            gam_new_tmp[,,i] = interp.surface.grid(obj, grid.list)$z
        }
        
    }

    gam_cum_tmp = array(0,dim=c(mt,nt,D))
    x2 = c(gam_tmp[,,2])
    y2 = sort(c(t(gam_tmp[,,1])))
    for (i in 1:D) {
        if (requireNamespace("akima", quietly = TRUE)) {
            tmp = akima::bicubic(Ut,Vt,gam_new_tmp[,,i],x2,y2)$z
        } else {
            grid.list<- cbind(x2, y2)
            obj<-list(x=Ut, y=Vt, z=gam_new_tmp[,,i])
            tmp = interp.surface(obj, grid.list)
        }
        gam_cum_tmp[,,i] = matrix(tmp,nrow=mt,byrow=F)
    }

    gam_cum = array(0,dim=c(m,n,D))
    for (i in 1:D) {
        if (requireNamespace("akima", quietly = TRUE)) {
            gam_cum[,,i] = akima::bicubic.grid(Ut,Vt,gam_cum_tmp[,,i],xlim,ylim,dx1,dy1)$z
        } else {
            grid.list<- list(x=U, y=V)
            obj<-list(x=Ut, y=Vt, z=gam_cum_tmp[,,i])
            gam_tmp[,,i] = interp.surface.grid(obj, grid.list)$z
        }
    }

    return(gam_cum)

}


apply_gam_to_imag <- function(img, gam){
    if (ndims(img)==3){
        m = dim(img)[1]; n = dim(img)[2]; d = dim(img)[3]
        img_new = array(0,dim(img))
    } else if (ndims(img) == 2){
        m = dim(img)[1]; n = dim(img)[2]
        d = 1;
        img_new = matrix(0,m,n)
    }

    U = seq(0,1,length.out=m)
    V = seq(0,1,length.out=n)
    x2 = c(gam[,,2])
    y2 = sort(c(t(gam[,,1])))
    if (d==1){
        if (requireNamespace("akima", quietly = TRUE)) {
            tmp = akima::bicubic(U,V,img,x2,y2)$z
        } else {
            grid.list<- cbind(x2, y2)
            obj<-list(x=U, y=V, z=img)
            tmp = interp.surface(obj, grid.list)
        }
        img_new = matrix(tmp,nrow=m,byrow=F)
    } else {
        for (i in 1:d) {
            if (requireNamespace("akima", quietly = TRUE)) {
                tmp = akima::bicubic(U,V,img[,,i],x2,y2)$z
            } else {
                grid.list<- cbind(x2, y2)
                obj<-list(x=U, y=V, z=img[,,i])
                tmp = interp.surface(obj, grid.list)
            }
            img_new[,,i] = matrix(tmp,nrow=m,byrow=F)
        }
    }

    return (img_new)
}


apply_gam_gamid <- function(gamid, gaminc){
    m = dim(gamid)[1]; n = dim(gamid)[2]; d = dim(gamid)[3]
    U = seq(0,1,length.out=m)
    V = seq(0,1,length.out=n)
    x2 = c(gaminc[,,2])
    y2 = sort(c(t(gaminc[,,1])))
    gam_cum = array(0,dim=c(m,n,d))
    for (i in 1:d) {
        if (requireNamespace("akima", quietly = TRUE)) {
            tmp = akima::bicubic(U,V,gamid[,,i],x2,y2)$z
        } else {
            grid.list<- cbind(x2, y2)
            obj<-list(x=U, y=V, z=gamid[,,i])
            tmp = interp.surface(obj, grid.list)
        }
        gam_cum[,,i] = matrix(tmp,nrow=m,byrow=F)
    }

    return (gam_cum)
}


compgrad2D <- function(f){
    dims = ndims(f)
    if (dims < 3){
        n = dim(f)[1]; t = dim(f)[2]
        d = 1
        dfdu = matrix(0,n,t)
        dfdv = matrix(0,n,t)
    } else {
        n = dim(f)[1]; t = dim(f)[2]; d = dim(f)[3]
        dfdu = array(0,dim=dim(f))
        dfdv = array(0,dim=dim(f))
    }

    out = .Call('find_grad_2D', PACKAGE = 'fdasrvf', dfdu, dfdv, f, n, t, d)

    return(list(dfdu=out$dfdu,dfdv=out$dfdv))
}


gram_schmidt_c <- function(b){
    m = dim(b)[1]; n = dim(b)[2]; tmp = dim(b)[3]; N = dim(b)[4]
    ds = 1./((m-1)*(n-1))

    cnt = 1
    G = array(0,dim=dim(b))
    G[,,,cnt] = b[,,,cnt]

    out = compgrad2D(G[,,,cnt])
    dvx = out$dfdu; dvy = out$dfdv
    l = (t(c(dvx))%*%c(dvx)*ds + t(c(dvy))%*%c(dvy)*ds)[1]

    for (i in 2:N) {
        G[,,,i] = b[,,,i]
        out = compgrad2D(G[,,,i])
        dv1x = out$dfdu; dv1y = out$dfdv
        for (j in 1:(i-1)) {
            out = compgrad2D(G[,,,j])
            dv2x = out$dfdu; dv2y = out$dfdv
            t1 = (t(c(dv1x))%*%c(dv2x)*ds+t(c(dv1y))%*%c(dv2y)*ds)[1]
            G[,,,i] = G[,,,i]-t1*G[,,,j]
        }

        v = G[,,,i]
        l = (t(c(v))%*%c(v)*ds)[1]

        if (l>0){
            cnt = cnt + 1
            G[,,,cnt] = G[,,,cnt]/sqrt(l)
        }
    }

    return (G)
}


jacob_imag <- function(F1){
    m = dim(F1)[1]; n = dim(F1)[2]; d = dim(F1)[3]

    out = compgrad2D(F1)
    mult_factor = matrix(0,m,n)
    if (d==2){
        mult_factor = out$dfdu[,,1]%*%out$dfdv[,,2] - out$dfdu[,,2]%*%out$dfdv[,,1]
        mult_factor = abs(mult_factor)
    } else if (d==3){
        mult_factor = (out$dfdu[,,2]%*%out$dfdv[,,3] - out$dfdu[,,3]%*%out$dfdv[,,2])^2
            + (out$dfdu[,,1]%*%out$dfdv[,,3] - out$dfdu[,,3]%*%out$dfdv[,,1])^2
            + (out$dfdu[,,1]%*%out$dfdv[,,2] - out$dfdu[,,2]%*%out$dfdv[,,1])^2
        mult_factor = sqrt(mult_factor)
    }

    return (mult_factor)
}


makediffeoid <- function(nrows,ncols){
    D = 2
    gamid = array(0,dim=c(nrows,ncols,D))

    out = meshgrid(seq(0,1,length.out=ncols),seq(0,1,length.out=nrows))

    gamid[,,1] = out$X
    gamid[,,2] = out$Y

    return (gamid)
}


imag_to_q <- function(F1){
    dims = ndims(F1)
    if (dims < 3){
        stop("Data dimension is wrong!")
    }
    d = dim(F1)[3]
    if (d<2){
        stop("Data dimension is wrong!")
    }
    q = F1;

    sqrtmultfact = sqrt(jacob_imag(F1))
    for (i in 1:d) {
        q[,,i] = sqrtmultfact*F1[,,i]
    }

    return (q)
}


comp_dist <- function(q1, q2){
    m = dim(q1)[1]
    n = dim(q1)[2]
    ds = 1./((m-1)*(n-1))

    tmp = q1-q2

    H = sum(sqrt(c(tmp)*c(tmp))*ds)

    return (H)
}


comp_energy <- function(q1, q2) {
    m = dim(q1)[1]; n = dim(q1)[2]; d = dim(q1)[3]
    ds = 1./(m-1)/(n-1)

    tmp = q1-q2
    H = (t(c(tmp))%*%c(tmp)*ds)[1]

    return (H)
}


update_gam <- function(qt, qm, b){
    v = qt-qm
    w = findphistar(qt, b)

    out = findupdategam(v, w, b)

    return (out$gamupdate)
}


findphistar <- function(q, b){
    n = dim(q)[1]; t = dim(q)[2]; d = dim(q)[3]
    K = dim(b)[4]

    w = array(0,dim=c(n,t,d,K))

    out = .Call('find_phistar', PACKAGE = 'fdasrvf', w, q, b, n, t, d, K)

    return (out)
}


findupdategam <- function(v, w, b) {
    m = dim(b)[1]; n = dim(b)[2]; D = dim(b)[3]; K = dim(b)[4]
    ds = 1./((m-1)*(n-1))

    innp = rep(0,K)

    gamupdate = array(0,dim=c(m,n,D))

    for (k in 1:K) {
        vt = w[,,,k]
        innp[k] = (t(c(v))%*%c(vt)*ds)[1]

        gamupdate = gamupdate + innp[k]*b[,,,k]
    }

    return(list(gamupdate=gamupdate,innp=innp))
}


check_crossing <- function(f) {
    n = dim(f)[1]; t = dim(f)[2]; D = dim(f)[3]

    if (D!=2)
        stop("Third dimension of first argument to be 2")

    diffeo = .Call('check_cross', PACKAGE = 'fdasrvf', f, n, t, D)

    if (diffeo == 0)
        is_diffeo = F
    else if (diffeo == 1)
        is_diffeo = T
    else
        stop("check_crossing error")

    return (is_diffeo)
}


formbasisTid <- function(M,m,n,basis_type="t"){
    out = meshgrid(seq(0,1,length.out=n),seq(0,1,length.out=m))

    idx = 1

    if (basis_type=="t"){
        b = array(0, dim=c(m,n,2,2*M))
        for (s in 1:M) {
            c1 = sqrt(2)*pi*s
            sPI2 = 2*pi*s

            b[,,1,idx] = matrix(0,m,n)
            b[,,2,idx] = (cos(sPI2*out$Y)-1)/c1


            b[,,1,idx+1] = (cos(sPI2*out$X)-1)/c1
            b[,,2,idx+1] = matrix(0,m,n)

            idx = idx + 2
        }
    }

    if (basis_type == "s"){
        b = array(0,dim=c(m,n,2,2*M))
        for (s in 1:M) {
            c1 = sqrt(2)*pi*s
            sPI2 = 2*pi*2

            b[,,1,idx] = matrix(0,m,n)
            b[,,2,idx] = sin(sPI2*out$X)/c1

            b[,,1,idx+1] = sin(sPI2*out$X)/c1
            b[,,2,idx+1] = matrix(0,m,n)

            idx = idx + 2
        }
    }

    if (basis_type == "i"){
        b = array(0, dim=c(m,n,2,M*M*8))
        for (s1 in 1:M) {
            s1PI2 = 2*pi*s1
            for (s2 in 1:M) {
                s2PI2 = 2*pi*s2
                c1 = pi*sqrt(s1^2+3*s2^2)
                b[,,1,idx] = (cos(s1PI2*out$X)-1)*(cos(s2PI2*out$Y))/c1
                b[,,2,idx] = matrix(0,m,n)
                b[,,1,idx+2] = ((cos(s1PI2*out$X)-1)*sin(s2PI2*out$Y))/c1
                b[,,2,idx+2] = matrix(0,m,n)
                c1 = pi*sqrt(s1^2+s2^2)
                b[,,1,idx+4] = sin(s1PI2*out$X)*(cos(s2PI2*out$Y))/c1
                b[,,2,idx+4] = matrix(0,m,n)
                b[,,1,idx+6] = (sin(s1PI2*out$X)*sin(s2PI2*out$Y))/c1
                b[,,2,idx+6] = matrix(0,m,n)

                c1 = pi*sqrt(s1^2+3*s2^2)
                b[,,1,idx+1] = matrix(0,m,n)
                b[,,2,idx+1] = (cos(s1PI2*out$Y)-1)*(cos(s2PI2*out$X))/c1
                b[,,1,idx+3] = matrix(0,m,n)
                b[,,2,idx+3] = ((cos(s1PI2*out$Y)-1)*sin(s2PI2*out$X))/c1
                c1 = pi*sqrt(s1^2+s2^2)
                b[,,1,idx+5] = matrix(0,m,n)
                b[,,2,idx+5] = sin(s1PI2*out$Y)*(cos(s2PI2*out$X))/c1
                b[,,1,idx+7] = matrix(0,m,n)
                b[,,2,idx+7] = (sin(s1PI2*out$Y)*sin(s2PI2*out$X))/c1

                idx = idx + 8
            }
        }
    }

    if (basis_type == "o"){
        b = array(0,dim=c(m,n,2,M*4))
        for (s in 1:M) {
            c1 = sqrt((4*pi^2*s^2+9)/6)
            sPI2 = 2*pi*s
            b[,,1,idx] = (cos(sPI2*out$X))*out$Y/c1
            b[,,2,idx] = matrix(0,m,n)
            b[,,2,idx+1] = (cos(sPI2*out$Y)-1)*out$X/c1
            b[,,1,idx+1] = matrix(0,m,n)
            b[,,1,idx+2] = sin(sPI2*out$X)*out$Y/c1
            b[,,2,idx+2] = matrix(0,m,n)
            b[,,2,idx+3] = sin(sPI2*out$Y)*out$X/c1
            b[,,1,idx+3] = matrix(0,m,n)
            idx = idx + 4
        }
    }

    return(list(b=b,U=out$X,V=out$Y))
}


run_basis <- function(Ft, M=10, basis_type="t", is_orthon=TRUE){
    m = dim(Ft)[1]
    n = dim(Ft)[2]

    gamid = makediffeoid(m,n)

    out = formbasisTid(M, m, n, basis_type)
    b = out$b
    if (is_orthon){
        b = gram_schmidt_c(b)
    }

    return(list(b=b,gamid=gamid))
}
