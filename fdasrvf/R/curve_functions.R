calculatecentroid <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)

    betadot = matrix(0, n, T1)
    for (i in 1:n){
        betadot[i,] = gradient(beta[i,], 1.0/(T1-1))
    }

    normbetadot = rep(0, T1)
    integrand = matrix(0, n, T1)
    for (i in 1:T1){
        normbetadot[i] = pvecnorm(betadot[,i], 2)
        integrand[,i] = beta[,i] * normbetadot[i]
    }
    scale = trapz(seq(0,1,length.out=T1), normbetadot)
    centroid = trapz(seq(0,1,length.out=T1), integrand, 2)/scale

    return(centroid)
}


innerprod_q2 <- function(q1, q2){
    T1 = ncol(q1)
    val = sum(sum(q1*q2))/T1

    return(val)
}


find_best_rotation <- function(q1, q2){
    eps = .Machine$double.eps
    n = nrow(q1)
    T1 = ncol(q1)
    A = q1%*%t(q2)
    out = svd(A)
    s = out$d
    U = out$u
    V = out$v
    if (abs(det(U)*det(V)-1) < (10*eps)){
        S = diag(1,n)
    } else {
        S = diag(1,n)
        S[,n] = -S[,n]
    }
    R = U*S*t(V)
    q2new = R%*%q2

    return(list(q2new=q2new, R=R))
}


calculate_variance <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)
    betadot = matrix(0, n, T1)
    for (i in 1:n){
        betadot[i,] = gradient(beta[i,], 1.0/(T1-1))
    }
    normbetadot = rep(0,T)
    centroid = calculatecentroid(beta)
    integrand = array(0, c(n,n,T1))
    time = seq(0,1,length.out=T1)
    for (i in 1:T1){
        normbetadot[i] = pvecnorm(betadot[,i],2)
        a1 = beta[,i] - centroid
        integrand[,,i] = a1 %*% t(a1) * normbetadot[i]
    }
    l = trapz(time, normbetadot)
    variance = trapz(time, integrand, 3)
    varaince = variance / l

    return(variance)
}


psi <- function(x, a, q){
    T1 = ncol(q)
    dim(a) = c(length(a),1)
    covmat = calculate_variance(x+repmat(a,1,T1))
    psi1 = covmat[1,1] - covmat[2,2]
    psi2 = covmat[1,2]
    psi3 = x[1,T1]
    psi4 = x[2,T1]

    return(list(psi1=psi1,psi2=psi2,psi3=psi3,psi4=psi4))
}


find_basis_normal <- function(q){
    n = nrow(q)
    T1 = ncol(q)

    f1 = matrix(0,n,T1)
    f2 = matrix(0,n,T1)
    for (i in 1:T1){
        f1[,i] = q[1,i]*q[,i]/pvecnorm(q[,i])+c(pvecnorm(q[,i]),0)
        f2[,i] = q[2,i]*q[,i]/pvecnorm(q[,i])+c(0,pvecnorm(q[,i]))
    }
    h3 = f1
    h4 = f2
    integrandb3 = rep(0,T1)
    integrandb4 = rep(0,T1)
    for (i in 1:T1){
        integrandb3[i] = t(q[,i])%*%h3[,i]
        integrandb4[i] = t(q[,i])%*%h4[,i]
    }
    b3 = h3 - q*trapz(seq(0,1,length.out=T1),integrandb3)
    b4 = h4 - q*trapz(seq(0,1,length.out=T1),integrandb4)

    basis = list(b3,b4)

    return(basis)
}


calc_j <- function(basis){
    b1 = basis[[1]]
    b2 = basis[[2]]
    T1 = ncol(b1)

    integrand11 = rep(0,T1)
    integrand12 = rep(0,T1)
    integrand22 = rep(0,T1)

    for (i in 1:T1){
        integrand11[i] = t(b1[,i])%*%b1[,i]
        integrand12[i] = t(b1[,i])%*%b2[,i]
        integrand12[i] = t(b2[,i])%*%b2[,i]
    }

    j = matrix(0,2,2)
    j[1,1] = trapz(seq(0,1,length.out=T1), integrand11)
    j[1,2] = trapz(seq(0,1,length.out=T1), integrand12)
    j[2,2] = trapz(seq(0,1,length.out=T1), integrand22)
    j[2,1] = j[1,2]

    return(j)
}


shift_f <- function(f, tau){
    n = nrow(f)
    T1 = ncol(f)
    fn = matrix(0, n, T1)
    fn[,1:(T1-1)] = circshift(f[,1:(T1-1)], c(0,tau))
    fn[,T1] = fn[,1]

    return(fn)
}


find_rotation_seed_coord <- function(beta1, beta2){
    n = nrow(beta1)
    T1 = ncol(beta1)
    q1 = curve_to_q(beta1)
    Ltwo = rep(0,T1)
    Rlist = array(0,c(n,n,T1))
    for (ctr in 1:T1){
        beta2n = shift_f(beta2, ctr)
        out = find_best_rotation(beta1, beta2n)
        q2new = curve_to_q(out$q2new)
        Ltwo[ctr] = innerprod_q2(q1-q2new,q1-q2new)
        Rlist[,,ctr] = out$R
    }

    tau = which.min(Ltwo)
    O_hat = Rlist[,,tau]
    beta2new = shift_f(beta2,tau)
    beta2new = O_hat %*% beta2new

    return(list(beta2new=beta2new,O_hat=O_hat,tau=tau))
}


find_rotation_and_seed_q <- function(q1,q2){
    n = nrow(q1)
    T1 = ncol(q1)
    Ltwo = rep(0,T1)
    Rlist = array(0,c(n,n,T1))
    for (ctr in 1:T1){
        q2n = shift_f(q2,ctr)
        out = find_best_rotation(q1,q2n)
        Ltwo[ctr] = innerprod_q2(q1-out$q2new,q1-out$q2new)
        Rlist[,,ctr] = out$R
    }

    tau = which.min(Ltwo)
    O_hat = Rlist[,,tau]
    q2new = shift_f(q2,tau)
    q2new = O_hat %*% q2new

    return(list(q2new=q2new,O_hat=O_hat,tau=tau))
}


group_action_by_gamma <- function(q, gamma){
    n = nrow(q)
    T1 = ncol(q)
    gammadot = gradient(gamma, 1.0/T1)
    qn = matrix(0, n, T1)
    timet = seq(0, 1, length.out = T1)

    for (j in 1:n){
        qn[j,] = spline(timet, q[j,], xout=gamma)$y * sqrt(gammadot)
    }

    qn = qn/sqrt(innerprod_q2(qn,qn))

    return(qn)
}


group_action_by_gamma_coord <- function(f, gamma){
    n = nrow(f)
    T1 = ncol(f)
    fn = matrix(0, n, T1)
    timet = seq(0, 1, length.out = T1)

    for (j in 1:n){
        fn[j,] = spline(timet, f[j,], xout=gamma)$y
    }

    return(fn)
}


project_curve <- function(q){
    T1 = ncol(q)
    n = nrow(q)
    if (n==2){
      dt = 0.35
    }
    
    if (n==3) {
      dt = 0.2
    }
    
    epsilon =- 1e-6
    
    e = diag(1,n)
    iter = 1
    res = rep(1,n)
    J = matrix(0,n,n)
    s = seq(0,1,length.out=T1)
    qnorm = rep(0,T1)
    G = rep(0,n)
    C = rep(0,301)
    
    qnew = q
    qnew = qnew / sqrt(innerprod_q2(qnew,qnew))
    
    while (pvecnorm(res) > epsilon){
        if (iter > 300){
            break
        }
        
        # Compute Jacobian
        for (i in 1:n){
            for (j in 1:n){
              J[i,j]  = 3 * trapz(s, qnew[i,]*qnew[j,])
            }
        }
        J = J + diag(1,n)
        
        for (i in 1:T1){
            qnorm[i] = pvecnorm(qnew[,i])
        }
        
        # Compute the residue
        for (i in 1:n){
            G[i] = trapz(s,qnew[i,]*qnorm)
        }
        
        res = -1*G
        
        if (pvecnorm(res)<epsilon)
            break
        
        x = solve(J,res)
        C[iter] = pvecnorm(res)
        
        delG = find_basis_normal(qnew)
        tmp = 0
        for (i in 1:n){
            tmp = tmp + x[i]*delG[[i]]*dt
        }
        qnew = qnew + tmp
        
        iter = iter + 1
    }
    
    qnew = qnew / sqrt(innerprod_q2(qnew,qnew))

    return(qnew)
}


pre_proc_curve <- function(beta, T1=100){
    beta = resamplecurve(beta,T1)
    q = curve_to_q(beta)
    qnew = project_curve(q)
    x = q_to_curve(qnew)
    a = -1*calculatecentroid(x)
    dim(a) = c(length(a),1)
    betanew = x + repmat(a,1,T1)
    A = diag(1,2)

    return(list(betanew=betanew,qnew=qnew,A=A))
}


inverse_exp_coord <- function(beta1, beta2, mode="O"){
    T1 = ncol(beta1)
    centroid1 = calculatecentroid(beta1)
    dim(centroid1) = c(length(centroid1),1)
    beta1 = beta1 - repmat(centroid1, 1, T1)
    centroid2 = calculatecentroid(beta2)
    dim(centroid2) = c(length(centroid2),1)
    beta2 = beta2 - repmat(centroid2, 1, T1)

    q1 = curve_to_q(beta1)
    
    if (mode=="C"){
      isclosed = TRUE
    }

    # Iteratively optimize over SO(n) x Gamma using old DP
    out = reparam_curve(beta1, beta2, isclosed=isclosed)
    beta2 = out$R %*% shift_f(beta2,out$tau)
    gamI = invertGamma(out$gam)
    beta2 = group_action_by_gamma_coord(beta2, gamI)
    out = find_rotation_seed_coord(beta1, beta2)
    q2n = curve_to_q(out$beta2new)
    
    if (mode=="C"){
      q2n = project_curve(q2n)
    }

    # Compute geodesic distance
    q1dotq2 = innerprod_q2(q1, q2n)
    dist = acos(q1dotq2)

    if (q1dotq2>1){
        q1dotq2 = 1.
    }

    u = q2n - q1dotq2 * q1
    normu = sqrt(innerprod_q2(u,u))

    if (normu > 1e-4){
        v = u*acos(q1dotq2)/normu
    } else {
        v = matrix(0, 2, T1)
    }

    return(list(v=v,dist=dist))
}



inverse_exp <- function(q1, q2, beta2){
    T1 = ncol(q1)
    centroid1 = calculatecentroid(beta2)
    dim(centroid1) = c(length(centroid1),1)
    beta2 = beta2 - repmat(centroid1, 1, T1)

    # Optimize over SO(n) x Gamma
    beta1 = q_to_curve(q1)
    out = reparam_curve(beta1, beta2)
    gamI = invertGamma(out$gam)
    beta2 = out$R %*% shift_f(beta2, out$tau)

    # Applying optimal re-parameterization to the second curve
    beta2 = group_action_by_gamma_coord(beta2, gamI)
    q2 = curve_to_q(beta2)

    # Optimize over SO(n)
    out = find_rotation_and_seed_q(q1, q2)
    q2 = out$q2new

    # compute geodesic distance
    q1dotq2 = innerprod_q2(q1, q2)
    dist = acos(q1dotq2)

    if (q1dotq2>1){
        q1dotq2 = 1.
    }

    u = q2 - q1dotq2 * q1
    normu = sqrt(innerprod_q2(u,u))

    if (normu > 1e-4){
        v = u*acos(q1dotq2)/normu
    } else {
        v = matrix(0, 2, T1)
    }

    return(v)
}



gram_schmidt <- function(basis){
    b1 = basis[[1]]
    b2 = basis[[2]]

    basis1 = b1 / sqrt(innerprod_q2(b1,b1))
    b2 = b2 - innerprod_q2(basis1,b2)*basis1
    basis2 = b2 / sqrt(innerprod_q2(b2,b2))

    basis_o = list(basis1, basis2)

    return(basis_o)
}


project_tangent <- function(w, q, basis){
    w = w - innerprod_q2(w,q)*q
    bo = gram_schmidt(basis)

    wproj = w - innerprod_q2(w, bo[[1]])*bo[[1]] - innerprod_q2(w,bo[[2]])*bo[[2]]

    return(wproj)
}


scale_curve <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)
    normbetadot = rep(0,T1)
    betadot = matrix(0, n, T1)
    for (i in 1:n){
        betadot[i,] = gradient(beta[i,], 1.0/T1)
    }
    for (i in 1:T1){
        normbetadot[i] = pvecnorm(betadot[,i])
    }
    scale = trapz(seq(0,1,length.out=T1), normbetadot)
    beta_scaled = beta / scale

    return(list(beta_scaled=beta_scaled, scale=scale))
}


parallel_translate <- function(w, q1, q2, basis, mode='O'){
    wtilde = w - 2*innerprod_q2(w,q2) / innerprod_q2(q1+q2,q1+q2) * (q1+q2)
    l = sqrt(innerprod_q2(wtilde, wtilde))

    if (mode == 'C'){
        wbar = project_tangent(wtilde, q2, basis)
        normwbar = sqrt(innerprod_q2(wbar, wbar))
        if (normwbar>1e-4){
            wbar = wbar * l / normwbar
        }
    } else {
        wbar = wtilde
    }

    return(wbar)
}


rot_mat <- function(theta){
    O = matrix(0,2,2)
    O[1,1] = cos(theta)
    O[1,2] = -1*sin(theta)
    O[2,1] = sin(theta)
    O[2,2] = cos(theta)

    return(O)
}


karcher_calc <- function(beta, q, betamean, mu, mode="O"){
    if (mode=="C"){
        basis = find_basis_normal(mu)
    }
    # Compute shooting vector form mu to q_i
    out = inverse_exp_coord(betamean, beta)

    # Project to tangent space of manifold to obtain v_i
    if (mode=="O"){
        v = out$v
    } else {
        v = project_tangent(out$v, q, basis)
    }

    return(list(v=v,d=out$dist))
}
