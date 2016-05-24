pij_markov <- function(y, current) {

    alpha <- current$alpha
    f <- current$f
    sigma2 <- current$sigma2
    
    phi_delta <- forward_backward.quantities(y, f, alpha, sigma2)

    phi <- phi_delta$phi ## backward quantities
    delta <- phi_delta$delta  ## forward quantities

    J <- ncol(f)
    pij <-  sapply(1:J, function(j) {
        (delta[,j]*phi[,j])/rowSums(delta*phi)
    })

    return(pij)
}


forward_backward.quantities <- function(y, f, alpha, sigma2){

    N <- nrow(f)
    J = ncol(f)
    
    ### calculating the forward quantities (delta_ij's)
    delta <- matrix(NA, nrow = N, ncol = J)
    delta[1,] <- sapply(1:J, function(j){
        dnorm(y[1], mean = f[1, j], sd = sqrt(sigma2[j]))*alpha$PI[j]
    })
    for(i in 2:N) {
        delta[i,] <- sapply(1:J, function(j){
            sum(delta[(i-1),]*alpha$A[,j])*dnorm(y[i], mean = f[i,j], sd = sqrt(sigma2[j]))
        })
    }

    ## calculating the backward quantities (phi_ij's)
    phi <- matrix(NA, nrow = N, ncol = J)
    phi[N,] <- rep(1, J)

    for(i in N:2){
        phi[i-1,] <- sapply(1:J, function(j){
            sum(alpha$A[j,] * dnorm(y[i], mean = f[i,],
                                    sd = sqrt(sigma2))*phi[i,])
        })
    }

    return(list(delta = delta, phi = phi))
}


pilj <- function(indice, i, y, pij, f, A, phi, sigma2){
    l <- indice[1]
    j <- indice[2]
    pilj <- (pij[(i-1), l]*A[l, j]*
             dnorm(y[i], mean = f[i, j],
                   sd = sqrt(sigma2[j]))* phi[i, j])/phi[(i-1), l]
    ## pilj <- (pij[i,l]*A[l,j]*
    ##          dnorm(y[(i+1)], mean = f[(i+1), j],
    ##                sd = sqrt(sigma2[j]))*phi[(i+1), j])/phi[i, l]
    return(pilj)
}


alpha_markov <- function(y, pij, current) {

    alpha <- current$alpha
    f <- current$f
    sigma2 <- current$sigma2
    
    phi_delta <- forward_backward.quantities(y, f, alpha, sigma2)
    phi <- phi_delta$phi

    J <- ncol(f)

    indices <- as.matrix(expand.grid(1:J, 1:J), colnames = NULL)
    colnames(indices) <- NULL

    N <- length(y)
    QsiMatrix <- matrix(NA, nrow = (N-1), ncol = ncol(f)^2)
    for (i in 2:N) {
        QsiMatrix[(i-1),] <- apply(indices, 1, pilj,
                                   i = i, y = y,
                                   pij = pij, f = f,
                                   A = alpha$A,
                                   phi = phi, sigma2 = sigma2)
    }

    #### QsiMatrix thas columns according to the lines of indices
    #### for J = 2 we have columns corresponding to pi11,pi21,pi12,pi22

    anew.vec <- rep(NA, nrow(indices))
    for(k in 1:nrow(indices)){
        anew.vec[k] <- sum(QsiMatrix[,k])/sum(pij[,indices[k, 1]][1:(N-1)])
    }


    Anew <- matrix(anew.vec, byrow = FALSE, nrow = J)

    PInew <- pij[1,]

    alphanew <- list(A = Anew, PI = PInew)

    return(alphanew)
}


criteria_alpha_markov <- function(current, update) {
    c.A <- as.vector(abs(current$A - update$A))
    ##c.PI <- abs(current$PI - update$PI)
    ## in this one realization case it is not
    ## worthy looking for PI as the estimation is very poor
}


function_cases1AND2 <- function(z_num, y, f, A, s2, Phi, Delta) {
    N <- length(y)
    ## z_num = (r, s, t)
    r_num <- z_num[1]
    s_num <- z_num[2]
    t_num <- z_num[3]

    case1_num <- ( Phi[,t_num][3:N]*dnorm(y[3:N], mean=f[,t_num][3:N], sd=sqrt(s2[t_num]))*A[s_num, t_num]
                  *dnorm(y[2:(N-1)], mean=f[,s_num][2:(N-1)], sd=sqrt(s2[s_num]))*A[r_num, s_num]*Delta[,r_num][1:(N-2)] )

    case1_den <- 0
    for(r in 1:2) {
        for(s in 1:2) {
            for(t in 1:2) {
                case1_den <- ( case1_den + ( Phi[,t][3:N]*dnorm(y[3:N], mean=f[,t][3:N], sd=sqrt(s2[t]))*A[s, t]
                                            *dnorm(y[2:(N-1)], mean=f[,s][2:(N-1)], sd=sqrt(s2[s]))*A[r, s]*Delta[,r][1:(N-2)] ) )
            }
        }
    }

    case1 <- sum(case1_num/case1_den)

    return(case1)

}


function_cases3AND4  <- function(z_num, y, f, A, s2, Phi, Delta) {
    N <- length(y)
    ## z_num = c(r, s, t, u)
    r_num <- z_num[1]
    s_num <- z_num[2]
    t_num <- z_num[3]
    u_num <- z_num[4]

    Case = 0
    for(D in 1:(N-3)) {

        case_num_D = ( Phi[,u_num][(2+D+1):N]*dnorm(y[(2+D+1):N], mean=f[,u_num][(2+D+1):N], sd=sqrt(s2[u_num]))*A[t_num, u_num]
                      *dnorm(y[(2+D):(N-1)], mean=f[,t_num][(2+D):(N-1)], sd=sqrt(s2[t_num]))*((A%^%D)[s_num, t_num])
                      *dnorm(y[2:(N-(D+1))], mean=f[,s_num][2:(N-(D+1))], sd=sqrt(s2[s_num]))*A[r_num, s_num]*Delta[,r_num][1:(N-D-2)] )

        case_den_D <- 0
        for(r in 1:2) {
            for(s in 1:2) {
                for(t in 1:2) {
                    for(u in 1:2) {

                        case_den_D <- ( case_den_D + ( Phi[,u][(2+D+1):N]*dnorm(y[(2+D+1):N], mean=f[,u][(2+D+1):N], sd=sqrt(s2[u]))*A[t, u]
                                                      *dnorm(y[(2+D):(N-1)], mean=f[,t][(2+D):(N-1)], sd=sqrt(s2[t]))*((A%^%D)[s, t])
                                                      *dnorm(y[2:(N-(D+1))], mean=f[,s][2:(N-(D+1))], sd=sqrt(s2[s]))*A[r, s]*Delta[,r][1:(N-D-2)] ) )

                    } } } }

        Case <- Case + sum(case_num_D/case_den_D)

    }
    return(Case)
}


stderr_markov <- function(y, current) {
    pij <- current$pij
    alpha <- current$alpha
    f <- current$f
    sigma2 <- current$sigma2
    
    N <- length(y)
    J <- ncol(f)
    
    if (J > 2) {
        warning('Standard error in Markov case is only available for J=2')
    }
    else {
        delta_phi <- forward_backward.quantities(y, f, alpha, sigma2)
        Delta <- delta_phi$delta
        Phi <- delta_phi$phi  
        
#######################################################################
        ## Calculating the observed information matrix ---> Ihat = EB - ES2  ##
#######################################################################

        ## EB is a diagonal matrix with entries:
                                        #EB <- matrix(0, ncol=J, nrow=J)
                                        #diag(EB) <- c((sum(pij[1:(N-1), 1])/(alpha$A[1, 2]*(1-alpha$A[1, 2]))), (sum(pij[1:(N-1), 2])/(alpha$A[2, 1]*(1-alpha$A[2, 1]))))

        ## ES2 has diagonal elements equal to diag(EB) + c(Ba, Bb) and, 
        ## therefore, diag(EB) will cancel out, no need to calculate it. 

        ## calculating the off diagonal elements of ES2, which I call B

        ## B is a sum of 4 sums ---> B = Si - Sii - Siii + Siv

        ## The first sum Si is the sum for all i diff j p(z_i-1=1, z_i=2, z_j-1=2, z_j=1|y)
        ## To calculate Si we need to calculate four other sums, cases


                                        # case 1 :
        ## sum i=3 till N p(z_i-2=1, z_i-1=2, z_i=1|y)
        case1 <- function_cases1AND2(c(1, 2, 1), y, f, A=alpha$A, sigma2, Phi, Delta)
                                        # case 2 :
        ## sum i=3 till N p(z_i-2=2, z_i-1=1, z_i=2|y)
        case2 <- function_cases1AND2(c(2, 1, 2), y, f, A=alpha$A, sigma2, Phi, Delta)

        ## case 3
        ## sum i=2 till N-2 p(z_i-1=1, z_i=2, z_i+D=2, z_i+D+1=1|y) for D > or equal to 1
        case3 <- function_cases3AND4(z_num=c(1, 2, 2, 1), y, f, A=alpha$A, sigma2, Phi, Delta)


        ## case 4
        ## sum i=2 till N-2 p(z_i-1=2, z_i=1, z_i+D=1, z_i+D+1=2|y)
        case4 <- function_cases3AND4(z_num=c(2, 1, 1, 2), y, f, A=alpha$A, sigma2, Phi, Delta)

                                        #Si is the sum for all i diff j p(z_i-1=1, z_i=2, z_j-1=2, z_j=1|y) divided by a12*a21
        Si = (case1 + case2 + case3 + case4)/(alpha$A[1, 2]*alpha$A[2, 1])

        ## Sii is the sum for all i diff j p(z_i-1=1, z_i=2 , z_j-1=2 , z_j=2 |y)
        ## Sii is the sum of  case1 + case3 + case4 
        Sii <- ( ( function_cases1AND2(c(1, 2, 2), y, f, A=alpha$A, sigma2, Phi, Delta)
                  + function_cases3AND4(z_num=c(1, 2, 2, 2), y, f, A=alpha$A, sigma2, Phi, Delta)
                  + function_cases3AND4(z_num=c(2, 2, 1, 2), y, f, A=alpha$A, sigma2, Phi, Delta) )/(alpha$A[1, 2]*(1-alpha$A[2, 1])) )

        ## Siii is the sum for all i diff j p(z_i-1=1, z_i=1 , z_j-1=2 , z_j=1 |y)
        ## Siii is the sum of  case2 + case3 + case4
        Siii <- ( ( function_cases1AND2(c(2, 1, 1), y, f, A=alpha$A, sigma2, Phi, Delta)
                   + function_cases3AND4(z_num=c(1, 1, 2, 1), y, f, A=alpha$A, sigma2, Phi, Delta)
                   +  function_cases3AND4(z_num=c(2, 1, 1, 1), y, f, A=alpha$A, sigma2, Phi, Delta))/((1-alpha$A[1, 2])*alpha$A[2, 1]) )

        ## Siv is the sum for all i diff j p(z_i-1=1, z_i=1 , z_j-1=2 , z_j=2 |y)
        ## Siv is the sum of  case3 + case4
        Siv <- ( ( function_cases3AND4(z_num=c(1, 1, 2, 2), y, f, A=alpha$A, sigma2, Phi, Delta)
                  + function_cases3AND4(z_num=c(2, 2, 1, 1), y, f, A=alpha$A, sigma2, Phi, Delta))/((1-alpha$A[1, 2])*(1-alpha$A[2, 1])) )

        B = Si - Sii - Siii + Siv

        ES2 <- matrix(B, nrow=J, ncol=J)

        ## now I will calculate the elements in the diagonal of ES2 ---> Ea2 and Eb2
        ## Ea2 = Aa + Ba, but Aa = diag(EB)[1], no need to calculate it

        ## Ba is the sum of Sai - 2*Saii + Saiv
        ## Sai is the sum for all i diff j p(z_i-1=1, z_i=2, z_j-1=1 , z_j=2 |y)
        ## case 3 + case 4 but they are equal
        Sai <- (2*function_cases3AND4(z_num=c(1, 2, 1, 2), y, f, A=alpha$A, sigma2, Phi, Delta))/(alpha$A[1, 2]*alpha$A[1, 2])

        ## Saii is the sum for all i diff j p(z_i-1=1, z_i=2 , z_j-1=1 , z_j=1 |y)
        ## case2 + case3 + case4
        Saii = (( function_cases1AND2(c(1, 1, 2), y, f, alpha$A, sigma2, Phi, Delta)
                 + function_cases3AND4(z_num=c(1, 2, 1, 1), y, f, alpha$A, sigma2, Phi, Delta)
                 + function_cases3AND4(z_num=c(1, 1, 1, 2), y, f, alpha$A, sigma2, Phi, Delta) )/(alpha$A[1, 2]*(1-alpha$A[1, 2])))


        ## Saiv is the sum for all i diff j p(z_i-1=1, z_i=1 , z_j-1=1 , z_j=1 |y)
        Saiv = ( ( 2*function_cases1AND2(c(1, 1, 1), y, f, alpha$A, sigma2, Phi, Delta) +
                  2*function_cases3AND4(z_num=c(1, 1, 1, 1), y, f, alpha$A, sigma2, Phi, Delta) ) / ((1-alpha$A[1, 2])*(1-alpha$A[1, 2])) )

        Ea2 = Sai -2*Saii + Saiv

        ## Eb2 = Ab + Bb, but again Ab = diag(EB)[2], no need to calculate it
        ## Bb is the sum of Sbi - 2*Sbii + Sbiv

        ## Sbi is the sum for all i diff j p(z_i-1=2, z_i=1 , z_j-1=2 , z_j=1 |y)
        Sbi <- (2*function_cases3AND4(z_num=c(2, 1, 2, 1), y, f, A=alpha$A, sigma2, Phi, Delta))/(alpha$A[2, 1]*alpha$A[2, 1])

        ## Sbii is the sum for all i diff j p(z_i-1=2, z_i=1 , z_j-1=2 , z_j=2 |y)
        Sbii = ( ( function_cases1AND2(c(2, 2, 1), y, f, A=alpha$A, sigma2, Phi, Delta)
                  + function_cases3AND4(z_num=c(2, 1, 2, 2), y, f, A=alpha$A, sigma2, Phi, Delta)
                  + function_cases3AND4(z_num=c(2, 2, 2, 1), y, f, A=alpha$A, sigma2, Phi, Delta) )/(alpha$A[2, 1]*(1-alpha$A[2, 1])) )

        ## Sbiv is the sum for all i diff j p(z_i-1=2, z_i=2 , z_j-1=2 , z_j=2 |y)
        Sbiv = (( 2*function_cases1AND2(c(2, 2, 2), y, f, A=alpha$A, sigma2, Phi, Delta) +
                 2*function_cases3AND4(z_num=c(2, 2, 2, 2), y, f, A=alpha$A, sigma2, Phi, Delta) ) / ((1-alpha$A[2, 1])*(1-alpha$A[2, 1])) )

        Eb2 <- Sbi - 2*Sbii + Sbiv

        diag(ES2) <- c(Ea2, Eb2)

        Ihat = - ES2  ## EB - ES2, but because we are considering the terms that cancel out, it is just our ES2
        Ihat_inv = solve(Ihat)
        se <- sqrt(diag(Ihat_inv))

        return(se)
    }

}
