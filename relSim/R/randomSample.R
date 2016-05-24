randomSample = function(nLoc, Freqs, rel = "UN", N = 10000){
    rel = toupper(rel)

    if(!grepl("(UN|FS|PC)",rel)){
        stop("rel must be one of 'UN', 'FS', or 'PC'")
    }

    f = Freqs$freqs[[nLoc]]

    if(rel == "UN"){
        U1 = matrix(sample(1:length(f), size = 2*N, replace = TRUE, prob = f),
                    ncol = 2, byrow = TRUE)
        U2 = matrix(sample(1:length(f), size = 2*N, replace = TRUE, prob = f),
                    ncol = 2, byrow = TRUE)

        if(any(U1[,1]>U1[,2])){
            i1 = which(U1[,1]>U1[,2])
            tmp = U1[i1,1]
            U1[i1,1] = S1[i1,2]
            U1[i1,2] = tmp
        }

        if(any(U2[,1]>U2[,2])){
            i1 = which(U2[,1]>U2[,2])
            tmp = U2[i1,1]
            U2[i1,1] = S1[i1,2]
            U2[i1,2] = tmp
        }

        return(cbind(U1,U2))
    }else if(rel == "FS"){
        S1 = matrix(sample(1:length(f), size = 2*N, replace = TRUE, prob = f),
                    ncol = 2, byrow = TRUE)

        if(any(S1[,1]>S1[,2])){
            i1 = which(S1[,1]>S1[,2])
            tmp = S1[i1,1]
            S1[i1,1] = S1[i1,2]
            S1[i1,2] = tmp
        }

        S2 = matrix(0, ncol = 2, nrow = N)

        A =  matrix(sample(1:length(f), size = 2*N, replace = TRUE, prob = f),
                    ncol = 2, byrow = TRUE)
        Idx = sample(1:4, size = N, replace = TRUE)

        for(j in 1:N){
            i = Idx[j]

            switch(Idx[j],
               {S2[j,] = S1[j,]},
               {S2[j,] = c(S1[j,1], A[j,1])},
               {S2[j,] = c(A[j,2], S1[j,2])},
               {S2[j,] = A[j,]}
                   )
        }

        if(any(S2[,1]>S2[,2])){
            i1 = which(S2[,1]>S2[,2])
            tmp = S2[i1,1]
            S2[i1,1] = S2[i1,2]
            S2[i1,2] = tmp
        }

        return(cbind(S1,S2))
    }else if(rel == "PC"){
        P1 = matrix(sample(1:length(f), size = 2*N, replace = TRUE, prob = f),
                    ncol = 2, byrow = TRUE)

        if(any(P1[,1]>P1[,2])){
            i1 = which(P1[,1]>P1[,2])
            tmp = P1[i1,1]
            P1[i1,1] = P1[i1,2]
            P1[i1,2] = tmp
        }

        C1 = matrix(0, ncol = 2, nrow = N)

        A = sample(1:length(f), size = N, replace = TRUE, prob = f)
        U = runif(N)

        for(j in 1:N){
            if(U[i] < 0.5){
                C1[j,] = c(P1[j,1],A[j])
            }else{
                C1[j,] = c(A[j], P1[j,2])
            }
        }

        if(any(C1[,1]>C1[,2])){
            i1 = which(C1[,1]>C1[,2])
            tmp = C1[i1,1]
            C1[i1,1] = C1[i1,2]
            C1[i1,2] = tmp
        }

        return(cbind(P1,C1))
    }
}
