KenDist <-
function(Data, Penalty){
            N <- nrow(Data);
            Features <- ncol(Data);
            X <- numeric(Features);
            Y <- numeric(Features);
            Dist <- numeric(N*N);
            D = .C("KendallDistance", as.double(as.vector(t(Data))), as.integer(N), as.integer(Features), as.double(Penalty), as.double(X), as.double(Y), Dist = as.double(Dist));
            A <- D$Dist;
            ## A IS THE N^2 X 1 VECTOR FORM OF THE DISTANCE MATRIX. NEED TO CONVERT IT TO THE DISTANCE MATRIX.
            Distance <- matrix(0,N,N);
            for (i in 1:N){
                    Distance[i,] <- A[((i-1)*N+1):(i*N)];
             }
             return(Distance);
 }
