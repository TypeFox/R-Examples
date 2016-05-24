RPGenerate <-
function(p=100 #higher dimension
                    ,d=10   #lower dimension
                    ,method = "Haar" # Random projection method
                     )
    {
     if (p < d) stop("d must be less than p") 
        if (method == "Haar")
            {
                Z <- matrix(rnorm(p*d,0,1),p,d)
                Q <- qr.Q(qr(Z))
                R <- qr.R(qr(Z))
                D <- diag(sign(diag(R)))
                Q <- Q%*%D
            }           
        if (method == "axis")
            {
                S <- sample(1:p,d)
                Q <- matrix(0,p,d)
                for (D in 1:d)
                    {
                        Q[S[D],D] <- 1
                    }
            }
       if (method == "other")
            {
                Q <- matrix(0,p,d)
            }
        return(Q)
    }
