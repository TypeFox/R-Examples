`Simtocl` <-
function(Sim){
                n <- ncol(Sim)
                cl <- rep(0,n)
                groupn <- 1
                i <- 1
                validind <- cl == 0
                while(any(validind) & i <= n){
                    validind <- cl == 0
                    if(sum(Sim[validind, i] == 1) >= 1) {
                        cl[validind] <- (Sim[validind, i] == 1) * groupn
                        groupn <- groupn + 1
                    }
                    i <- i+1
                }
                cl
            }
