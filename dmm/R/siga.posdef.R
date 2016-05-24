siga.posdef <-
function(siga, am, ctable,varopt="both",covopt="bound"){
#  siga.posdef() - make each matrix which is a row of siga positive definite (if a var)
    smalleig <- 1e-09
    for (i in 1:am$v) {

      if(any(!is.na(match(rownames(siga)[i],ctable$allvar)))){
        # variance - check if pd
         if(varopt == "eigen" || varopt == "both") {
          eigsiga <- eigen(matrix(siga[i, ], am$l, am$l),symmetric=TRUE)
          neg <- FALSE
          for (j in 1:am$l) {
              if (eigsiga$values[j] <= 0) {
                  eigsiga$values[j] <- smalleig
                  neg <- TRUE
              }
          }
          if (neg) {
              siga[i, ] <- matrix(eigsiga$vectors %*% diag(eigsiga$values, 
                  nrow = length(eigsiga$values)) %*% solve(eigsiga$vectors), 
                  1, am$l * am$l)
          }
         }
         if(varopt == "nearPD" || varopt == "both"){
           pdmat <- nearPD(matrix(siga[i, ], am$l, am$l),ensureSymmetry=T)
          #overwrite siga
           siga[i, ] <- matrix(pdmat$mat, 1, am$l * am$l)
         }
      }

      else  {
        # covariance - keep correlation in bounds
        if(rownames(siga)[i] == "CovG(Ia,Ma)") {
          # work out the c1x & c2y nos for Ia and Ma - ie for VarG(Ia) & VarG(Ma)
          c1 <- match("VarG(Ia)", rownames(siga))
          c2 <- match("VarG(Ma)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Ma,Ia)") {
          c1 <- match("VarG(Ma)", rownames(siga))
          c2 <- match("VarG(Ia)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Id,Md)") {
          c1 <- match("VarG(Id)", rownames(siga))
          c2 <- match("VarG(Md)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Md,Id)") {
          c1 <- match("VarG(Md)", rownames(siga))
          c2 <- match("VarG(Id)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Ia:a,Ma:a)") {
          c1 <- match("VarG(Ia:a)", rownames(siga))
          c2 <- match("VarG(Ma:a)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Ma:a,Ia:a)") {
          c1 <- match("VarG(Ma:a)", rownames(siga))
          c2 <- match("VarG(Ia:a)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Ia:d,Ma:d)") {
          c1 <- match("VarG(Ia:d)", rownames(siga))
          c2 <- match("VarG(Ma:d)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Ma:d,Ia:d)") {
          c1 <- match("VarG(Ma:d)", rownames(siga))
          c2 <- match("VarG(Ia:d)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Id:d,Md:d)") {
          c1 <- match("VarG(Id:d)", rownames(siga))
          c2 <- match("VarG(Md:d)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovG(Md:d,Id:d)") {
          c1 <- match("VarG(Md:d)", rownames(siga))
          c2 <- match("VarG(Id:d)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovGs(Ia,Ma)") {
          c1 <- match("VarGs(Ia)", rownames(siga))
          c2 <- match("VarGs(Ma)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovGs(Ma,Ia)") {
          c1 <- match("VarGs(Ma)", rownames(siga))
          c2 <- match("VarGs(Ia)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovE(I,M)") {
          c1 <- match("VarE(I)", rownames(siga))
          c2 <- match("VarE(M)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovE(M,I)") {
          c1 <- match("VarE(M)", rownames(siga))
          c2 <- match("VarE(I)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovE(I,M&!C)") {
          c1 <- match("VarE(I)", rownames(siga))
          c2 <- match("VarE(M&!C)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovE(M&!C,I)") {
          c1 <- match("VarE(M&!C)", rownames(siga))
          c2 <- match("VarE(I)", rownames(siga))
        } 
        if(rownames(siga)[i] == "CovE(I,M&C)") {
          c1 <- match("VarE(I)", rownames(siga))
          c2 <- match("VarE(M&C)", rownames(siga))
        }
        if(rownames(siga)[i] == "CovE(M&C,I)") {
          c1 <- match("VarE(M&C)", rownames(siga))
          c2 <- match("VarE(I)", rownames(siga))
        } 
#  cat("i = ",i," c1 = ",c1," c2 = ",c2,"\n")
         for( j in 1 : am$l) {
           jj <- (j-1)*am$l + j
           for(k in 1 : am$l) {
             kk <- (k-1)*am$l + k
             jk <- (j-1) * am$l + k
             if(covopt == "nearPD") {
               covmat <- matrix(c(siga[c1,jj],siga[i,jk],siga[i,jk],siga[c2,kk]),2,2)
               covmat <- nearPD(covmat,keepDiag=T)$mat
               # overwrite siga[i,jk],siga[c1,jj],siga[c2,kk]
               siga[i,jk] <- covmat[1,2]
#              siga[c1,jj] <- covmat[1,1]
#              siga[c2,kk] <- covmat[2,2]
             }
             else if(covopt == "bound"){
               clim <- sqrt(siga[c1,jj]*siga[c2,kk])
               if(siga[i,jk] < 0) {
                 siga[i,jk] <- max(siga[i,jk],-clim)
               }
               else if(siga[i,jk] > 0) {
                 siga[i,jk] <- min(siga[i,jk],clim)
               }
             }
           }
         }
      } # end else
    } # end for i
    return(siga)
}
