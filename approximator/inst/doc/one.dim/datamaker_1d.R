# This file creates a function that creates data suitable for analysis
# with the approximator package.  It is designed to be called by file
# apprex_1d.R

"generate.1d.observations" <-
function (D1, subsets, basis.fun, hpa, betas = NULL, export.truth=FALSE)
{

      if(is.null(betas)){
        betas <-
          rbind(c(1, 2),
                c(1, 1),
                c(1, 3))
        colnames(betas) <- c("const", "x")
        rownames(betas) <- paste("level", 1:3, sep = "")
      }
      
      if(export.truth){
        return(list(
                    hpa=hpa,
                    betas=betas
                    )
               )
      }


    sigma_squareds <- hpa$sigma_squareds
    B <- hpa$B
    rhos <- hpa$rhos
    delta <- function(i) {
      out <- rmvnorm(n = 1,
                     mean = basis.fun(D1[subsets[[i]], , drop =
                       FALSE]) %*% betas[i, ],
                     sigma = sigma_squareds[i] * corr.matrix(xold = D1[subsets[[i]], , drop = FALSE], pos.def.matrix = B[[i]])
                     )
      out <- drop(out)
      names(out) <- rownames(D1[subsets[[i]], , drop = FALSE])
      return(out)
    }
    
    use.clever.but.untested.method <- FALSE
    
    if(use.clever.but.untested.method){
      z1 <- delta(1)
      z2 <- delta(2) + rhos[1] * z1[match(subsets[[2]], subsets[[1]])]
      z3 <- delta(3) + rhos[2] * z2[match(subsets[[3]], subsets[[2]])]
      return(list(z1 = z1, z2 = z2, z3 = z3))
    } else {
      out <- NULL
      out[[1]] <- delta(1)
      for(i in 2:length(subsets)){
        out[[i]] <- delta(i) + rhos[i-1] *
    out[[i-1]][match(subsets[[i]], subsets[[i-1]])]
      }
      return(out)
    }
  }
