
# Normalization rules for A0.
# Method = "DistanceMLA"
#        = "DistanceMLAhat"
#        = "Euclidean"
#        = "PositiveDiagA"
#        = "PositiveDiagAinv"

"normalize.svar" <- function(A0unnormalized, A0mode,
                           method=c("DistanceMLA", "DistanceMLAhat",
                             "Euclidean", "PositiveDiagA",
                             "PositiveDiagAinv", "Unnormalized"),
                           switch.count=0)
  {
  A0normalized <- A0unnormalized   # Default starting case.

  if(method=="DistanceMLA")
    { indx <- which(diag(solve(A0unnormalized)%*%A0mode)<0)
      if(length(indx)!=0)
        { A0normalized[,indx] <- -A0unnormalized[,indx]
            switch.count <- switch.count+1
        }
    }
  else if(method=="DistanceMLAhat")
    { indx <- which(diag(solve(A0mode)%*%A0unnormalized)<0)
      if(length(indx)!=0)
         { A0normalized[,indx] <- -A0unnormalized[,indx]
           switch.count <- switch.count+1
         }
    }
  else if(method=="Euclidean")
    { Adiff <- (A0unnormalized - A0mode)^2
      Adiffn <- (-A0unnormalized - A0mode)^2
      cAdiff <- colSums(Adiff)
      cAdiffn <- colSums(Adiffn)

      # find the shorter distance
      indx <- which(cAdiffn < cAdiff)
      if(length(indx)!=0)
        { A0normalized[,indx] <- -A0unnormalized[,indx]
          switch.count <- switch.count+1
        }
    }
    else if(method=="PositiveDiagA")
      { indx <- which(diag(A0unnormalized)<0)
        if(length(indx)!=0)
          { A0normalized[,indx] <- -A0unnormalized[,indx]
            switch.count <- switch.count+1
          }
      }
    else if(method=="PositiveDiagAinv")
      { indx <- which(diag(solve(A0unnormalized))<0)
        if(length(indx)!=0)
          { A0normalized[,indx] <- -A0unnormalized[,indx]
            switch.count <- switch.count+1
          }
      }
    else if(method=="Unnormalized")
      {
        switch.count <- switch.count
      }
    else stop("No valid normalization rule selected.  Check 'method' argument.")

    return(list(A0normalized=A0normalized, switch.count=switch.count))
}
