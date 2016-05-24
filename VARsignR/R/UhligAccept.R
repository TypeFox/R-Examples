UhligAccept <-
function(Q, first, last, constrained, impulses){#ok
  for(k in first:last){#ok
    ik <- impulses[k, , ]%*%Q#ok
    for(i in 1:length(constrained)){#ok
      if(constrained[i]<0){#ok
        value <- ik[-1.0 * constrained[i]]#ok
      }else{#ok
        value <- -1.0 * ik[constrained[i]]#ok
      }#ok
      if(value>0.0){#ok
        if(k==first & i==1){#ok
          Q <- -1.0 * Q#ok
          ik <- -1.0 * ik#ok
        }else{#ok
          acc <- 0
          uar <- list(Q=Q, acc=acc, ika=ik)
          return(uar)
        }#ok
      }#ok
    }#end i #ok
  }#end k #ok
  acc <- 1
  uar <- list(Q=Q, acc=acc, ika=ik)
  return(uar)
}
