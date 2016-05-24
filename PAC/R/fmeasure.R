#' Compute the F measure between the ground truth and the estimated label
#' 
#' @param  g  the ground truth
#' @param  t  estimated labels
#' @return f  the F measure 
#' @export
  fmeasure <- function(g, t) {
  level = unique(g);
  C = length(level);
  if (length(unique(t)) != C) {
    print("Warning: number of cluster does not match with number of class")
  }
  S = matrix(0, C, 1)
  for (i in 1:C) {
    S[i] = sum(t == i)
  }
  
  P = matrix(0, C, 1)
  N = 0
  
  for (i in 1:C) {
    if (S[i] > 0) {
      P[i] = 0.5*S[i]*(S[i]-1)
    } else {
      P[i] = 0;
    }
  }
  
  for (i in 1:C) {
    for (j in (i+1):C) {
      if (j > C) break;
      N = N + S[i]*S[j]
    }
  }
  
  
  TP = matrix(0, C, 1)
  FP = matrix(0, C, 1)
  FN = matrix(0, C, 1)
  TN = 0
  
  
  for (i in 1:C) {
    idx = (t == i)
    subcluster = g[idx] 
    sublevel = unique(subcluster)
    numlabel = length(sublevel);
    subnum = matrix(0, numlabel,1);
    if (numlabel == 0) {
      next
    }
    for (j in 1:numlabel) {
      subnum[j] = sum(subcluster == sublevel[j]);
      if (subnum[j] > 0) {
        TP[i,1] = TP[i,1] + subnum[j]*(subnum[j] - 1)/2;
      }
      
    }
    
    for (j in 1:numlabel) {
      for (k in (j+1):numlabel) {
        if (k > numlabel) break;
        FP[i,1] = FP[i,1]+ subnum[j]*subnum[k];
      }
    }
  }
  
  
  NG = matrix(0, C, C)
  for(k in 1:C) {
    idx = t == k
    subcluster = g[idx];
    for (i in 1:C) {
      NG[i,k] = sum(subcluster == i);
      
    }
  }
  
  
  for (i in 1:C) {
    for (j in 1:C) {
      for (k in (j+1):C) {
        if (k > C) break;
        FN[i,1] = FN[i,1]+ NG[i,j] * NG[i,k];
      }
    }
    
    for (j in (i+1):C) {
      if (j > C) break
      for (k in 1:C) {
        for (l in 1:C) {
          if (k == l) next
          TN = TN + NG[k,i]*NG[l,j];
        }
      }
    }
  }
  
  Precision = sum(TP)/ sum(P) ;
  Recall = sum(TP) / (sum(TP) + sum(FN));
  
  f = 2*Precision*Recall / (Precision + Recall);
  return(f)
}
