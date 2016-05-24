#' Computes abod measure for outlier detection
#' @param data - Dataframe in which to compute angle-based outlier factor.
#' @param method - Method to perform. 'complete' will use the entire dataset (O(n³) complexity) to
#' compute abof. 'randomized' will use a random sample of the data of size 'n_sample_size'. 'knn' will
#' compute abof among 'k' nearest neighbours. Please note that 'knn' has to compute an euclidean distance
#' matrix before computing abof.
#' @param n_sample_size - Number of random observations to choose in randomized method.
#' @param k - Number of nearest neighbours to choose in knn method.
#' @return Returns angle-based outlier factor for each observation. A small abof respect the others would indicate
#' presence of an outlier.
#' @examples
#' abod(mtcars[,1:4], method = "complete")
#' abod(faithful, method = "randomized", n_sample_size = 50)
#' abod(faithful, method = "knn", k = 20)
#' @references
#' [1] Angle-Based Outlier Detection in High-dimensional Data. KDD 2008.  Hans-Peter Kriegel, Matthias Schubert, Arthur Zimek. (http://www.dbs.ifi.lmu.de/Publikationen/Papers/KDD2008.pdf)
#' @export


abod <- function(data, method = "complete", n_sample_size = trunc(nrow(data)/10), k = 15){

  n <- nrow(data)
  
  if(method == "complete"){
    # Preallocate space for abof for each point
    abof = double(n)
    # Iterate through all (i,j,k)
    for(i in 1:n){
      # Preallocate space for angles and norms
      index = 1
      angles = c()
      # Get vertice
      a = as.numeric(data[i, ])
      for(j in 1:n){
        b = as.numeric(data[j, ])
        # Skip if identical
        if(identical(a, b)){
          next
        }
        for(k in 1:n){
          c = as.numeric(data[k, ])
          # Skip if identical (may cause trouble)
          if(identical(a, c) | identical(b, c)){
            next
          }
          # Some operations for clarity
          ab = b - a
          ac = c - a
          norm_ab = sqrt(sum(ab^2))
          norm_ac = sqrt(sum(ac^2)) # Using l2 norm
          dot_prod = sum(ab * ac)

          # Compute angle between vertice i, points j, k. Store it.
          angles = c(angles, dot_prod/(norm_ab * norm_ac))
          index = index + 1


        }
      }
      abof[i] = var(angles)
      if(i %% 10 == 0){
        cat("Done: ", i, "/", n, "\n")
      }
    }
    return(abof)
  }

  else if(method == "randomized"){
    aprox_abof = double(n)
    
    for(i in 1:n){
        a = as.numeric(data[i, ])
        selected_indexes = sample(1:n, n_sample_size)
        selected_data = data[selected_indexes, ]
        index = 1
        angles = c()
        for(j in 1:n_sample_size){
          b = as.numeric(selected_data[j, ])
          if(identical(a, b)){
            next
          }
          for(k in 1:n_sample_size){
            c = as.numeric(selected_data[k, ])
            if(identical(a, c) | identical(b, c)){
              next
            }
            # Perform operations
            ab = b - a
            ac = c - a
            norm_ab = sqrt(sum(ab^2))
            norm_ac = sqrt(sum(ac^2)) # Using l2 norm
            dot_prod = sum(ab * ac)
            
            # Compute angle between vertice i, points j, k. Store it.
            angles = c(angles, dot_prod/(norm_ab * norm_ac))
            index = index + 1            
          }
        }
        aprox_abof[i] = var(angles)
        if(i %% 10 == 0){
          cat("Done: ", i, "/", n, "\n")
        }
    }
    return(aprox_abof)
  }
  
  else if(method == "knn"){
    aprox_abof = double(n)
    distances = as.matrix(daisy(data))
    
    for(i in 1:n){
      # Find k nearest neighbours and subset data
      a = as.numeric(data[i, ])
      selected_indexes = sort(as.numeric(distances[i, ]), decreasing = T, index.return = TRUE)$ix[2:(k+1)]
      selected_data = data[selected_indexes, ]
      index = 1
      angles = c()
      for(j in 1:k){
        b = as.numeric(selected_data[j, ])
        if(identical(a, b)){
          next
        }
        for(k in 1:k){
          c = as.numeric(selected_data[k, ])
          if(identical(a, c) | identical(b, c)){
            next
          }
          # Perform operations
          ab = b - a
          ac = c - a
          norm_ab = sqrt(sum(ab^2))
          norm_ac = sqrt(sum(ac^2)) # Using l2 norm
          dot_prod = sum(ab * ac)
          
          # Compute angle between vertice i, points j, k. Store it.
          angles = c(angles, dot_prod/(norm_ab * norm_ac))
          index = index + 1            
        }
      }
      aprox_abof[i] = var(angles)
      if(i %% 10 == 0){
        cat("Done: ", i, "/", n, "\n")
      }
    }
    return(aprox_abof)
  }
    
}
