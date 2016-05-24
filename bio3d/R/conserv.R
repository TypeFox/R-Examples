`conserv` <-
function(x,
         method=c("similarity","identity","entropy22","entropy10"),
         sub.matrix=c("bio3d","blosum62","pam30","other"),
         matrix.file=NULL,
         normalize.matrix = TRUE) {

  method <- match.arg(method)
  sub.matrix <- match.arg(sub.matrix)


  ##  cat(paste("Options are, method=",method,", matrix=",sub.matrix,
  ##            ", norm=",normalize.matrix),"\n" )

  if(is.list(x)) x=x$ali

  aa <- c("V","I","L","M",  "F","W","Y",  "S","T",
          "N","Q",  "H","K","R",  "D","E",
          "A","G",  "P",  "C",  "-","X")
  composition <- table(x)
  unk <- composition[!names( composition  ) %in% aa]
  if(length(unk) > 0) {
    warning(paste("non standard residue code:",names(unk),"mapped to X\n  "))
    for(i in 1:length(unk))
      x[x==names(unk[i])]="X"
  }
    
  if(method == "entropy10") { return(entropy(x)$H.10.norm) }
  if(method == "entropy22") { return(entropy(x)$H.norm) }
  if(method == "identity") {
    ## Identity (exact matchs score 1)
    freq.aa <- apply(x,2,
                     function(i){
                       i.freq <- table(i[i!="-"])
                       if(length(i.freq)==0) { return(0) } else {
                         return( max(table(i[i!="-"])) ) }
                     } )
    return( freq.aa/nrow(x) )
  }


  if(method == "similarity") {
    #####cat(sub.matrix)
    ## Pairwise matches are assigned score from a 'similarity matrix'
    if(sub.matrix=="other") {
      if(is.null(matrix.file)) stop("Missing argument: similarity requires a 'matrix.file'")
      mat.file <- matrix.file
    } else {
      ##mat.file <- system.file("matrices/similarity.mat", package="bio3d")
      mat.file <- system.file(paste("matrices/",sub.matrix,".mat",sep=""), package="bio3d")
      ##mat.file <- paste("matrices/",sub.matrix,".mat",sep="")
    }
    mat <- read.table(mat.file)
    colnames(mat)[24]="-"

    if(normalize.matrix) {
      ## Karlin Normalize
      o.mat <- mat
      n <- nrow(o.mat)
      for(a.ind in 1:n) {
        for(b.ind in 1:n) {
          ab <- o.mat[a.ind,b.ind]
          aa <- o.mat[a.ind,a.ind]
          bb <- o.mat[b.ind,b.ind]
          aabb <- aa*bb
          if(ab==0 && aabb==0) {
            mat[a.ind,b.ind] <- 0
          } else {
            if(aabb<0) {
              mat[a.ind,b.ind] <- ab / -sqrt(abs(aabb))
            } else {
              mat[a.ind,b.ind] <- ab / sqrt(aabb)
            }
          }
        }
      }
    }

    scorecol <- function(col, mat) {
      freq.aa   <- table(col)
      unique.aa <- names(freq.aa)
      missing.aa <- unique.aa[!unique.aa %in% colnames(mat)]
      
      count <- 0; score <- 0

      for(i in 1:length(unique.aa)) {
        aa.i <- unique.aa[i]; freq.i <- freq.aa[i]
        for(j in i:length(unique.aa)) {
          aa.j <- unique.aa[j]; freq.j <- freq.aa[j]
          ##sim <- mat[aa.i,aa.j]
          if(length(missing.aa)>0) {
            if(i==missing.aa || j==missing.aa) {
              sim <- 0
            } else {
              sim <- mat[aa.i,aa.j]
            }
          } else { sim <- mat[aa.i,aa.j] }

          ## number of comparisons
          if(aa.i == aa.j) {
            ncmp <- freq.i * (freq.i - 1)/2
          } else {
            ncmp <- freq.i * freq.j
          }
          count <- count + ncmp
          score <- score + (ncmp * sim)
        }
      }
      return(score/count)
    }


    return( apply(x, 2, scorecol, mat=mat) )
  }

}

