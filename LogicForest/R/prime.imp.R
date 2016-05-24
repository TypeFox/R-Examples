prime.imp <-
function (tree, data, Xs) 
{
  mat<-TTab(data, tree, Xs=Xs)
  mat[mat == 0] <- -1
  n.var <- ncol(mat)
  if (is.null(colnames(mat))) 
       colnames(mat) <- paste("X", 1:n.var, sep = "")
  ia <- p.combos(n.var)
  ia.rowS <- rowSums(ia)
  vec.primes <- character(0)
  list.pimps<-list()
  list.cover <- list()
  mat.in <- NULL
  name.paste <- function(x) 
     {
     x <- x[x != ""]
     paste(x, collapse = " & ")
     }
  tmp.list<-vector("list",n.var)
  for (i in 1:n.var) 
     {
     pairt <- p.combos(i, conj = -1)
     n.p <- nrow(pairt)
     ia2 <- matrix(ia[ia.rowS == i, ], ncol = n.var)
     tmp <- matrix(0, nrow(ia2) * n.p, n.var)
     for (j in 1:nrow(ia2)) tmp[((j - 1) * n.p + 1):(j * n.p), 
            ia2[j, ] == 1] <- pairt
        if (length(vec.primes) > 0) {
            tmp9 <- tmp %*% t(mat.in)
            tmp10 <- diag(mat.in %*% t(mat.in))
            tmp11 <- t(tmp10) %x% rep(1, nrow(tmp)) == tmp9
            tmp.in <- which(rowSums(tmp11) == 0)
            tmp <- tmp[tmp.in, ]
        }
     tmp2 <- tmp %*% t(mat) == i
     ids <- which(rowSums(tmp2) == 2^(n.var - i))
     tmp <- tmp[ids, ]
     tmp.list[[i]]<-tmp
     if (length(ids) > 0) {
            mat.in <- rbind(mat.in, tmp)
            for (k in ids) list.cover[[length(list.cover) + 1]] <- which(tmp2[k, ])
            mat.names <- matrix(rep(colnames(mat), e = length(ids)), 
                ncol = n.var)
            mat.names[tmp == 0] <- ""
    mat.pimps<-vector("list", nrow(mat.names))
    for (m in 1:nrow(mat.names)){
    mat.pimps[[m]]<-which(names(data)%in%mat.names[m,])}
    list.pimps<-append(list.pimps, mat.pimps)
       mat.names[tmp == -1] <- paste("!", mat.names[tmp == -1], sep = "")
            tmp.prime <- apply(mat.names, 1, name.paste)
            vec.primes <- c(vec.primes, tmp.prime)
        }
      cover <- unique(unlist(list.cover))
      if (length(cover) == nrow(mat)) 
            break
      }
    tmp.mat<-NULL
    for (i in 1:length(tmp.list)){
    tmp.mat<-rbind(tmp.mat, tmp.list[[i]])
    }
  vec.pimpvars<-sort(unique(unlist(list.pimps)))
  n.prime <- length(vec.primes)
  listPI <- list(vec.primes=vec.primes, tmp.mat=tmp.mat, vec.pimpvars=vec.pimpvars, list.pimps=list.pimps)
  class(listPI) <- "primeImp"
  listPI
}
