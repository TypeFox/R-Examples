##
## Function to group means
##

make.TukeyC.groups <- function(x)
{
  j <- 0
  i <- 1

  if(sum(apply(x,
               1,
               function(x,y) sum(x|x))) == 0) {
    res <- rep(letters[1],
               length(x[i, ]))

    names(res) <- rownames(x)

    return(res)
  }  

  res <- matrix('',
                nrow=nrow(x))

  while(i <= length(x[1, ])) {
    j <- j + 1

    if(i == 1) {
      ga <- x[i, ]           # Primeira linha de x

      gr <- rep(letters[1],
                length(ga))  # Vetor com a letra "a"

      res[1:length(gr[!ga]), 1] <- gr[!ga]  # Marca o grupo "a"

      b <- apply(x,
                 1,
                 function(x) sum(ga*x) == sum(ga*ga))

      i <- length(gr[!ga]) + 1                    # O próximo tratamento usado para marcar o grupo "b"

      j <- j + 1
    }

    g  <- x[i,]                                   # Linha do primeiro tratamento significativo com o grupo anterior

    gr <- rep(letters[j],
              length(x[i, ]))                     # Vetor com a letra do próximo tratamento

    g[(sum(b*b) + 1) : length(g)] <- FALSE

    g[length(g[g]) + 1] <- TRUE                   # Marcador para encontrar o final do grupo

    b <- apply(x,
               1,
               function(x) sum(g*x) == sum(g*g))  # Marcador do tratamento onde

    # Termina o grupo a ser marcado seguindo a lógica citada
    b <- !b

    g[length(g[g])] <- FALSE

    res <- cbind(res,
                 c(rep("",
                       length(gr[g])),
                   gr[(length(gr[g]) + 1) : sum(b*b)],
                   rep("",
                       (length(g) - sum(b*b)))))  # Coluna de res com o grupo marcado

    i <- sum(b*b) + 1
  }

  rownames(res) <- rownames(x)

  return(res)
}
