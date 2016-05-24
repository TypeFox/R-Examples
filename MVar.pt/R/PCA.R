PCA <- function(Data, Type = 1) {
  # Funcao Executa a Analise dos Componentes Principais - PCA 
  # Desenvolvida por Paulo Cesar Ossani em 07/2013
  
  # Entrada:
  # Data - Dados a serem a analizados
  # Type - 1 para analise utilizando a matriz de covariancia (default)
  #        2 para analise utilizando a matriz de correlacao
  
  # Retorna:
  # MatrixMC      - Matriz de Covariancia ou de  Correlacao conforme Type
  # MatrixAutoVlr - Matriz de Autovalores (Variancias) com as proporcoes e proporcoes acumuladas
  # MatrixAutoVec - Matriz de Autovetores - Componentes Principais
  # MatrixVCP     - Matriz da Covariancia dos Componentes Principais com as Variaveis Originais
  # MatrixCCP     - Matriz da Correlacao dos Componentes Principais com as Variaveis Originais
  # MatrixEsc     - Matriz com os escores dos Componentes Principais
  
  if (!is.data.frame(Data)) 
     stop("Entrada 'Data' esta incorreta, deve ser do tipo dataframe. Verifique!")
  
  if (Type!=1 && Type!=2) 
     stop("Entrada para 'Type' esta incorreta, deve ser numerica, sendo 1 ou 2. Verifique!")
  
  if (Type == 1)     # Considera a Matriz de Covariancia para a decomposicao
     MC <- cov(Data) # Matriz de Covariancia
  
  if (Type == 2)     # Considera a Matriz de Correlacao para a decomposicao
     MC <- cor(Data) # Matriz de Correlacao
  
  ## Encontrando a Matriz de Decomposicao Expectral
  MAV <- eigen(MC) # Encontra a matriz de autovalor e autovetor
  MAutoVlr <- MAV$values  # Matriz de Autovalores - Variancias
  MAutoVec <- MAV$vectors # Matriz de Autovetores - Componentes Principais
  
  ## Matriz das Variancias
  MEigen <- as.data.frame(matrix(NA, length(MAutoVlr), 3))
  rownames(MEigen) <- paste("Comp", 1:length(MAutoVlr))
  colnames(MEigen) <- c("Autovalor", "% da variancia","% acumulada da variancia")
  MEigen[, "Autovalor"] <- MAutoVlr
  MEigen[, "% da variancia"] <- (MAutoVlr/sum(MAutoVlr)) * 100
  MEigen[, "% acumulada da variancia"] <- cumsum(MEigen[,"% da variancia"])
  
  ## Matriz de Autovetores,ou seja, os Componentes Principais
  colnames(MAutoVec) <- paste("Comp.", 1:nrow(MC), sep = " ")
  rownames(MAutoVec) <- colnames(Data)  
  
  ## Covariancia dos Componentes Principais com as Variaveis Originais
  VCP <- diag(MAutoVlr,nrow(MC),ncol(MC))%*%t(MAutoVec)
  rownames(VCP) <- paste("Comp", 1:length(MAutoVlr))
  
  ## Correlacao dos Componentes Principais com as Variaveis Originais
  CCP <- diag(sqrt(MAutoVlr),nrow(MC),ncol(MC))%*%t(MAutoVec)%*%diag(1/sqrt(diag(MC)),nrow(MC),ncol(MC))
  colnames(CCP) <- colnames(Data) # Nomeia as linhas
  rownames(CCP) <- paste("Comp", 1:length(MAutoVlr))
  
  Esc = as.matrix(Data)%*%MAutoVec # Escores do componentes principais
  rownames(Esc) <- rownames(Data)
  
  Lista <- list(MatrixMC = MC, MatrixAutoVlr = MEigen,
                MatrixAutoVec = MAutoVec, MatrixVCP = VCP, 
                MatrixCCP = CCP, MatrixEsc = Esc)
  
  return(Lista)
}
