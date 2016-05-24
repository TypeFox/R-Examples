MFA <- function(Data, Grupo, TipoGrupo = rep("n",length(Grupo)), NomeGrupos = NULL) {
  # Rotina para softwre R para uso do Metodo MFA para dados Quantitativos,
  # Categoricos e de Frequencia desenvolvida por Paulo Cesar Ossani em 
  # 03/2014

  # Entrada:
  # Data - Dados a serem analisados
  # Grupo - Numero de colunas para cada Grupo em 
  #         ordem seguindo a ordem dos Dados de 'Data'
  # TipoGrupo - "n" para dados Numericos (default)
  #             "c" para dados Categoricos
  #             "f" para dados de Frequencia
  # NomeGrupos - Nomes para cada Grupo
  
  # Retorna:
  # MatrixG  - Matriz com os tamanhos de cada grupo
  # MatrixNG - Matriz com os nomes de cada grupo
  # MatrixPLin - Matriz com os valores usados para balancear as linhas da Matriz Z
  # MatrixPCol - Matriz com os valores usados para balancear as colunas da Matriz Z
  # MatrixZ  - Matriz Concatenada e Balanceada
  # MatrixA  - Matriz com Autovalores (Variancias)
  # MatrixU  - Matriz U da decomposicao singular da Matriz Z
  # MatrixV  - Matriz V da decomposicao singular da Matriz Z
  # MatrixF  - Matriz Global dos Escores dos Fatores
  # MatrixEFG - Matriz dos Escores dos Fatores por Grupo
  # MatrixCCP - Matriz com a Correlacao dos Componentes Principais com os Grupos
  # MatrixEscVar - Matriz das Inercias Parciais/Escores das Variareis
 
  CA_MFA <- function(Data) {
    # Funcao que executa Analise de Correspondencia - CA 
    # nos dados e retorna o primeiro autovalor     
    # Esta funcao e usada na funcao que balanceia dados Categoricos
    
    # Entrada:
    # Data - Dados a serem analisados
    
    # Retorna:
    # Inercia - Primeiro auto valor
    
    SDados <- sum(Data) # Soma Total dos Dados
    
    MP <- as.matrix(Data/SDados) # Matriz da frequencia relativa
    
    r = apply(MP,1,sum) # Soma das Linhas
    
    c = apply(MP,2,sum) # Soma das Colunas
    
    Dr = diag(r) # Matriz diagonal de r
    
    Dc = diag(c) # Matriz diagonal de c
    
    MZ = diag(1/sqrt(diag(Dr)))%*%(MP - r%*%t(c))%*%diag(1/sqrt(diag(Dc))) # Matriz Z
    
    Mdvs <- svd(MZ) # Matriz de Decomposicao Valor Singular
    
    Md = diag(Mdvs$d) # Matriz diagonal Lambda
    
    Inercia = diag(Md%*%Md) # Calculo das inercias - Autovalores
    
    return(Inercia[1]) 
  }
  
  ICA <- function(Data,SomaLin) {
    # Funcao que retorna a Analise de Correspondencia Interna (ICA)
    # Esta funcao e usada na funcao que balanceia dados de Frequencia
    
    # Entrada:
    # Data    - Tabela de Frequencia dos Dados a serem analisados
    # SomaLin - Matriz com a soma geral das linhas dos dados de frequencia
    
    # Retorna:
    # MatriFR    - Matriz com as Frequencias Relativas
    # MatrixICA  - Matriz da Analise de Correspondencia Interna (ICA)
    # MatrixSLin - Matriz com as Somas das Linhas
    # MatrixSCol - Matriz com as Somas das Colunas
    
    SomaTot <- sum(SomaLin)
    
    MFR <- as.matrix(Data/SomaTot) # Matriz da Frequencia Relativa
    
    rTotal = SomaLin/SomaTot # Soma Geral das Linhas
    
    cTotal = apply(MFR,2,sum) # Soma Geral das Colunas
    
    rGrupo <- apply(MFR,1,sum) # Soma das Linhas do Grupo i
    
    cGrupo <- apply(MFR,2,sum) # Soma das Colunas do Grupo i
    
    ICA    <- MFR  # Matriz ICA do Grupo i
    
    P..t   <- sum(ICA) # Soma Total do Grupo i
    
    for (col in 1:ncol(ICA)) {    
      
      P.jt = cGrupo[col] # Soma Geral da coluna j da Tabela t
      
      for (lin in 1:nrow(ICA)) {
        
        Pijt = ICA[lin,col] # Elemento ij da Tabela t
        
        Pi.. = rTotal[lin]  # Soma Geral da linha i da Tabela Concatenada
        
        Pi.t = rGrupo[lin]  # Soma geral da linha i da Tabela t
        
        ICA[lin,col] =  1 / Pi.. * ( Pijt / P.jt - Pi.t / P..t) # Matriz ICA do Grupo i
      }      
    }
    
    Lista <- list(MatrixFR = MFR, MatrixICA = ICA, MatrixSLin = rTotal, MatrixSCol = cTotal)
    
    return(Lista)
  }
   
  MBQ <- function(DataQ,PondGeral) {  
    # Funcao que balanceia Dados quantitativos
    
    # Entrada:
    # DataQ - Dados a serem balanceados
    # PondGeral - usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
    
    # Retorna:
    # MZ   - Matriz Balanceada
    # PLin - Pesos das Linhas
    # PCol - Pesos das Colunas
    
    MZ <- NULL    # cria uma matriz Z nula
    
    PLin <- NULL   # Matriz com os pesos das linhas nula
    
    PCol <- NULL   # Matriz com os pesos das coluna nula
    
    ### INICIO - Centraliza na Media e Padroniza os dados por coluna,  ###
    ### assim teremos media zero e soma quadrado igual ao numero de linhas ###
    MC <- as.matrix(DataQ) # Matriz dados por grupo de variaveis
    
    if (sum(PondGeral)!=0) # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias 
      Media <- apply(sweep(MC,1,PondGeral,FUN="*"),2,sum) # Matriz com as medias por colunas poderada pelas linhas ponderadas da tabela de frequencia
    else
      Media <- apply(MC,2,mean) # Matriz com as medias por colunas
    
    MC <- sweep(MC, 2, Media, FUN = "-") # Centraliza na media
    
    if (sum(PondGeral)!=0) # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
      SqSum <- sqrt(colSums(sweep(as.matrix(MC^2),1,PondGeral,FUN="*")))  # raiz quadrada da soma ao quadrado dos elementos de MC dividido pelas linhas ponderadas da tabela de frequencia
    else
      SqSum <- sqrt(colSums(MC^2)/nrow(MC))
    
    MC <- sweep(MC, 2, SqSum, FUN = "/") # Normaliza os dados ou seja as somas dos quadrados e o numero de linhas  
    ### FIM - Centraliza na Media e Padroniza os dados  ###  
    
    if (sum(PondGeral)!=0) # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
      PLin <- PondGeral  # raiz quadrada da soma ao quadrado dos elementos de MC dividido pelas linhas ponderadas da tabela de frequencia
    else {
      MC <- as.matrix(DataQ)
      
      PLin <- rep(1,nrow(MC))
      
      NLin <- nrow(MC)
      
      SCol1 <- colSums(MC) / NLin
      
      MC <- sweep(MC, 2, SCol1, FUN = "-")
      
      SCol2 <- sqrt(colSums(MC^2)/NLin)
      
      MC <- sweep(MC, 2, SCol2, FUN = "/")
    }
    
    Pe <- GSVD(MC,PLin,rep(1,ncol(MC)))$d[1]^2   # Encontra o primeiro auto valor de MC
    
    PCol <- cbind(PCol,t(rep(1/Pe,ncol(MC)))) # Matriz com os pesos das colunas
    
    Lista <- list(MZ=MC, PLin=PLin, PCol=PCol)
    
    return(Lista)
  }
   
  MBC <- function(DataC,PondGeral) {  
    # Funcao que balanceia Dados Categoricos
    
    # Entrada:
    # DataQ - Dados a serem balanceados
    # PondGeral - usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
    
    # Retorna:
    # MZ   - Matriz Balanceada
    # PLin - Pesos das Linhas
    # PCol - Pesos das Colunas
    
    MZ   <- NULL   # cria uma matriz Z nula
    
    PLin <- NULL   # Matriz com os pesos das linhas nula
    
    PCol <- NULL   # Matriz com os pesos das coluna nula
    
    IM <- NULL     # Matriz Indicadora
    
    DB <- IM(DataC)  # Matriz dados binarios
    
    IM <- cbind(IM,DB) # Matriz Indicadora
    
    PVS <- CA_MFA(DB)  # Encontra o primeiro Valor Singular
    
    NL  <- nrow(DB)    # numero de linhas
    
    if (sum(PondGeral)!=0) # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias      
      PRL <- as.vector(PondGeral)  # pondera as linhas de acordo com os pesos das linhas da tabela de frequencia
    else  
      PRL <- as.vector(rep(1/NL,NL)) # probabilidade de ocorrencia de cada elemento da linha
    
    MB1 <- sweep(DB,1,PRL,FUN="*") # matriz pre-balanciada 1
    
    PVS <- CA_MFA(MB1) # Encontra o primeiro Valor Singular de MC
    
    SLI <- apply(MB1,2,sum)  # soma das colunas
    
    if (sum(PondGeral)!=0) # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
      MCO <- colSums(MB1)  # media das colunas
    else 
      MCO <- apply(MB1,2,mean) # media das colunas
    
    DIF <- 1-SLI             # 1 menos soma das colunas
    
    if (sum(PondGeral)!=0) # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
      MC <- sweep(DB,2,MCO,FUN="-") # matriz pre-balanciada 2 
    else 
      MC <- sweep(MB1,2,MCO,FUN="-") # matriz pre-balanciada 2 - subtrai MCO(media) de MB1
    
    if (sum(PondGeral)!=0) # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
      VET <- sqrt(colSums(sweep(as.matrix(MC^2),1,PRL,FUN="*")))  # raiz quadrada da soma ao quadrado dos elementos de MC dividido pelas linhas ponderadas da tabela de frequencia
    else
      VET <- sqrt(colSums(MC^2)/NL)  # raiz quadrada da soma ao quadrado dos elementos de MC dividido pelo numero de linhas
    
    MB  <- sweep(MC,2,VET,FUN="/") # matriz balanciada - divide MB2 por VET
    
    QVC <- sum(MB1)      # quantidade de categorias de variaveis
    
    Pe  <- DIF/(PVS*QVC) # valor usado para a ponderacao do PCA
    
    PCol <- rbind(PCol,as.matrix(Pe)) # Matriz com os pesos por coluna
    
    PLin <- rep(1/nrow(MB),nrow(MB)) # Matriz com os pesos das linhas  
    
    Lista <- list(MZ=MB, PLin=PLin, PCol=PCol)
    
    return(Lista)   
  }
  
  
  MBF <- function(DataF,SumLin) {  
    # Funcao que balanceia Dados de Frequencia
    
    # Entrada:
    # DataF  - Dados a serem balanceados
    # SumLin - Matriz com a soma geral das linhas dos dados de frequencia
    
    # Retorna:
    # MZ   - Matriz Balanceada
    # PLin - Pesos das Linhas
    # PCol - Pesos das Colunas
    
    FACI <- ICA(DataF,SumLin) # Retorna dados da funcao ICA - Analise de Correspondencia Interna
    
    MACI <- FACI$MatrixICA    # Matriz da Analise de Correspondencia Interna do Grupo i
    
    PCol <- NULL  # Matriz com os pesos das colunas 
    
    SLin <- as.matrix(FACI$MatrixSLin) # Matriz com as Somas das linhas da Matriz de Frequencia
    
    SCol <- as.matrix(FACI$MatrixSCol) # Matriz com as Somas das colunas da Matriz de Frequencia
    
    MPVS <- NULL  # Matriz com os primeiros Valores Singulares 
    
    PVS <- GSVD(MACI, SLin, SCol)$d[1]^2 # Encontra o primeiro Auto Valor do Grupo i
    
    MPVS <- rbind(MPVS,as.matrix(rep(PVS,ncol(MACI)))) # Matriz com os primeiros Valores Singulares 
    
    PCol <- sweep(SCol,1,MPVS,FUN="/") # Matriz com os Pesos das Colunas - divide cada soma das linhas pelos primeiros Auto Valores de cada grupo
    
    Lista <- list(MZ=MACI, PLin=SLin, PCol=PCol)
    
    return(Lista)  
  }
  
  if (is.null(NomeGrupos)) # Cria nomes para as variaveis caso nao exista
     NomeGrupos <- paste("Variavel", 1:length(TipoGrupo), sep = " ")
  
  if (!is.numeric(Grupo))
     stop("A entrada para 'Grupo' esta incorreta, deve ser do tipo numerico. Verifique!")
  
  if (!is.character(TipoGrupo))
     stop("A entrada para 'TipoGrupo' esta incorreta, deve ser do tipo caracter. Verifique!")
  
  if (!is.character(NomeGrupos))
     stop("A entrada para 'NomeGrupos' esta incorreta, deve ser do tipo caracter ou string. Verifique!")

  if (length(TipoGrupo)!=length(Grupo))
     stop("O numero de componetes da entrada 'TipoGrupo' difere da entrada 'Grupo'. Verifique!")

  if (length(NomeGrupos)!=length(Grupo))
     stop("O numero de componetes da entrada 'NomeGrupos' difere da entrada 'Grupo'. Verifique!")
  
  if (is.null(NomeGrupos)) # Cria nomes para as variaveis caso nao exista
     NomeGrupos <- paste("Variavel", 1:length(TipoGrupo), sep = " ")
  
  for (i in 1:length(TipoGrupo)) # transforma em maiusculo
    TipoGrupo[i]= ifelse(TipoGrupo[i]=="n","N",ifelse(TipoGrupo[i]=="c","C",ifelse(TipoGrupo[i]=="f","F",TipoGrupo[i])))   
    
  for (i in 1:length(TipoGrupo)) 
    if (TipoGrupo[i]!="N" && TipoGrupo[i]!="C" && TipoGrupo[i]!="F")
       stop("A entrada 'TipoGrupo' esta incorreta, deve ser: n, c, ou f. Verifique!")
  
  ### Inicio - Balanceia os valores dos grupos de variaveis ###
  NumGrupos = length(Grupo) # numero de grupos formados

  MZG   <- NULL  # cria uma matriz Geral Z nula
  
  PLinG <- NULL  # Matriz Geral com os pesos das linhas nula
  
  PColG <- NULL  # Matriz Geral com os pesos das coluna nula
  
  PondGeral <- 0 # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
  
  ### Inicio - Encontra as somas totais dos dados de frequencia ###
  if("F"%in%TipoGrupo) {
    
    SomaLinhas <- 0 # Matriz com a soma geral das linhas dos dados de frequencia
    
    j  <- 1        # coluna inicial do grupo de variaveis
    
    k  <- Grupo[1] # coluna final do grupo de variaveis
    
    for (i in 1:NumGrupos) {
      
      if (TipoGrupo[i]=="F") # Dados de Frequencia
         SomaLinhas <- SomaLinhas + apply(Data[,j:k],1,sum) # Matriz com a soma geral das linhas dos dados de frequencia
      
      j <- j + Grupo[i] # coluna inicial do grupo de variaveis
      
      k <- k + Grupo[i+ifelse(i!=NumGrupos,1,0)]  # coluna final do grupo de variaveis  
    }
    PondGeral <- SomaLinhas/sum(SomaLinhas) # usado para equilibrar os conjuntos quantitativos e categoricos, quando ha tabelas de frequencias
  } 
  ### Fim - Encontra as somas totais dos dados de frequencia ###
  
  j  <- 1       # coluna inicial do grupo de variaveis

  k  <- Grupo[1] # coluna final do grupo de variaveis
  
  for (i in 1:NumGrupos) {
      
     if (TipoGrupo[i]=="N"){  # Dados Quantitativos
        MB   <- MBQ(Data[,j:k],PondGeral)
        MZ   <- MB$MZ
        PLin <- MB$PLin
        PCol <- MB$PCol
        colnames(PCol) <- colnames(Data[,j:k])
     }
     
     if (TipoGrupo[i]=="C") { # Dados Categoricos
        MB   <- MBC(Data[,j:k],PondGeral)
        MZ   <- MB$MZ
        PLin <- MB$PLin
        PCol <- t(MB$PCol)
     }  
     
     if (TipoGrupo[i]=="F") {  # Dados de Frequencia
        MB   <- MBF(Data[,j:k],SomaLinhas)
        MZ   <- MB$MZ
        PLin <- t(MB$PLin)
        PCol <- t(MB$PCol)
     }  

     PLinG <- PLin  # Matriz Geral com os pesos das linhas
     
     PColG <- cbind(PColG,PCol) # Matriz Geral com os pesos das colunas
     
     MZG   <- cbind(MZG,MZ)     # Matriz Geral Balanceada
     
     j <- j + Grupo[i]    # coluna inicial do grupo de variaveis
      
     k <- k + Grupo[i+ifelse(i!=NumGrupos,1,0)]  # coluna final do grupo de variaveis  
     
     if (TipoGrupo[i]=="C")  # Dados Categoricos
        Grupo[i] <- ncol(MZ) # Como houve expansao da matriz de dados recoloca novo valor para o tamanho do Grupo

  }
  
  PColG <- t(PColG)  
  ### Fim - Balanceia os valores dos grupos de variaveis ###
  
  ### Inicio - Encontra os Autovetores e Autovalores ###
  MDS <- GSVD(MZG, PLinG, PColG) # Encontra a matriz de autovalor e autovetor
  MAutoVlr  <- MDS$d  # Matriz de Autovalores
  MAutoVecU <- MDS$u  # Matriz de Autovetores
  MAutoVecV <- MDS$v  # Matriz de Autovetores
  
  ## Matriz das Variancias
  MEigen <- as.data.frame(matrix(NA, length(MAutoVlr), 3))
  rownames(MEigen) <- paste("Comp", 1:length(MAutoVlr))
  colnames(MEigen) <- c("Autovalor", "% da variancia","% acumulada da variancia")
  MEigen[, "Autovalor"] <- MAutoVlr^2
  MEigen[, "% da variancia"] <- (MAutoVlr^2/sum(MAutoVlr^2)) * 100
  MEigen[, "% acumulada da variancia"] <- cumsum(MEigen[,"% da variancia"])
  
  NumAutoVlr <- length(MAutoVlr) # Numero de auto valores
  
  NE <- length(MAutoVlr[MAutoVlr>1e-10]) # Numero de elementos sigificativos dos Autovalores considerados somente valores acima de 10xe^(-9), isto e importante para calculos das inversas
  ### Fim - Encontra os Autovetores e Autovalores ###
  
  ### INICIO - Matriz Glogal Escore ###
  MF <-  MAutoVecU[,1:NE]%*%diag(MAutoVlr[1:NE],NE) # Matriz F - Matriz Glogal dos Escores de Fatores
  rownames(MF) <- rownames(Data) # Nomeia as linhas
  colnames(MF) <- paste("Comp.", 1:ncol(as.matrix(MF)), sep = " ") # Nomeia as colunas
  ### FIM - Matriz Glogal Escore ###
  
  ### INICIO - Matriz dos Escores dos Fatores por Grupo ###
  j  <- 1        # coluna inicial do Grupo de variaveis
  
  k  <- Grupo[1] # coluna final do Grupo de variaveis
  
  LMFGrupo <- as.list(1:NumGrupos) # cria lista vazia para a matriz de escores dos fatores dos Grupos
  
  for (i in 1:NumGrupos) {       
    
    MFG <- NumGrupos * MZG[,j:k]
    
    MFG <- sweep(MFG, 2, PColG[j:k], FUN="*")
 
    LMFGrupo[[i]] <- MFG%*%MAutoVecV[j:k,] # cria Matriz dos Escores dos Fatores por Grupo
    
    colnames(LMFGrupo[[i]]) <- paste("Comp.", 1:ncol(as.matrix(LMFGrupo[[i]])), sep = " ") # Nomeia as colunas
 
    j <- j + Grupo[i]      # coluna inicial do Grupo de variaveis
    
    k <- k + Grupo[i+ifelse(i!=NumGrupos,1,0)]  # coluna final do Grupo de variaveis  
  }
  
  names(LMFGrupo) <- paste("Grupo", 1:NumGrupos, sep = "") # nomeia os Grupos
  ### FIM - Matriz dos Escores dos Fatores por Grupo ###
  
  ### INICIO -  Correlacao dos Componentes Principais com as Variaveis Originais ###
  CCP <- sweep(as.matrix(MAutoVecV), 2, MAutoVlr, FUN = "*")  
  CCP <- t(CCP)
  rownames(CCP) <- paste("Comp.", 1:NumAutoVlr, sep = " ")
  colnames(CCP) <- colnames(MZG)
  ### FIM -  Correlacao dos Componentes Principais com as Variaveis Originais ###
  
  ### INICIO - Matriz das Inercias Parciais/Escores das Variareis ###
  CoordVar <- sweep(as.matrix(MAutoVecV), 2, sqrt(MAutoVlr), FUN = "*")  # Coordenadas das variaveis
  
  ContrVar <- sweep(as.matrix(CoordVar^2), 2, MAutoVlr, "/") # Contribuicao das variaveis
 
  ContrVar <- sweep(as.matrix(ContrVar), 1, PColG, "*")
  
  ContrGru <- matrix(data = NA, nrow = NumGrupos, ncol = NumAutoVlr) # Matriz com Contribuicoes dos Grupos
  
  j  <- 1        # coluna inicial do Grupo de variaveis
  
  k  <- Grupo[1] # coluna final do Grupo de variaveis
 
  for (i in 1:NumGrupos) {
    
    ContrGru[i,] <- apply(ContrVar[j:k, ], 2, sum) # Matriz com Contribuicoes dos Grupos
    
    j <- j + Grupo[i]      # coluna inicial do Grupo de variaveis
    
    k <- k + Grupo[i+ifelse(i!=NumGrupos,1,0)]  # coluna final do Grupo de variaveis  
    
  }
  
  EscVar <- sweep(ContrGru, 2, MAutoVlr^2, "*") # cria Matriz de Escores das variaveis/Inercias parciais
  
  colnames(EscVar) <- paste("Comp.", 1:ncol(as.matrix(EscVar)), sep = " ") # Nomeia as colunas
  
  rownames(EscVar) = NomeGrupos # Nomeias as linhas
  ### FIM - Matriz das Inercias Parciais/Escores das Variareis ###
    
  Lista <- list(MatrixG = Grupo, MatrixNG = NomeGrupos, MatrixPLin = PLinG,
                MatrixPCol = PColG, MatrixZ = MZG, MatrixA = MEigen,
                MatrixU = MAutoVecU, MatrixV = MAutoVecV, MatrixF = MF, 
                MatrixEFG = LMFGrupo, MatrixCCP = CCP, MatrixEscVar = EscVar)
  
  return(Lista)
}