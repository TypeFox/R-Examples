FA <- function(Data, Method = "PC", Type = 2, NFactor = 1, Rotation = "None", ScoresObs = "Bartlett", Screeplot = TRUE, Converg = 1e-5, Iteracao = 1000, TestFit = TRUE) {
   # Funcao executa a Analise Fatorial.
   # Desenvolvida por Paulo Cesar Ossani em 22/06/2013 e adapitada em 25/03/2016
   
   # Entrada:
   # Data      - Dados a serem analisados
   # Method    - Tipo de analises:
   #             Componentes Principais - PC (Principal Components) (default)
   #             Fator Principal - PF (Principal Factor)
   #             Maxima Verossimilhanca - ML (Maximum Likelihood)
   # Type      - 1 para analise utilizando a matriz de covariancia
   #             2 para analise utilizando a matriz de correlacao - default
   # Rotation  - Tipo de rotacao: "None" (default) e "Varimax" 
   # NFactor   - Numero de fatores (default = 1)
   # ScoresObs - Tipo de scores para as observacoes: "Bartlett" (default) ou "Regression"
   # Screeplot - Gera o grafico Screeplot para as variancias dos fatores (defaut = TRUE), somente para Rotation="None"
   # Converg   - Valor limite para convergencia para soma do quadrado dos residuos para metodo de Maxima Verossimilhanca (default = 1e-5)
   # Iteracao  - Numero maximo de iteracoes para metodo de Maxima Verossimilhanca (default = 1000)
   # TestFit   - Testa o ajusto do modelo para o metodo de Maxima Verossimilhanca (default = TRUE)
  
   # Saida:
   # MatrixMC      - Matriz de Correlacao/Covariancia
   # MatrixAutoVlr - Matriz de autovalores
   # MatrixAutoVec - Matriz de autovetores
   # MatrixVar     - Matriz de variancias e proporcoes
   # MatrixCarga   - Matriz de cargas fatoriais
   # MatrixVarEsp  - Matriz das variancias especificas
   # MatrixComuna  - Matriz das comunalidades
   # MatrixResiduo - Matriz dos residuos
   # VlrSQRS       - Valor limite superior para a soma do quadrados dos residuos
   # VlrSQR        - Soma dos Quadrados dos Residuos
   # MatrixResult  - Matriz com todos os resultados associados
   # MatrixScores  - Matriz com os escores das observarcoes

   if (!is.data.frame(Data)) 
      stop("Entrada 'Data' esta incorreta, deve ser do tipo dataframe. Verifique!")
  
   if (Method!="PC" && Method!="PF" && Method!="ML") 
      stop("Entrada 'Method' esta incorreta, deve ser 'PC', 'PF' ou 'ML'. Verifique!")
  
   if (Type!=1 && Type!=2) 
      stop("Entrada para 'Type' esta incorreta, deve ser numerica, sendo 1 ou 2. Verifique!")
  
   if (!is.numeric(NFactor)) 
      stop("Entrada para 'NFactor' esta incorreta, deve ser numerica. Verifique!")

   if (NFactor>ncol(Data)) 
      stop("Entrada para 'NFactor' esta incorreta, deve ser igual ou inferior ao numero de variaveis de 'Data'. Verifique!")
 
   if (NFactor<=0) 
      stop("Entrada para 'NFactor' esta incorreta, deve ser numero inteiro maior ou igual a 1. Verifique!")
 
   if (Rotation!="None" && Rotation!="Varimax" && Rotation!="Quartimax") 
      stop("Entrada para 'Rotation' esta incorreta, deve ser 'None', 'Varimax', 'Quartimax' ou ?????. Verifique!")
  
   if (Rotation!="None" && NFactor<2)
      stop("Para a rotacao, he necessario mais do que um fator. Altere o numero de fatores (NFactor) para continuar!")

   if (ScoresObs!="Bartlett" && ScoresObs!="Regression") 
      stop("Entrada para 'ScoresObs' esta incorreta, deve ser 'Bartlett' ou 'Regression'. Verifique!")
   
   if (!is.logical(Screeplot))
      stop("Entrada para 'Screeplot' esta incorreta, deve ser TRUE ou FALSE. Verifique!")
       
   if (Type == 1)     # Considera a Matriz de Covariancia para a decomposicao
      MC <- cov(Data) # Matriz de Covariancia
  
   if (Type == 2)     # Considera a Matriz de Correlacao para a decomposicao
      MC <- cor(Data) # Matriz de Correlacao
  
   #library("stats")  
   #library("MASS")

   Rotacao <- function(MData,Type=NULL,Normalise=TRUE) {
   # Funcao que executa as rotacoes
     if (Type=="Varimax") {
        Var <- varimax(MData,normalize=Normalise)
        Res <- Var$loadings[,]
     }
     
     return(Res)
   }
   
   if (Method=="PC") { # Metodo dos Componentes Principais
      
      # Encontrando a Matriz de Decomposicao Expectral
      MAV <- eigen(MC) # Encontra a matriz de autovalor e autovetor
      MAutoVlr <- MAV$values  # Matriz de Autovalores 
      MAutoVec <- MAV$vectors # Matriz de Autovetores
  
      Gama = MAutoVec%*%diag(sqrt(abs(MAutoVlr)),nrow(MC),ncol(MC)) # Matriz de Cargas Fatoriais
      if (Rotation!="None") {
        Gama <- Rotacao(Gama[,1:NFactor],Rotation)
        MAutoVlr <- colSums(Gama^2)
      }
      rownames(Gama) <- colnames(Data)
      colnames(Gama) <- paste("Fator",1:ncol(Gama))
      
      Psi = diag(MC - Gama[,1:NFactor]%*%t(Gama[,1:NFactor])) # Matriz de Variancias Especificas
      
      Comun = diag(MC - Psi) # Matriz de Comunalidades
 
      # Valor Limite Superior para a Soma de Quadrados de Residuos
      SQRS = MAutoVlr[(NFactor+1):nrow(as.matrix(MAutoVlr))]%*%(MAutoVlr[(NFactor+1):nrow(as.matrix(MAutoVlr))])
     
      M = MC - (Gama[,1:NFactor]%*%t(Gama[,1:NFactor]) + diag(Psi)) # Matriz dos residuos 
      
      SQR = sum(diag(M%*%t(M))) # Soma dos Quadrados dos Residuos
      
      # Matriz das Variancias
      MEigen <- as.data.frame(matrix(NA, length(MAutoVlr), 3))
      rownames(MEigen) <- paste("Fator", 1:length(MAutoVlr))
      colnames(MEigen) <- c("Autovalor", "% da variancia","% acumulada da variancia")
      MEigen[, "Autovalor"] <- MAutoVlr
      MEigen[, "% da variancia"] <- (MAutoVlr/sum(MAutoVlr)) * 100
      MEigen[, "% acumulada da variancia"] <- cumsum(MEigen[,"% da variancia"]) 
      
      # Matriz com todos os resultados associados
      Result <- as.matrix(cbind(Gama[,1:NFactor],Comun,Psi))
      Result <- rbind(Result,t(rbind(as.matrix(MEigen[1:NFactor,1]),sum(Comun),NA)))
      Result <- rbind(Result,t(rbind(as.matrix(MEigen[1:NFactor,2]/100),MEigen[NFactor,3]/100,NA)))
      colnames(Result) <- c(paste("Carga Fator",1:NFactor),"Comunalidade","Variancias especificas")
      rownames(Result) <- c(colnames(Data),"Variancia","% Variancia")
      
   }
    
   if (Method=="PF") { # Metodo dos Fatores Principais
     
      Psi0 <- (solve(diag(diag(solve(MC))))) # Encontrando a Matriz Psi

      Sr <- MC - Psi0 # Encontrando a Matriz Sr

      # Encontrando a Matriz de Decomposicao Expectral
      MAV <- eigen(Sr) # Encontra a matriz de autovalor e autovetor
      MAutoVlr <- MAV$values  # Matriz de Autovalores 
      MAutoVec <- MAV$vectors # Matriz de Autovetores

      Gama = MAutoVec%*%diag(sqrt(abs(MAutoVlr)),nrow(MC),ncol(MC)) # Matriz de Cargas Fatoriais
      if (Rotation!="None") {
        Gama <- Rotacao(Gama[,1:NFactor],Rotation)
        MAutoVlr <- colSums(Gama^2)
      }
      rownames(Gama) <- colnames(Data)
      colnames(Gama) <- paste("Fator",1:ncol(Gama))
      
      Psi = diag(MC - Gama[,1:NFactor]%*%t(Gama[,1:NFactor])) # Matriz de Variancias Especificas
     
      Comun = diag(MC - Psi) # Matriz de Comunalidades
      
      ## Valor Limite Superior para a Soma de Quadrados de Residuos
      SQRS = MAutoVlr[(NFactor+1):nrow(as.matrix(MAutoVlr))]%*%(MAutoVlr[(NFactor+1):nrow(as.matrix(MAutoVlr))])
   
      # Soma dos Quadrados dos Residuos
      M = MC - (Gama[,1:NFactor]%*%t(Gama[,1:NFactor]) + diag(Psi))
      SQR = sum(diag(M%*%t(M)))
      if (Rotation!="None") {
        Gama <- Rotacao(Gama[,1:NFactor],Rotation)
        MAutoVlr <- colSums(Gama^2)
      }
      rownames(Gama) <- colnames(Data)
      colnames(Gama) <- paste("Fator",1:ncol(Gama))
      
      # Matriz das Variancias
      MEigen <- as.data.frame(matrix(NA, length(MAutoVlr), 3))
      rownames(MEigen) <- paste("Comp", 1:length(MAutoVlr))
      colnames(MEigen) <- c("Autovalor", "% da variancia","% acumulada da variancia")
      MEigen[, "Autovalor"] <- MAutoVlr
      MEigen[, "% da variancia"] <- (MAutoVlr/sum(MAutoVlr)) * 100
      MEigen[, "% acumulada da variancia"] <- cumsum(MEigen[,"% da variancia"])
      
      # Matriz com todos os resultados associados
      Result <- as.matrix(cbind(Gama[,1:NFactor],Comun,Psi))
      Result <- rbind(Result,t(rbind(as.matrix(MEigen[1:NFactor,1]),sum(Comun),NA)))
      Result <- rbind(Result,t(rbind(as.matrix(MEigen[1:NFactor,2]/100),MEigen[NFactor,3]/100,NA)))
      colnames(Result) <- c(paste("Carga Fator",1:NFactor),"Comunalidade","Variancias especificas")
      rownames(Result) <- c(colnames(Data),"Variancia","% Variancia")
   }
   
   if (Method=="ML") { # Metodo de maxima verossimilhanca
   
      n <- ncol(Data)*nrow(Data) # numero de elementos amostrais
      MC <- (n-ncol(Data))/n*MC  # Matriz de Covariancia/Correlacao Maximizada para o teste
      
      # Encontrando a Matriz de Decomposicao Expectral
      MAV <- eigen(MC) # Encontra a matriz de autovalor e autovetor
      MAutoVlr <- MAV$values  # Matriz de Autovalores 
      MAutoVec <- MAV$vectors # Matriz de Autovetores

      Gama = MAutoVec%*%diag(sqrt(abs(MAutoVlr)),nrow(MC),ncol(MC)) # Matriz de Cargas Fatoriais para Inicializacao da Iteracao

      Psi = (diag(MC - Gama[,1:NFactor]%*%t(Gama[,1:NFactor]))) # Matriz das Variancias Especificas
   
      M = MC - (Gama[,1:NFactor]%*%t(Gama[,1:NFactor]) + diag(Psi)) # Matriz dos residuos
      
      SQRi= sum(diag(M%*%t(M))) # Soma dos Quadrados dos Residuos

      ### INICIO DA ITERAcaO ###
      i = 1 # inicializa o contador de iteracoes
      while (1) {
         MC_new = diag(1/sqrt(Psi))%*%(MC - diag(Psi))%*% diag(1/sqrt(Psi)) # nova matriz para iteracao
   
         # Encontrando a Matriz de Decomposicao Expectral
         MAV <- eigen(MC_new) # Encontra a matriz de autovalor e autovetor
         MAutoVlr1 <- MAV$values  # Matriz de Autovalores 
         MAutoVec1 <- MAV$vectors # Matriz de Autovetores
         
         # Matriz das Cargas Fatoriais
         Gama_new = diag(sqrt(Psi))%*%MAutoVec1%*%diag(sqrt(abs(MAutoVlr1)),nrow(MC_new),ncol(MC_new))
   
         Psi = (diag(MC - Gama_new[,1:NFactor]%*%t(Gama_new[,1:NFactor]))) # Matriz das Variancias Especificas
   
         # Valor Limite Superior para a Soma de Quadrados de Residuos
         SQRS = MAutoVlr1[(NFactor+1):nrow(as.matrix(MAutoVlr1))]%*%(MAutoVlr[(NFactor+1):nrow(as.matrix(MAutoVlr1))])
         
         M = MC - (Gama_new[,1:NFactor]%*%t(Gama_new[,1:NFactor]) + diag(Psi)) # Matriz dos Residuos
         
         SQR = sum(diag(M%*%t(M))) # Soma dos Quadrados dos Residuos
        
         if (SQR <= Converg) break # sai do loop quando atingir a convergencia
       
         if (i >= Iteracao) break # sai do loop apos esse limite de iteracoes
       
         i = i + 1 # incrementa o contador de iteracoes
         
      }
      ### FIM DA ITERAcaO ###
      
      Gama = Gama_new # Matriz com as cargas fatoriais
  
      if (Rotation!="None") 
         Gama <- Rotacao(Gama[,1:NFactor],Rotation,Normalise=TRUE)

      rownames(Gama) <- colnames(Data)
      colnames(Gama) <- paste("Fator",1:ncol(Gama))
      
      if (Type == 1) {# Considera a Matriz de Covariancia para a decomposicao
         Gama <- diag(1/sqrt(diag(MC)))%*%Gama[,1:NFactor] # Matriz com as cargas fatoriais
         Comun = rowSums(Gama^2)#apply(Gama,1,function(Gama) Gama^2)) # Matriz de Comunalidades
      }
      
      MAutoVlr <- colSums(Gama^2)
      
      if (Type == 2)     # Considera a Matriz de Correlacao para a decomposicao
         Comun = diag(MC - Psi) # Matriz de Comunalidades
      
      # Matriz das Variancias
      MEigen <- as.data.frame(matrix(NA, length(MAutoVlr), 3))
      rownames(MEigen) <- paste("Comp", 1:length(MAutoVlr))
      colnames(MEigen) <- c("Autovalor", "% da variancia","% acumulada da variancia")
      MEigen[, "Autovalor"] <- MAutoVlr
      MEigen[, "% da variancia"] <- (MAutoVlr/sum(MAutoVlr)) * 100
      MEigen[, "% acumulada da variancia"] <- cumsum(MEigen[,"% da variancia"])
      
      print(paste("Numero de iteracoes:",i))
      
      # Matriz com todos os resultados associados
      Result <- as.matrix(cbind(Gama[,1:NFactor],Comun,Psi))
      Result <- rbind(Result,t(rbind(as.matrix(MEigen[1:NFactor,1]),sum(Comun),NA)))
      Result <- rbind(Result,t(rbind(as.matrix(MEigen[1:NFactor,2]/100),MEigen[NFactor,3]/100,NA)))
      colnames(Result) <- c(paste("Carga Fator",1:NFactor),"Comunalidade","Variancias especificas")
      rownames(Result) <- c(colnames(Data),"Variancia","% Variancia")  
      
      ### INICIO - Teste da falta de ajusto do modelo fatorial - teste Qui-quadrado ###
      if (TestFit) {
         p <- nrow(Gama)  # numero de parametros
      
         gl <- ((p - NFactor)^2 - NFactor - p)/2 # grau de liberdade
    
         cat("### TESTE DO AJUSTO DO MODELO ###\n")
      
         cat(paste("Grau de liberdade observado:", round(gl,5)),"\n")
        
         if (gl < 0) 
            cat("Nao foi possivel realizar o teste de ajuste do modelo, pois grau de libertade foi negativo, aconselha-se a mudar os parametros, para processeguir com o teste. Exemplo: numero de fatores ou mesmo 'Type'.\n")
   
         if (det(MC)<=0) 
            cat("Nao foi possivel realizar o teste de ajuste do modelo, pois o determinante da matriz de variancia/covariancia deve ser diferente de zero, para processeguir com o teste mude os parametros.\n")
         
         if (gl>=0 && det(MC)>0) {
          
            Ps_i = diag(diag(MC - Gama[,1:NFactor]%*%t(Gama[,1:NFactor])))
          
            Chi.Quad.Observado <- (n - 1 - (2*p + 5)/6 - 2*NFactor/3)*log(det(Gama[,1:NFactor]%*%t(Gama[,1:NFactor])+Ps_i)/det(MC))

            qt = qchisq(0.95,gl,ncp=0)
    
            cat(paste("Valor da estatistica do teste Qui-quadrado (Chiq1):", round(Chi.Quad.Observado,3)),"\n")
         
            cat(paste("Valor Qui-quadrado observado (Chiq2) com 5% de significancia:", round(qt,3)),"\n")
        
            if (Chi.Quad.Observado<=qt) cat("Como Chiq1 <= Chiq2, verifica-se que o numero de fatores FORAM suficientes.\n")
        
            if (Chi.Quad.Observado>qt) cat("Como Chiq1 > Chiq2, verifica-se que o numero de fatores NAO FORAM suficientes.\n")
            
            cat("Valor-p:", pchisq(Chi.Quad.Observado,gl,ncp=0, lower.tail = F))
         } 
      }
      ### FIM - Teste da falta de ajusto do modelo fatorial - teste Qui-quadrado ###
   }
   
   ### INICIO - Scree-plot dos fatores ####
   if (Screeplot && Rotation=="None")
      plot(1:length(MEigen[,1]), MEigen[,1], type = "b", 
           xlab = "Ordem dos fatores", 
           ylab = "Variancia dos fatores",
           main = "Scree-plot das variancias dos fatores sem rotacao")
   ### FIM - Scree-plot dos fatores
   
   ### INICIO - encontrar os scores das observacoes ###
   if (Type == 1)  {   # Considera a Matriz de Covariancia para os calculos
      Media  <- apply(Data, 2, mean)
      DataNorm <- sweep(as.matrix(Data), 2, Media, FUN = "-") # Centraliza na media por colunas
   }
   
   if (Type == 2) { # Considera a Matriz de Correlacao para os calculos
      # Centraliza na media por colunas e divide pelo desvio padrao de cada coluna
      Media  <- apply(Data, 2, mean) # DataNorm com as medias por colunas
      DataNorm <- sweep(Data, 2, Media, FUN = "-")   # Centraliza na media
      Desvio <- sqrt(colSums(DataNorm^2)/(nrow(DataNorm)-1)) # raiz da soma do quadrado - desvio padrao amostral
      DataNorm <- sweep(DataNorm, 2, Desvio, FUN = "/")  # Divide pelo desvio padrao
   }
   
   if (ScoresObs=="Bartlett") { # Metodo Bartlett (minimos quadrados)
      # foi necessario usar a inversa generalizada pois algumas vezes a matriz he singular, assim nao tem inversa normal
      Scores <- MASS::ginv(t(Gama)%*%solve(diag(Psi))%*%Gama)%*%(t(Gama)%*%solve(diag(Psi)))%*%t(DataNorm) # Matriz com os escores das observacoes
      #Scores <- solve(t(Gama)%*%solve(diag(Psi))%*%Gama)%*%(t(Gama)%*%solve(diag(Psi)))%*%t(DataNorm) # Matriz com os escores das observacoes
      #Scores <- DataNorm%*%solve(MC)%*%Gama # outro modo de encontrar a solucao acima
   }
   
   if (ScoresObs=="Regression") { # Metodo de Regressao
      Media <- mean(as.matrix(Data))
      DataNorm <- sweep(as.matrix(Data), 2, Media, FUN = "-") # Centraliza na media geral todos os dados
      I <- diag(rep(ncol(Gama)))
      Scores <- solve(I + t(Gama)%*%solve(diag(Psi))%*%Gama)%*%(t(Gama)%*%solve(diag(Psi)))%*%t(DataNorm) # Matriz com os escores das observacoes
      #Scores <- t(Gama)%*%solve(Gama%*%t(Gama)+diag(Psi))%*%t(DataNorm) # outro modo de encontrar a solucao acima
   }
   Scores <- t(Scores)
   colnames(Scores) <- colnames(Gama)
   rownames(Scores) <- rownames(Data)
   ### FIM - encontrar os scores das observacoes ###  
   
   Lista <- list(MatrixMC = MC, MatrixAutoVlr = MAutoVlr,
                 MatrixAutoVec = MAutoVec, MatrixVar = MEigen,
                 MatrixCarga = Gama[,1:NFactor], MatrixVarEsp = Psi,
                 MatrixComuna = Comun, MatrixResiduo = M, VlrSQRS = SQRS,
                 VlrSQR = SQR, MatrixResult = Result, MatrixScores=Scores[,1:NFactor])

   return(Lista)
}
