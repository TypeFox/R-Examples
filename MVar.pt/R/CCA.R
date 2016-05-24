CCA <- function(X = NULL, Y = NULL, Type = 1, Test = "Bartlett", Sign = 0.05) {
   # Esta funcao executa a Analise de Correlacao Canonica
   # Desenvolvida por Paulo Cesar Ossani em 07/2013
  
   # Entrada:
   # X - Primeiro grupo de variaveis de um conjunto dados.
   # Y - Segundo grupo de variaveis de um conjunto dados.
   # Type  - 1 para analise utilizando a matriz de covariancia (default),
   #         2 para analise utilizando a matriz de correlacao,
   # Test  - Teste significancia da relacao entre o grupo X e Y:
   #        "Bartlett" (default) ou "Rao".
   # Sign - Grau de significancia do teste (default 5%)

   # Retorna:
   # Cxx     - Matriz de Covariancia ou Correlacao Cxx,
   # Cyy     - Matriz de Covariancia ou Correlacao Cyy,
   # Cxy     - Matriz de Covariancia ou Correlacao Cxy,
   # Cyx     - Matriz de Covariancia ou Correlacao Cyx,
   # Var.UV  - Matriz com autovalores (variancias) dos pares cononicos U e V,
   # Corr.UV - Matriz de Correlacao dos pares cononicos U e V,
   # Coef.X  - Matriz dos Coeficientes canonicos do grupo X,
   # Coef.Y  - Matriz dos Coeficientes canonicos do grupo Y,
   # Coor.X  - Matriz das Correlacoes entre as variaveis canonicas e as variaveis originais do grupo X,
   # Coor.Y  - Matriz das Correlacoes entre as variaveis canonicas e as variaveis originais do grupo Y,
   # Score.X - Matriz com os scores do grupo X,
   # Score.Y - Matriz com os scores do grupo Y,
   # SigTest - Restorna o teste significancia da relacao entre o grupo X e Y, 'Bartlett' (default) ou "Rao".
  
   if (!is.data.frame(X)) 
      stop("Entrada 'X' esta incorreta, deve ser do tipo dataframe. Verifique!")

   if (!is.data.frame(Y)) 
      stop("Entrada 'Y' esta incorreta, deve ser do tipo dataframe. Verifique!")
  
   if (nrow(X)!=nrow(Y)) 
      stop("O numero de observacoes de 'X' deve ser igual a de 'Y'. Verifique!")
  
   if (Type!=1 && Type!=2) 
      stop("Entrada para 'Type' esta incorreta. Verifique!")
  
   if (Test!="Bartlett" && Test!="Rao")
      stop("Entrada para 'Test' esta incorreta, deve ser 'Bartlett' ou 'Rao'. Verifique!")
   
   if (!is.numeric(Sign)) 
      stop("Entrada para 'Sign' esta incorreta, deve ser numerica com valores entre 0 e 1. Verifique!")
  
   if (Sign<=0 || Sign>1) 
      stop("Entrada para 'Sign' esta incorreta, deve ser valores entre 0 e 1. Verifique!")
  
   if (Type == 1) { # Considera a Matriz de Covariancia para a decomposicao
      MC  = cov(cbind(X,Y))  # Matriz de Covariancia
      S11 = cov(X)   # Matriz de Covariancia dos primeiros dados
      S22 = cov(Y)   # Matriz de Covariancia dos segundos dados
      S12 = cov(X,Y) # Matriz de Covariancia de X e Y
      S21 = cov(Y,X) # Matriz de Covariancia de Y e X
   }
    
   if (Type == 2) { # Considera a Matriz de Correlacao para a decomposicao
      MC  = cor(cbind(X,Y)) ## Matriz de Correlacao
      S11 = cor(X)   # Matriz dos Coeficientes dos primeiros dados
      S22 = cor(Y)   # Matriz dos Coeficientes dos segundos dados
      S12 = cor(X,Y) # Matriz dos Coeficientes de X e Y
      S21 = cor(Y,X) # Matriz dos Coeficientes de Y e X
   }

   ### Calculo da Matriz S11^(-1/2) e inversa de S22 
   M1 <- eigen(S11) # Calcula a matriz de autovalor e de autovetor de S11
   MAutoVlr1 <- M1$values  # Matriz de Autovalores 
   MAutoVec1 <- M1$vectors # Matriz de Autovetores
   S11_12 = MAutoVec1%*%diag(1/sqrt(MAutoVlr1))%*%t(MAutoVec1)
   S22_Inv = solve(S22) # Calcula a inversa da matriz S22

   ### Calcula a matriz S11_12 * S12 * S22_Inv * S21 * S11_12
   M_aux = S11_12%*%S12%*%S22_Inv%*%S21%*%S11_12
   M2 <- eigen(M_aux) # Calcula a matriz de autovalor e de autovetor de S11
   MAutoVlr2 <- M2$values  # Matriz de Autovalores 
   MAutoVec2 <- M2$vectors # Matriz de Autovetores

   ### Matriz das variancias dos pares canonicos
   Mr = MAutoVlr2 # Matriz r das variancias dos pares canonicos
   MEigen <- as.data.frame(matrix(NA, length(Mr), 3))
   rownames(MEigen) <- paste("U",1:length(Mr),"V", 1:length(Mr),sep="")
   colnames(MEigen) <- c("Variancia", "% da Variancia","% acumulada da Variancia")
   MEigen[, "Variancia"] <- Mr
   MEigen[, "% da Variancia"] <- (Mr/sum(Mr)) * 100
   MEigen[, "% acumulada da Variancia"] <- cumsum(MEigen[,"% da Variancia"])
   
   ## Matriz de correlacao dos pares canonicos
   CorrUV <- as.matrix(sqrt(Mr),ncol=length(CorrUV),nrow=1)
   rownames(CorrUV) <- paste("U",1:length(CorrUV),"V", 1:length(CorrUV),sep="")
   colnames(CorrUV) <- c("Correlacao")
   
   ### Matriz dos coeficientes canonicos do primeiro conjunto
   Coef_A = S11_12%*%MAutoVec2
   rownames(Coef_A) <- colnames(X)
   colnames(Coef_A) <- paste("U",1:ncol(Coef_A), sep="")
  
   ### Matriz dos coeficientes canonicos do segundo conjunto
   Coef_B = S22_Inv%*%S21%*%Coef_A%*%solve(diag(sqrt(MAutoVlr2)))
   colnames(Coef_B) <- paste("V",1:ncol(Coef_B), sep="")
   
   ### INICIO - Determina as Correlacoes amostrais entre as variaveis ###
   ### canonicas e as variaveis originais de cada grupo ###
   ## Calculo da Matriz D11^(-1/2) 
   M3 <- eigen(diag(diag(S11))) # Calcula a matriz de autovalor e de autovetor de D11^(-1/2)
   MAutoVlr3 <- M3$values  # Matriz de Autovalores 
   MAutoVec3 <- M3$vectors # Matriz de Autovetores
   D11_12 = MAutoVec3%*%diag(1/sqrt(MAutoVlr3))%*%t(MAutoVec3)
   
   ## Calculo da Matriz D22^(-1/2) 
   M4 <- eigen(diag(diag(S22))) # Calcula a matriz de autovalor e de autovetor de D22^(-1/2)
   MAutoVlr4 <- M4$values  # Matriz de Autovalores 
   MAutoVec4 <- M4$vectors # Matriz de Autovetores
   
   D22_12 = MAutoVec4%*%diag(1/sqrt(MAutoVlr4))%*%t(MAutoVec4)
   
   Rux = t(Coef_A)%*%S11%*%D11_12 # Matriz de correlacao U X
   Ruy = t(Coef_A)%*%S12%*%D22_12 # Matriz de correlacao U Y
   Rvy = t(Coef_B)%*%S22%*%D22_12 # Matriz de correlacao V Y
   Rvx = t(Coef_B)%*%S21%*%D11_12 # Matriz de correlacao V X

   Rux <- t(Rux) # Grupo X para U
   rownames(Rux) <- colnames(X)
   Ruy <- t(Ruy) # Grupo Y para U
   rownames(Ruy) <- colnames(Y)
   Rvx <- t(Rvx) # Grupo X para V
   rownames(Rvx) <- colnames(X)
   Rvy <- t(Rvy) # Grupo Y para V
   rownames(Rvy) <- colnames(Y) 

   CVU <- cbind(Rux,Rvx)
   CVV <- cbind(Ruy,Rvy)
   ### FIM - Determina as Correlacoes amostrais entre as variaveis ###
   ### canonicas e as variaveis originais de cada grupo ###
   
   ### INICIO - Calculo dos scores dos grupos ###
   X.Aux = scale(X, center = TRUE, scale = FALSE) # Centraliza na media
   Y.Aux = scale(Y, center = TRUE, scale = FALSE) # Centraliza na media
   X.Aux[is.na(X.Aux)] = 0
   Y.Aux[is.na(Y.Aux)] = 0
   X.Scores = X.Aux %*% Coef_A # Scores co grupo X
   Y.Scores = Y.Aux %*% Coef_B # Scores co grupo Y
   ### INICIO - Calculo dos scores dos grupos ###
   
   ### INICIO - Teste Qui-Quadrado de Bartlett ###
   ## Este teste he usado para saber quantos pares de variaveis canonica sao significativos
   if (Test=="Bartlett") {
      n <- nrow(X) # numero de observacoes
      p <- ncol(X) # numero de variaveis de X
      q <- ncol(Y) # numero de variaveis de Y 
      QtdF <- length(CorrUV) # quantidade de fatores
     
      T.Bartlett <- as.data.frame(matrix(NA, QtdF, 6))
      colnames(T.Bartlett) <- c("Pares Canonicos", "Lambda de Wilks","Qui-Quadrado aproximado","Densidade Qui-Quadrado","Grau de Liberdade","Valor-p")
      T.Bartlett[,1]<-paste("U",1: QtdF,"V", 1: QtdF,sep="")
      i=1
      for (i in 1:QtdF) {
          Lambda <- prod(1-CorrUV[i:QtdF]^2) # lambida de Wilks
          
          Chi.Observado <- -((n - 1) - (p + q + 1)/2)*log(Lambda) # Estatistica do Teste
          
          gl  <- (p -i +1)*(q -i +1) # grau de libardade
          
          Chi.Encontrado <- qchisq(1 - Sign,gl,ncp=0)
          
          pValor <- pchisq(Chi.Observado,gl,ncp=0, lower.tail = F) 
          
          T.Bartlett[i, "Lambda de Wilks"] <- round(Lambda,5)
          T.Bartlett[i, "Qui-Quadrado aproximado"] <- round(Chi.Observado,5)
          T.Bartlett[i, "Densidade Qui-Quadrado"] <- round(Chi.Encontrado,5)
          T.Bartlett[i, "Grau de Liberdade"] <- gl
          T.Bartlett[i, "Valor-p"] <- round(pValor,5)
      }
      Teste <- T.Bartlett
   }
   ### FIM - Teste Qui-Quadrado de Bartlett ###
   
   ### INICIO - Teste F de RAO ###
   ## Este teste he usado para saber quantos pares de variaveis canonica sao significativos
   if(Test=="Rao") {
     n  <- nrow(X) # numero de observacoes
     p1 <- ncol(X) # numero de variaveis de X
     q1 <- ncol(Y) # numero de variaveis de Y 
     QtdF <- length(CorrUV) # quantidade de fatores
   
     T.Rao <- as.data.frame(matrix(NA, QtdF, 7))
     colnames(T.Rao) <- c("Pares Canonicos", "Lambda de Wilks","F aproximado","Densidade F","Grau de Liberdade 1","Grau de Liberdade 2","Valor-p")
     T.Rao[,1]<-paste("U",1: QtdF,"V", 1: QtdF,sep="")
     
     for (i in 1:QtdF) {
       p <- p1 - i + 1
       q <- q1 - i + 1
  
       t <- (n - 1) - (p + q + 1)/2
       
       s <- ifelse((p^2 + q^2) <= 5, 1, sqrt((p^2 * q^2 -4)/(p^2 + q^2 -5)))
       
       Lambda <- prod(1-CorrUV[i:QtdF]^2) # lambida de Wilks
       
       gl1 <- p * q # grau de liberdade 1
       gl2 <- (1 + t*s - p*q/2) # grau de liberdade 2
       F.Observado <-((1 - Lambda^(1/s))/Lambda^(1/s)) * gl2/gl1 # Estatistica do Teste
       
       F.Encontrado <- qf(1-Sign,gl1,gl2,ncp=0)
       
       pValor <- pf(F.Observado,gl1,gl2,ncp=0, lower.tail = FALSE)
       
       T.Rao[i, "Lambda de Wilks"] <- round(Lambda,5)
       T.Rao[i, "F aproximado"] <- round(F.Observado,5)
       T.Rao[i, "Densidade F"] <- round(F.Encontrado,5)
       T.Rao[i, "Grau de Liberdade 1"] <- gl1
       T.Rao[i, "Grau de Liberdade 2"] <- round(gl2,5)
       T.Rao[i, "Valor-p"] <- round(pValor,5)
     }
     Teste <- T.Rao
   }
   ### FIM - Teste F de RAO ###

   Lista <- list(MC = MC, Cxx = S11, Cyy = S22, Cxy = S12, 
                 Cyx = S21, Var.UV = MEigen, Corr.UV =  CorrUV, 
                 Coef.X = Coef_A, Coef.Y = Coef_B, Coor.X = CVU, 
                 Coor.Y = CVV, Score.X = X.Scores, Score.Y = Y.Scores,
                 SigTest = Teste)
   
   return(Lista)
}
