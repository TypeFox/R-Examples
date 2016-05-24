ccF <-
function(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE,
    main = NULL){
    trans.fator = function(matriz){
      trt <- matriz[,1]
      y <- matriz[,2]
      list(trt = trt,y = y)
    }
    anova.dados = function(dados){  
      dados$trt <- as.factor(dados$trt)  
      anovadados <- lm(dados$y ~ dados$trt)
      anovadad <- anova(anovadados)
      anovadad
    }
    anova.medias = function(dados){
      medias <- tapply(dados$y,list(dados$trt),mean)
    }
    sq.f = function(medias,r,k){
      S <- matrix(0,k,k)
      for (i in 1:(k-1)){
        seq <- (i+1):k
        for (j in seq){
          g <- medias[i:j]
          S[i,j] <- r*(j-i)*var(g)
        }
      }
      S  
    }
    valores.anava = function(dados,r){
      anava <- anova.dados(dados)
      S2 <- anava$"Mean Sq"[2]
      f <- anava$"Df"[2]
      raiz <- sqrt(S2/r)
      list(S2 = S2,f = f,raiz = raiz)
    }
    particiona = function(med,SQ.m,S2,f){
      result <- matrix(1,k,k+2)
      result[1:k,k+1] <- 0
      result[k,1:k] <- 1:k
      p <- 1              
      result[1,k+1] <- SQ.m[1,k]
      C <- SQ.m[1,k]/(S2*(k-1))
      valorp <- 1 - pf(C,k-1,f)
      result[p,k+2] <- valorp
      for(p in 2:(k-1)){
        kk <- kmeans(med,p,algorithm="Hartigan-Wong",nstart=10)
        particao <- as.matrix(kk$cluster)
        result[p,1:k] <- particao
        soma <- 0
        aux <- 1 
        for(i in 1:k){
          if(result[p,i] != result[p,i+1]){
            soma <- soma + SQ.m[aux,i]
            aux <- i + 1
          }
        }
        result[p,k+1] <- soma
        C <- result[p,k+1]/(S2*(k-1))
        valorp <- 1 - pf(C,k-1,f)
        result[p,k+2] <- valorp
      }
      result
    }
    anova.dados = function(dados){  
      dados$trt <- as.factor(dados$trt)  
      anovadados <- lm(dados$y ~ dados$trt)
      anovadad <- anova(anovadados)
      anovadad
    }
    k <- length(tapply(y, trt, length))
    r <- length(trt)/k
    dados <- cbind(trt,y)
    dados <- trans.fator(dados)
    medias <- anova.medias(dados)
    med.ord <- sort(medias)
    S2f <- cbind(SSerror/DFerror,DFerror)
    valores <- as.vector(S2f)
    names(valores) <- c("S2","f")
    nomes <- as.numeric(rownames(med.ord))
    somaq <- sq.f(medias,r,k)
    cluster <- particiona(med.ord,somaq,valores[1],valores[2])
    grupos <- vector(mode = "integer", length = k)
    i <- 1
    while(cluster[i,k+2] < alpha){
      grupos <- cluster[i+1,1:k]
      i <- i + 1
    }
    grupos   
    agrupamento <- matrix(0,k,2)
    agrupamento[,1] <- nomes
    agrupamento[,2] <- grupos
    j = k
    group <- matrix(0,k,2)
    for(i in 1:k){
      group[i,1] <- agrupamento[j,1]
      group[i,2] <- agrupamento[j,2]
      j = j - 1
    }
    mediasd <- sort(medias, decreasing = TRUE)
    resultado <- data.frame(r = 0, f = group[,1], m = mediasd)
    res <- 1
    for (i in 1:(k-1)) {
      if (group[i,2] != group[i + 1,2]) {
        resultado$r[i] <- letters[res]
        res <- res + 1
        if (i == (length(resultado$r) - 1)) {
          resultado$r[i + 1] <- letters[res]
        }
      }
      else {
        resultado$r[i] <- letters[res]
        if (i == (length(resultado$r) - 1)) {
          resultado$r[i + 1] <- letters[res]
        }
      }
    }
    trat.n <- as.vector(levels(trt))
    resultado[,2] <- trat.n[resultado[,2]]
    names(resultado) <- c("Grupos", "Tratamentos", "Medias")
    cat("\nTeste de Calinski & Corsten baseado na distribuicao F\n------------------------------------------------------------------------\n")
    print(resultado)    
    cat("------------------------------------------------------------------------\n")
}
