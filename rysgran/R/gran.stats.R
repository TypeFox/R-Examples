gran.stats <-function(data, output = "phi", method="folk", verbal=FALSE, lang="en-US"){
  
  if (is.null(data)) # Verifica de existe uma matriz
    stop("Please provide a data frame")
  if (any(sapply(1:ncol(data), function(i) !is.numeric(data[,i])))) # Verifica se todos os valores sao numericos
    stop("Please check your data. Data must be numeric only. The sample labels should be used as row names of your data file.")
  
  tab.ver <- data[,order(data[1,], decreasing=FALSE)] # ordena data baseado nas classes (do menor para o maior)
  if (sum(data[1,]) > 550){ #se a soma das classes for maior que 550 entao as classes nao estao em phi 
    tab.ver[1,]<-(-log2(tab.ver[1,]/1000)) #transforma em phi
    tab.ver <- tab.ver[,order(tab.ver[1,], decreasing=FALSE)]  # reordena baseado nas classes (do menor para o maior)
  }
  
  clz <- as.numeric(tab.ver[1,]) # armazena as classes granulometricas
  tab <- tab.ver[2:nrow(data),] # retira as classes da matrix
  tab.res <- as.data.frame(tab[,1],row.names=rownames(tab)) #cria a matrix que ira receber os resultados
  
  if (method!="moment" && method!="folk" && method!="mcA" && method!="mcB" && method!="trask" && method!="otto")
    stop("The methods supported are moment, folk, mcA, mcB, trask and otto. See help(gran.stats) for further details.")
  
if (method=="moment") # se for utilizar o metodo das medidas dos momentos
{ 
 if (verbal==FALSE) # nao vai ter as descricoes das classes
 {
  tab.res <- as.data.frame(cbind(tab.res,NA,NA,NA,NA,NA)) 
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p") # portugues, tudo codificado para linguagem ascii
   colnames(tab.res) <- c("M\u00E9dia","Sele\u00E7\u00E3o","Assimetria","Curtose","Quinto.Momento","Sexto.Momento")
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e") #ingles
   colnames(tab.res) <- c("Mean","Sorting","Skewness","Kurtosis","Fifth.Moment","Sixth.Moment")
 }
  else # vai ter as descricoes das classes
 {
  tab.res <- as.data.frame(cbind(tab.res,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")# portugues, tudo codificado para linguagem ascii
   colnames(tab.res) <- c("M\u00E9dia","Class.M\u00E9dia","Sele\u00E7\u00E3o","Class.Sele\u00E7\u00E3o","Assimetria","Class.Assimetria","Curtose","Class.Curtose","Quinto.Momento","Sexto.Momento")
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")#ingles
   colnames(tab.res) <- c("Mean","Verbal.Mean","Sorting","Verbal.Sorting","Skewness","Verbal.Skewness","Kurtosis","Verbal.Kurtosis","Fifth.Moment","Sixth.Moment")
 }
}
 else # para qualquer outro metodo (folk, trask, etc)
{
 if (verbal==FALSE)# nao vai ter as descricoes das classes
 {
  tab.res <- as.data.frame(cbind(tab.res,NA,NA,NA,NA))
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")# portugues, tudo codificado para linguagem ascii
   colnames(tab.res) <- c("M\u00E9dia","Mediana","Sele\u00E7\u00E3o","Assimetria","Curtose")
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")#ingles
   colnames(tab.res) <- c("Mean","Median","Sorting","Skewness","Kurtosis")
 }
 else # vai ter as descricoes das classes
 {
  tab.res <- as.data.frame(cbind(tab.res,NA,NA,NA,NA,NA,NA,NA,NA))
  if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")# portugues, tudo codificado para linguagem ascii
   colnames(tab.res) <- c("M\u00E9dia","Class.M\u00E9dia","Mediana","Sele\u00E7\u00E3o","Class.Sele\u00E7\u00E3o","Assimetria","Class.Assimetria","Curtose","Class.Curtose")
  if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")#ingles
   colnames(tab.res) <- c("Mean","Verbal.Mean","Median","Sorting","Verbal.Sorting","Skewness","Verbal.Skewness","Kurtosis","Verbal.Kurtosis")
 }
}
for (j in 1:nrow(tab)){# para cada amostra/linha
  vbz <- as.numeric(tab[j,]) #armazena os pesos/porcentagens da amostra
  if (method=="moment"){ # se for utilizar o metodo das medidas dos momentos
    table<-moment(clz, vbz, output=output) # chama a funcao moment e armazena em table
    if (verbal==FALSE){ # se nao quiser descricao pega todos os resultados da funcao moment
      tab.res[j,1] <- table$mean
      tab.res[j,2] <- table$sort
      tab.res[j,3] <- table$skew
      tab.res[j,4] <- table$kurt
      tab.res[j,5] <- table$A5
      tab.res[j,6] <- table$A6
      }
    if (verbal ==TRUE){ # se qusier descricao das estatisticas chama as funcoes de descricoes
      tab.res[j,1] <- table$mean
      tab.res[j,2] <- class.mean(table$mean, lang=lang, output=output)
      tab.res[j,3] <- table$sort
      tab.res[j,4] <- class.sort(table$sort, method=method, lang=lang, output=output)
      tab.res[j,5] <- table$skew
      tab.res[j,6] <- class.skew(table$skew, method=method, lang=lang, output=output)
      tab.res[j,7] <- table$kurt
      tab.res[j,8] <- class.kurt(table$kurt, method=method, lang=lang, output=output)
      tab.res[j,9] <- table$A5
      tab.res[j,10] <- table$A6
      }
    }
  if (method=="folk"|method=="mcA"|method=="mcB"|method=="trask"|method=="otto"){ #para todos os outros metodos, menos de medida dos momentos
    table <- other.methods(clz, vbz, method=method, output=output) # chama a funcao other.methods e armazena em table
    if (verbal==FALSE){# se nao quiser descricao pega todos os resultados da funcao other.methods
      tab.res[j,1] <- table$mean
      tab.res[j,2] <- table$median
      tab.res[j,3] <- table$sort
      tab.res[j,4] <- table$skew
      tab.res[j,5] <- table$kurt
      }
    if (verbal ==TRUE){# se qusier descricao das estatisticas chama as funcoes de descricoes
      tab.res[j,1] <- table$mean
      tab.res[j,2] <- class.mean(table$mean, lang=lang, output=output)
      tab.res[j,3] <- table$median
      tab.res[j,4] <- table$sort
      tab.res[j,5] <- class.sort(table$sort, method=method, lang=lang, output=output)
      tab.res[j,6] <- table$skew
      tab.res[j,7] <- class.skew(table$skew, method=method, lang=lang, output=output)
      tab.res[j,8] <- table$kurt
      tab.res[j,9] <- class.kurt(table$kurt, method=method, lang=lang, output=output)
      }
    }
  }
  return (tab.res) #retorna um data.frame com as estatisticas granulometricas (nas colunas) para cada amostra (linhas)
}
