Plot.PCA <- function(PC, Titles = matrix(NA,1,2), Color = "s", LinLab = NULL) {
  # Rotina para Plotar Graficos do Metodo PCA desenvolvida 
  # por Paulo Cesar Ossani em 11/2014
  
  # Entrada:
  # PC     - Dados da funcao PCA
  # Titles - Titulos para os graficos
  # Color  - "s" para graficos coloridos (default)
  #          "n" para graficos em preto e branco
  # LinLab - Vetor com o rotulo para as linhas, se nao
  #          informado retorna o padrao dos dados.
  
  # Retorna:
  # Varios graficos
  
  ##### INICIO - Informacoes usadas nos Graficos #####
  # Cria Titulos para os graficos caso nao existam
  if (!is.character(Titles[1]) || is.na(Titles[1])) Titles[1] = c("Grafico Correspondente as Linhas (Observacoes)")
  if (!is.character(Titles[2]) || is.na(Titles[2])) Titles[2] = c("Grafico Correspondente as Colunas (Variaveis)")
  
  Color  = ifelse(Color=="s","S",ifelse(Color=="n","N",Color))    # transforma em maiusculo
  
  if (Color!="S" && Color!="N")
    stop("Entrada para 'Color' esta incorreta, deve ser do tipo caracter, sendo 's' ou 'n'. Verifique!")
  
  if (!is.null(LinLab) && length(LinLab)!=nrow(PC$MatrixEsc))
    stop("O numero elementos do rotulo para linhas 'LinLab' difere do numero de linhas da base de dados. Verifique!")
  
  if (is.null(LinLab))
    LinLab <- rownames(PC$MatrixEsc)
  
  DescEixo1  = paste("Primeira Coordenada Principal (",round(PC$MatrixAutoVlr[1,2],2),"%)",sep="")
  DescEixo2  = paste("Segunda Coordenada Principal (",round(PC$MatrixAutoVlr[2,2],2),"%)",sep="")
  #####   FIM - Informacoes usadas nos Graficos  #####
  
  ##### INICIO - Plotagem dos Autovalores #####
  mp <- barplot(PC$MatrixAutoVlr[,1],names.arg=paste(round(PC$MatrixAutoVlr[,2],2),"%",sep=""),main = "Autovalor")
  ##### FIM - Plotagem dos Autovalores #####
  
  ##### INICIO - Plotagem dos Dados das linhas #####
  plot(PC$MatrixEsc, # cria grafico para as coordenadas principais das linhas
       xlab = DescEixo1,  # Nomeia Eixo X
       ylab = DescEixo2,  # Nomeia Eixo Y
       main = Titles[1],  # Titulo
       asp = 1,           # Aspecto do Grafico
       pch = 15,          # Formato dos pontos 
       cex = 1,           # Tamanho dos pontos
       xlim=c(min(PC$MatrixEsc[,1])-0.05,max(PC$MatrixEsc[,1])+0.05), # Dimensao para as linhas do grafico
       ylim=c(min(PC$MatrixEsc[,2])-0.05,max(PC$MatrixEsc[,2])+0.05), # Dimensao para as colunas do grafico
       col = ifelse(Color=="S","red","black"))  # Cor dos pontos
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  #text(PC$MatrixEsc, cex = 1, pos = 3, LinLab)  # Coloca os nomes dos pontos das coordenadas principais das linhas
  LocLab(PC$MatrixEsc, cex = 1, LinLab)
  ##### FIM - Plotagem dos Dados das linhas #####
  
  ##### INICIO - Plotagem das Correlacoes dos Componentes Principais com as Variaveis Originais #####
  plot(0,0, # cria grafico para as coordenadas das Correlacoes dos Componentes Principais com as Variaveis Originais
       xlab = DescEixo1, # Nomeia Eixo X
       ylab = DescEixo2, # Nomeia Eixo Y
       main = Titles[2], # Titulo
       asp = 1,          # Aspecto do Grafico
       cex=0,            # Tamanho dos pontos
       xlim=c(-1.1,1.1), # Dimensao para as linhas do grafico
       ylim=c(-1.1,1.1)) # Dimensao para as colunas do grafico
  
  symbols(0, 0, circles = 1, inches = FALSE, fg = 1, add = TRUE) # cria um circulo
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  arrows(0,0,PC$MatrixCCP[1,],PC$MatrixCCP[2,], lty=1, code = 2, angle = 10, col = ifelse(Color=="S","Blue","Black")) # cria a seta apontando para cada coordenada principal
  
  #text(t(PC$MatrixCCP), cex=1, colnames(PC$MatrixCCP) , col = ifelse(Color=="S","Blue","Black"), pos = 3, xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais
  LocLab(t(PC$MatrixCCP), cex=1, colnames(PC$MatrixCCP) , col = ifelse(Color=="S","Blue","Black"), xpd = TRUE)
  ##### FIM - Plotagem das Correlacoes dos Componentes Principais com as Variaveis Originais #####
}