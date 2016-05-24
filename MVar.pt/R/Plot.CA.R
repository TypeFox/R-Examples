Plot.CA <- function(AC, Titles = matrix(NA,1,3), Color = "s", LinLab = NULL) {
  # Rotina para Plotar Graficos do Metodo AC desenvolvida 
  # por Paulo Cesar Ossani em 11/2014
  
  # Entrada:
  # AC     - Dados da funcao CA
  # Titles - Titulos para os graficos
  # Color  - "s" para graficos coloridos (default)
  #          "n" para graficos em preto e branco
  # LinLab - Vetor com o rotulo para as linhas para dados de frequencia,
  #          se nao informado retorna o padrao dos dados.
  
  # Retorna:
  # Varios graficos
  
  ##### INICIO - Informacoes usadas nos Graficos #####
  # Cria Titulos para os graficos caso nao existam
  if (!is.character(Titles[1]) || is.na(Titles[1])) Titles[1] = c("Grafico Correspondente as Linhas (Observacoes)")
  if (!is.character(Titles[2]) || is.na(Titles[2])) Titles[2] = c("Grafico Correspondente as Colunas(Variaveis)")
  if (!is.character(Titles[3]) || is.na(Titles[3])) Titles[3] = c("Grafico Correspondente as Observacoes e Variaveis")
  
  Color  = ifelse(Color=="s","S",ifelse(Color=="n","N",Color))  # transforma em maiusculo
  
  if (Color!="S" && Color!="N")
    stop("Entrada para 'Color' esta incorreta, deve ser do tipo caracter, sendo 's' ou 'n'. Verifique!")
  
  if (!is.null(LinLab) && length(LinLab)!=nrow(AC$MatrixX) && AC$TypData=="F")
    stop("O numero elementos do rotulo para linhas 'LinLab' difere do numero de linhas da base de dados. Verifique!")
  
  if (is.null(LinLab) && AC$TypData=="F")
    LinLab <- rownames(AC$MatrixX)
  
  DescEixo1  = paste("Primeira Coordenada Principal (",round(AC$MatrixAutoVlr[1,2],2),"%)",sep="")
  DescEixo2  = paste("Segunda Coordenada Principal (",round(AC$MatrixAutoVlr[2,2],2),"%)",sep="")
  #####   FIM - Informacoes usadas nos Graficos  #####
  
  ##### INICIO - Plotagem dos Autovalores #####
  mp <- barplot(AC$MatrixAutoVlr[,1],names.arg=paste(round(AC$MatrixAutoVlr[,2],2),"%",sep=""),main = "Autovalor")
  ##### FIM - Plotagem dos Autovalores #####
  
  ##### INICIO - Plotagem dos Dados das linhas #####
  if (AC$TypData=="F") { # plota se nao for analise de correspondencia multipla
    plot(AC$MatrixX, # cria grafico para as coordenadas principais das linhas
         xlab = DescEixo1,  # Nomeia Eixo X
         ylab = DescEixo2,  # Nomeia Eixo Y
         main = Titles[1],  # Titulo
         asp = 1,           # Aspecto do Grafico
         pch = 15,          # Formato dos pontos 
         cex=1,             # Tamanho dos pontos
         xlim=c(min(AC$MatrixX[,1])-0.1,max(AC$MatrixX[,1])+0.1), # Dimensao para as linhas do grafico
         ylim=c(min(AC$MatrixX[,2]-0.1),max(AC$MatrixX[,2])+0.1), # Dimensao para as colunas do grafico
         col = ifelse(Color=="S","red","black"))             # Cor dos pontos
    
    abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
    
    #text(AC$MatrixX,cex=1, pos=3, LinLab)  # Coloca os nomes dos pontos das coordenadas principais das linhas
    LocLab(AC$MatrixX,cex=1, LinLab)
  }
  ##### FIM - Plotagem dos Dados das linhas #####
  
  ##### INICIO - Plotagem dos Dados das colunas #####
  plot(AC$MatrixY, # cria grafico para as coordenadas principais das linhas
       xlab = DescEixo1,  # Nomeia Eixo X
       ylab = DescEixo2,  # Nomeia Eixo Y
       main = Titles[2],  # Titulo
       asp = 1,           # Aspecto do Grafico
       pch = 16,          # Formato dos pontos 
       cex=1.2,           # Tamanho dos pontos
       xlim=c(min(AC$MatrixY[,1])-0.1,max(AC$MatrixY[,1])+0.1), # Dimensao para as linhas do grafico
       ylim=c(min(AC$MatrixY[,2]-0.1),max(AC$MatrixY[,2])+0.1), # Dimensao para as colunas do grafico
       col = ifelse(Color=="S","blue","black"))             # Cor dos pontos
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  #text(AC$MatrixY, cex=1, pos=3, rownames(AC$MatrixY))  # Coloca os nomes dos pontos das coordenadas principais das colunas
  LocLab(AC$MatrixY, cex=1, rownames(AC$MatrixY))
  ##### FIM - Plotagem dos Dados das colunas #####
  
  ##### INICIO - Plotagem dos Dados das linhas e colunas conjuntamente #####
  if (AC$TypData=="F") {     # plota se nao for analise de correspondencia multipla
    plot(AC$MatrixX,        # cria grafico para as coordenadas principais das linhas
         xlab = DescEixo1,  # Nomeia Eixo X
         ylab = DescEixo2,  # Nomeia Eixo Y
         main = Titles[3],  # Titulo
         asp = 1,           # Aspecto do Grafico
         pch = 15,          # Formato dos pontos 
         cex=1,             # Tamanho dos pontos
         xlim=c(min(AC$MatrixX[,1])-0.1,max(AC$MatrixX[,1])+0.1), # Dimensao para as linhas do grafico
         ylim=c(min(AC$MatrixX[,2]-0.1),max(AC$MatrixX[,2])+0.1), # Dimensao para as colunas do grafico
         col = ifelse(Color=="S","red","black"))             # Cor dos pontos
    
    points(AC$MatrixY, pch = 16, cex = 1.2, col = ifelse(Color=="S","blue","black")) # adiciona ao grafico as coordenadas principais das colunas
    
    abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
    
    #text(AC$MatrixX, cex=1,  pos=3, LinLab)  # Coloca os nomes dos pontos das coordenadas principais das linhas
    LocLab(AC$MatrixX, cex=1,LinLab)
    
    #text(AC$MatrixY, cex=1, pos=3, rownames(AC$MatrixY))  # Coloca os nomes dos pontos das coordenadas principais das colunas
    LocLab(AC$MatrixY, cex=1, rownames(AC$MatrixY))
  }
  ##### FIM - Plotagem dos Dados das linhas e colunas conjuntamente #####
}