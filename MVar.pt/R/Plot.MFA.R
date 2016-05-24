Plot.MFA <- function(MFA,Titles = matrix(NA,1,3), PosLeg=2, BoxLeg="s", Color="s",NamArr="n") {
  # Rotina para Plotar Graficos do Metodo MFA desenvolvida 
  # por Paulo Cesar Ossani em 09/2013 a 01/2014
  
  # Entrada:
  # MF - Dados da funcao MFA
  # Titles - Titulos para os graficos. Se nao for definido assume texto padrao.
  # PosLeg - 1 para legenda no canto superior esquerdo
  #          2 para legenda no canto superior direito (default)
  #          3 para legenda no canto inferior direito
  #          4 para legenda no canto inferior esquerdo
  # BoxLeg - "s" para colocar moldura na legenda (default)
  #          "n" nao coloca moldura na legenda
  # Color  - "s" para graficos coloridos (default)
  #          "n" para graficos em preto e branco
  # NamArr - "s" para colocar nomes pontos na nuvem ao redor do
  #              centroide no Grafico Correspondente a 
  #              Analise Global dos Individuos e Variaveis
  #          "n" Caso contrario (default)
  
  # Retorna:
  # Varios graficos
  
  ##### INICIO - Informacoes usadas nos Graficos #####
  # Cria Titulos para os graficos caso nao existam
  if (!is.character(Titles[1]) || is.na(Titles[1])) Titles[1] = c("Grafico Correspondente a Analise Global dos Individuos")
  if (!is.character(Titles[2]) || is.na(Titles[2])) Titles[2] = c("Grafico Correspondente a Analise Global dos Individuos e Variaveis")
  if (!is.character(Titles[3]) || is.na(Titles[3])) Titles[3] = c("Grafico das Inercias dos Grupos de Variaveis")
  
  Color  = ifelse(Color=="s","S",ifelse(Color=="n","N",Color))    # transforma em maiusculo
  BoxLeg = ifelse(BoxLeg=="s","S",ifelse(BoxLeg=="n","N",BoxLeg)) # transforma em maiusculo
  NamArr = ifelse(NamArr=="s","S",ifelse(NamArr=="n","N",NamArr)) # transforma em maiusculo
  
  if (PosLeg<1 || PosLeg>4)
     stop("Entrada para posicao da legenda 'PosLeg' esta incorreta. Verifique!")
  
  if (BoxLeg!="S" && BoxLeg!="N") 
     stop("Entrada para moldura da legenda 'BoxLeg' esta incorreta, deve ser do tipo caracter, sendo 's' ou 'n'. Verifique!")
  
  if (Color!="S" && Color!="N") 
     stop("Entrada para 'Color' esta incorreta, deve ser do tipo caracter, sendo 's' ou 'n'. Verifique!")

  if (NamArr!="S" && NamArr!="N") 
     stop("Entrada para 'NamArr' esta incorreta, deve ser do tipo caracter, sendo 's' ou 'n'. Verifique!")
  
  Grupos     = MFA$MatrixG  # tamanho de cada grupo
  NomeGrupos = MFA$MatrixNG # nomes de cada grupo
  NomeLinhas = rownames(MFA$MatrixF) # nomes das linhas que formam os dados
  NumGrupos  = length(NomeGrupos) # Numero de Grupos
  cor        = 1 # cor inicial
  DescEixo1  = paste("Primeiro Componente Principal (",round(MFA$MatrixA[1,2],2),"%)",sep="")
  DescEixo2  = paste("Segundo Componente Principal (",round(MFA$MatrixA[2,2],2),"%)",sep="")
  
  if (PosLeg==1) PosLeg="topleft"     # posicao das legendas nos graficos
  if (PosLeg==2) PosLeg="topright"
  if (PosLeg==3) PosLeg="bottomright"
  if (PosLeg==4) PosLeg="bottomleft"
  
  BoxLeg = ifelse(BoxLeg=="S","o","n") # moldura nas legendas, "n" sem moldura, "o" com moldura
  
  Color_a = ifelse(Color=="S","red","black") # cores nos pontos dos graficos
  Color_b = cor # coreas para letras das legendas e suas representacoes no grafico
  if (Color=="S") Color_b = (cor+1):(cor+NumGrupos)
  #####   FIM - Informacoes usadas nos Graficos  #####
  
  ##### INICIO - Plotagem dos Autovalores #####
  mp <- barplot(MFA$MatrixA[,1],names.arg=paste(round(MFA$MatrixA[,2],2),"%",sep=""),main = "Autovalor")
  ##### FIM - Plotagem dos Autovalores #####
  
  ##### INICIO - Plotagem da Analise Global #####
  plot(MFA$MatrixF, # cria grafico para as coordenadas principais da Analise Global
       xlab = DescEixo1,  # Nomeia Eixo X
       ylab = DescEixo2,  # Nomeia Eixo Y
       main = Titles[1],  # Titulo
       asp = 2,           # Aspecto do Grafico
       pch = 15,          # Formato dos pontos 
       cex=1,             # Tamanho dos pontos
       xlim=c(min(MFA$MatrixF[,1])-0.1,max(MFA$MatrixF[,1])+0.1), # Dimensao para as linhas do grafico
       ylim=c(min(MFA$MatrixF[,2]-0.1),max(MFA$MatrixF[,2])+0.1), # Dimensao para as colunas do grafico
       col = ifelse(Color=="S","red","black"))       # Cor dos pontos
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  LocLab(MFA$MatrixF[,1:2], NomeLinhas)  # Coloca os nomes dos pontos das coordenadas principais das linhas
  #text(MFA$MatrixF, cex=1, NomeLinhas, pos=3, xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais das linhas
  ##### FIM - Plotagem da Analise Global #####
  
  ##### INICIO - Plotagem da Analise por Grupo Juntamente com a Analise Global #####
  ## INICIO - Encontra as dimensoes maximas e minimas para as colunas e linhas ##
  MLC <- MFA$MatrixF[,1:2]
  for (i in 1:length(MFA$MatrixEFG)) 
    MLC <- rbind(MLC,MFA$MatrixEFG[[i]][,1:2])
  maxX = max(MLC[,1]) # Dimenssoes maximas das linhas do grafico
  minX = min(MLC[,1]) # Dimenssoes minimas das linhas do grafico
  maxY = max(MLC[,2]) # Dimenssoes maximas das colunas do grafico
  minY = min(MLC[,2]) # Dimenssoes minimas das colunas do grafico
  ## FIM - Encontra as dimensoes maximas e minimas para as colunas e linhas ##
  
  plot(MFA$MatrixF, # cria grafico para as coordenadas principais da Analise por Grupo
       xlab = DescEixo1,  # Nomeia Eixo X
       ylab = DescEixo2,  # Nomeia Eixo Y
       main = Titles[2], # Titulo
       asp = 1,           # Aspecto do grafico
       pch = 15,          # Formato dos pontos 
       cex=1.2,           # Tamanho dos pontos
       xlim=c(minX,maxX), # Dimensao para as linhas do grafico
       ylim=c(minY,maxY), # Dimensao para as colunas do grafico
       col = Color_a)     # Cor dos pontos
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  LocLab(MFA$MatrixF[,1:2], NomeLinhas)  # Coloca os nomes dos pontos das coordenadas principais da analise global
  #text(MFA$MatrixF, cex=1,NomeLinhas, pos=3, xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais da analise global
  ## Acrescenta no grafico da Analise Global as coordenadas principais da Analise por Grupo
  NumObserv = 4 # numero de centroides a considerar para plotagem das orbitas
  NumLinhas = nrow(MFA$MatrixEFG[[1]]) # numero de linhas
  if (NumObserv<NumLinhas) {
    Position = floor(NumLinhas/NumObserv)
    Observ = as.vector(c(rep(1,NumObserv))) # observacoes a serem plotadas orbitando os centroides
    for (i in 1:(length(Observ)-2)) {
      Observ[i+1] = Position*i
    }     
    Observ[length(Observ)]=NumLinhas # observacoes a serem plotadas orbitando os centroides
  }
  
  if (NumObserv>=NumLinhas)
    Observ = 1:NumLinhas  # observacoes a serem plotadas orbitando os centroides
  
  for (i in 1:length(MFA$MatrixEFG)) {
    if (NamArr=="N") 
      points(MFA$MatrixEFG[[i]][Observ,1:2], pch = (2 + ifelse(Color=="S",i,0)), cex = 1.2, col = 1 + ifelse(Color=="S",i,0)) # adiciona ao grafico as coordenadas principais dos Grupos
    else
      LocLab(MFA$MatrixEFG[[i]][Observ,1:2],NomeGrupos[i], col = 1 + ifelse(Color=="S",i,0)) # Coloca os nomes dos pontos das coordenadas principais dos Grupos
      #text(MFA$MatrixEFG[[i]][Observ,1:2], pos=3, cex=1, NomeGrupos[i], col = 1 + ifelse(Color=="S",i,0),xpd = TRUE) # Coloca os nomes dos pontos das coordenadas principais dos Grupos
  }
  
  ## liga os pontos de cada Analise Global com cada ponto da Analise por Grupo
  for (j in 1:length(MFA$MatrixEFG)) 
    segments(MFA$MatrixF[Observ,1], MFA$MatrixF[Observ,2], MFA$MatrixEFG[[j]][Observ,1], MFA$MatrixEFG[[j]][Observ,2], lty = cor + j, col = ifelse(Color=="S",cor + j,cor))
  
  if (NamArr=="N")
    legend(PosLeg, NomeGrupos, lty = (cor+1):(cor+NumGrupos), col = Color_b, text.col = Color_b,
           bty=BoxLeg, text.font = 6, y.intersp = 0.8,xpd = TRUE) # cria a legenda
  ##### FIM - Plotagem de Analise por Grupo Juntamento com a Analise Global #####
  
  ##### INICIO - Plotagem das Correlacoes dos Componentes Principais com as Variaveis Originais #####
  plot(0,0, # cria grafico para as coordenadas das Correlacoes dos Componentes Principais com as Variaveis Originais
       xlab = DescEixo1, # Nomeia Eixo X
       ylab = DescEixo2, # Nomeia Eixo Y
       main = "Circulo de Correlacao", # Titulo
       asp = 1,           # Aspecto do grafico
       cex=0,             # Tamanho dos pontos
       xlim=c(-1.1,1.1),  # Dimensao para as linhas do grafico
       ylim=c(-1.1,1.1))  # Dimensao para as colunas do grafico
  
  symbols(0, 0, circles = 1, inches = FALSE, fg = 1, add = TRUE) # cria um circulo
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  j  <- 1         # coluna inicial do Grupo de variaveis
  k  <- Grupos[1] # coluna final do Grupo de variaveis
  for (i in 1:NumGrupos) {  # foi necessario criar este for para poder colocar cores diferentes para cada Grupo de variaveis
    
    arrows(0,0,MFA$MatrixCCP[1,j:k],MFA$MatrixCCP[2,j:k], lty=i, code = 2, angle = 10, col = ifelse(Color=="S",cor + i,cor)) # cria a seta apontando para cada coordenada principal
    
    if (is.null(colnames(MFA$MatrixCCP[,j:k])))
      NomeVar<- paste("Comp.", 1:Grupos[i], sep = "") # Nomeia as colunas
    else
      NomeVar<- colnames(MFA$MatrixCCP[,j:k])
    
    LocLab(t(MFA$MatrixCCP[,j:k]), NomeVar, col = ifelse(Color=="S",cor + i,cor)) # Coloca os nomes dos pontos das coordenadas principais
    #text(t(MFA$MatrixCCP[,j:k]), cex=1, pos=3, NomeVar, col = ifelse(Color=="S",cor + i,cor), xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais
    
    j <- j + Grupos[i]  # coluna inicial do Grupo de variaveis
    
    k <- k + Grupos[i+ifelse(i!=NumGrupos,1,0)]  # coluna final do Grupo de variaveis  
    
  }
  
  legend(PosLeg, NomeGrupos, lty = cor:(cor+NumGrupos), col = Color_b, text.col = Color_b,
         bty=BoxLeg, text.font = 6, y.intersp = 0.8,xpd = TRUE) # cria a legenda
  ##### FIM - Plotagem das Correlacoes dos Componentes Principais com as Variaveis Originais #####
  
  ##### INICIO - Plotagem das Inercias Parciais/Escores das Variareis #####
  VlrMinX = ifelse(min(MFA$MatrixEscVar[,1])>0,-0.01,min(MFA$MatrixEscVar[,1])) # Valor minimo para a linha X
  VlrMinY = ifelse(min(MFA$MatrixEscVar[,2])>0,-0.01,min(MFA$MatrixEscVar[,2])) # Valor minimo para a linha Y
  VlrMaxX = 1.01 # Valor maximo para a linha X
  VlrMaxY = 1.01 # Valor maximo para a linha Y
  plot(MFA$MatrixEscVar, # cria grafico para as coordenadas Inercias Parciais/Escores das Variareis
       xlab = DescEixo1,  # Nomeia Eixo X
       ylab = DescEixo2,  # Nomeia Eixo Y
       main = Titles[3], # Titulo
       asp = 1,           # Aspecto do grafico
       pch = 15,          # Formato dos pontos 
       cex=1,             # Tamanho dos pontos
       xlim=c(VlrMinX,VlrMaxX), # Dimensao para as linhas do grafico
       ylim=c(VlrMinY,VlrMaxY), # Dimensao para as colunas do grafico
       col = Color_a)       # Cor dos pontos
  
  abline(h = 0, v=0, cex = 1.5, lty=2) # cria o eixo central
  
  LocLab(MFA$MatrixEscVar[,1:2],rownames(MFA$MatrixEscVar))  # Coloca os nomes dos pontos das coordenadas principais das linhas
  #text(MFA$MatrixEscVar,cex=1, rownames(MFA$MatrixEscVar), pos=3, xpd = TRUE)  # Coloca os nomes dos pontos das coordenadas principais das linhas
  ##### FIM - Plotagem das Inercias Parciais/Escores das Variareis #####
}