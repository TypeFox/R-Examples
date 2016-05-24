Biplot <- function(Data, alfa=0.5, Title=NA, Label_x=NA, Label_y=NA, Color="s") {
  # Rotina para gerar Biplot desenvolvida 
  # por Paulo Cesar Ossani em 20/06/2015
  
  # Entrada:
  # Data  - Dados para plotagem.
  # alfa  - Representatividade dos individuos (alfa), 
  #         representatividade das variaveis (1-alfa). 
  #         Sendo 0.5 o default.
  # Title  - Titulo para o grafico. Se nao for definido assume texto padrao.
  # Label_x - Rotulo do eixo X. Se nao for definido assume texto padrao.
  # Label_y - Rotulo do eixo Y. Se nao for definido assume texto padrao.
  # Color  - "s" para graficos coloridos (default)
  #          "n" para graficos em preto e branco.
  
  # Retorna:
  # Grafico Biplot.
  # Md - Matriz autovalores.
  # Mu - Matriz U (autovetores).
  # Mv - Matriz V (autovetores).
  # Coor_I - Coordenadas dos individuos.
  # Coor_V - Coordenadas das variaveis.
  # PVar   - Proporcao dos componentes principais.
  
  ##### INICIO - Informacoes usadas nos Graficos #####
  
  if (!is.data.frame(Data)) 
     stop("Entrada para 'Data' esta incorreta, deve ser do tipo dataframe. Verifique!")
  
  if (!is.numeric(alfa) || alfa < 0 || alfa > 1)
     stop("Entrada para 'alfa' esta incorreta, deve ser numerica, com valor entre 0 e 1. Verifique!")
  
  if (!is.character(Title) && !is.na(Title))
     stop("Entrada para 'Title' esta incorreta, deve ser do tipo caracter ou string. Verifique!")
  
  if (!is.character(Label_x) && !is.na(Label_x))
     stop("Entrada para 'Label_x' esta incorreta, deve ser do tipo caracter ou string. Verifique!")
  
  if (!is.character(Label_y) && !is.na(Label_y))
     stop("Entrada para 'Label_y' esta incorreta, deve ser do tipo caracter ou string. Verifique!")

  Color  = ifelse(Color=="s","S",ifelse(Color=="n","N", Color))    # transforma em maiusculo
  
  if (Color!="S" && Color!="N") 
     stop("Entrada para 'Color' esta incorreta, deve ser do tipo caracter, sendo 's' ou 'n'. Verifique!")
  
  if (is.na(Title)) Title = "Grafico Biplot" 
  
  MData = as.matrix(Data) # transforma dados em matriz
  
  ### Centraliza os dados na media
  Media <- apply(MData, 2, mean) # medias por colunas
  MData <- sweep(MData, 2, Media, FUN = "-") # Centraliza na media
    
  ### Decompondo Singularmente a Matriz de Dados
  dim  <- 2 # dimenssao 
  Mdvs <- svd(MData) # Matriz de Decomposicao Valor Singular
  Md = Mdvs$d # Matriz autovalores
  Mu = Mdvs$u # Matriz U (autovetores)
  Mv = Mdvs$v # Matriz V (autovetores)
  
  Coor_I <- Mu[,1:dim]%*%diag(Md[1:dim])^alfa     # coordenadas individuos
  Coor_V <- Mv[,1:dim]%*%diag(Md[1:dim])^(1-alfa) # coordenadas variaveis
  
  PVar <- (Md^2/sum(Md^2)) * 100 # Proporcao dos primeiros (dim) componentes principais
  
  if (is.na(Label_x))
     Label_x = paste("Primeiro Componente (",round(PVar[1],2),"%)",sep="")

  if (is.na(Label_y))
     Label_y = paste("Segundo Componente (",round(PVar[2],2),"%)",sep="")
  
  MaxX <- max(Coor_I[,1],Coor_V[,1])+1 # Dimenssoes maximas das linhas
  MinX <- min(Coor_I[,1],Coor_V[,1])-1 # Dimenssoes minimas das linhas
  MaxY <- max(Coor_I[,2],Coor_V[,2])+1 # Dimenssoes maximas das colunas
  MinY <- min(Coor_I[,2],Coor_V[,2])-1 # Dimenssoes minimas das colunas
  
  cor  <- c(1) # cor inicial
  
  ##### INICIO - Grafico Biplot #####  
  plot(0,0, # Plota as variaveis
       xlab = Label_x,  # Nomeia Eixo X
       ylab = Label_y,  # Nomeia Eixo Y
       main = Title,    # Titulo
       asp = 1,           # Aspecto do grafico
       cex=0,             # Tamanho dos pontos
       xlim=c(MinX,MaxX), # Dimensao para as linhas do grafico
       ylim=c(MinY,MaxY)) # Dimensao para as colunas do grafico
  
  abline(h = 0, v = 0, cex = 1.5, lty = 2) # cria o eixo central
  
  NomeVar <- colnames(MData) # nomes das variaveis
  for (i in 1:nrow(Coor_V)) {  # foi necessario criar este for para poder colocar cores diferentes para cada variavel
    arrows(0,0,Coor_V[i,1],Coor_V[i,2], lwd = 2, code = 2, angle = 10, col = ifelse(Color=="S", cor + i, 1)) # cria a seta apontando para cada variavel  
    #text(Coor_V[i,1], Coor_V[i,2], cex = 1, pos = 4, NomeVar[i], col = ifelse(Color=="S", cor + i, 1), xpd = TRUE)  # Coloca os nomes das variaveis
  }
  
  if (Color=="S") cor <- c((cor+1):(length(NomeVar)+1))
  LocLab(Coor_V[,1:2], NomeVar, col = cor)  # Coloca os nomes das variaveis
  
  NomeVar <- rownames(MData) # nomes os individuos
  LocLab(Coor_I[,1:2], NomeVar) # Coloca os nomes dos individuos
  #for (i in 1:nrow(Coor_I)) 
    #text(Coor_I[i,1], Coor_I[i,2], cex = 1, pos = 3, NomeVar[i], xpd = TRUE)  # Coloca os nomes dos individuos
  
  points(Coor_I,    # Coloca pontos nas posicoes dos individuos
         asp = 1,   # Aspecto do grafico
         pch = 15,  # Formato dos pontos 
         cex = 1.2, # Tamanho dos pontos         
         col = ifelse(Color=="S",2,1))
 
  ##### FIM - Grafico Biplot #####
  
  Lista <- list(Md=Md, Mu=Mu, Mv=Mv, Coor_I=Coor_I,
                Coor_V=Coor_V, PVar=PVar)
  
  return (Lista) 
  
}

