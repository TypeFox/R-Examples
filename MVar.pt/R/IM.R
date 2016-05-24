IM <- function(Data, Names = "s") {
  # Converte para variaveis Dummy para execucao da Analise
  # de Correspondencia Multipla, ou seja, em 0 e 1, caso dados nominais
  # Esta funcao e usada na funcao que balanceia dados Categoricos
  
  # Entrada:
  # Data  - Dados Categoricos 
  # Names - "s" para incluir os nomes das variaveis nos niveis da Matriz Indicadora (default)
  #         "n" nao inclui
  
  # Retorna:
  # Dados - Dados convertidos em Matriz Indicadora
  
  if (!is.data.frame(Data)) 
     Data = as.data.frame(Data)
  
  Names = ifelse(Names=="s","S",ifelse(Names=="n","N",Names)) # transforma em maiusculo
  
  if (Names!="S" && Names!="N") 
     stop("Entrada para 'Names' esta incorreta, deve ser do tipo caracter, sendo 's' ou 'n'. Verifique!")
  
  NumLinha  <- nrow(Data)  # Numero de linhas na tabela
  
  for (k in 1:ncol(Data)) {
    
    MConver   <- factor(Data[,k]) # Matriz com os dados para a conversao
    
    Nivel     <- levels(MConver)  # Nomes dos Niveis
    
    Qtd_Nivel <- nlevels(MConver) # Quantidade de Niveis
    
    MDummy = matrix(0,NumLinha,Qtd_Nivel) # Cria Matriz Vazia com elementos zero
    
    for (i in 1:Qtd_Nivel)
      
      for ( j in 1:NumLinha)
        
        if (MConver[j]==Nivel[i]) MDummy[j,i] <- 1
    
    if (Names=="S")
      colnames(MDummy) <- paste("(",Nivel,") ",colnames(Data[k]),sep="")# Nomeia as colunas 
    
    if (Names=="N")
      colnames(MDummy) <- Nivel # Nomeia as colunas  
    
    if (k==1) MFinal <- MDummy
    
    else
      
      MFinal <- cbind(MFinal,MDummy)
    
  }
  
  return(Dados=MFinal)
}