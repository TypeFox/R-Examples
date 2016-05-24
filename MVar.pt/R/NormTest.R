NormTest <- function(Data, Sign = 0.05) {
  # Executa teste para verificar a normalidade dos dados baseado
  # no teste de coeficiente de assimetria
  # Desenvolvida por Paulo Cesar Ossani em 22/06/2013
  # Ver Livro Daniel Furtado pg. 115 e Rencher pg. 114
  
  # Entrada
  # Data - Dados a serem analisados
  # Sign - Grau de significancia do teste (default 5%)
  
  # Saida: Resultado do teste
  
  if (!is.data.frame(Data) && !is.matrix(Data)) 
     stop("Entrada 'Data' esta incorreta, deve ser do tipo dataframe ou matriz. Verifique!")
  
  if (!is.numeric(Sign)) 
     stop("Entrada para 'Sign' esta incorreta, deve ser numerica com valores entre 0 e 1. Verifique!")
  
  if (Sign<=0 || Sign>1) 
     stop("Entrada para 'Sign' esta incorreta, deve ser valores entre 0 e 1. Verifique!")
  
  n <- ncol(Data)*nrow(Data) # numero de elementos amostrais
  
  p <- ncol(Data)  # numero de parametros
  
  gl =  p*(p+1)*(p+2)/6 # grau de liberdade
  
  Media = as.vector(apply(Data,2,mean))  # Data medias das colunas
  
  G     = t(t(Data)-Media)%*%solve(cov(Data))%*%(t(Data)-Media)
  
  B1p   = sum((diag(G))^3/n^2)
  
  Chi.Quad.Observado <- n*B1p/6 # Estatistica do Teste
  
  qt = qchisq(1-Sign,gl,ncp=0)
  
  cat(paste("Grau de liberdade observado:", round(gl,7)),"\n")
  
  cat(paste("Valor da estatistica do teste Qui-quadrado (Chiq1):", round(Chi.Quad.Observado,7)),"\n")
  
  cat(paste("Valor Qui-quadrado observado (Chiq2) com", Sign*100,"% de significancia:", round(qt,7)),"\n")
  
  if (Chi.Quad.Observado<=qt) cat("Como Chiq1 <= Chiq2, VERIFICA-SE a normalidade dos dados.\n")
  
  if (Chi.Quad.Observado>qt) cat("Como Chiq1 > Chiq2, NAO VERIFICA-SE a normalidade dos dados.\n")
  
  cat("Valor-p:", pchisq(Chi.Quad.Observado,gl,ncp=0, lower.tail = F))
  
}