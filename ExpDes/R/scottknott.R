scottknott <-
function(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE, main = NULL)
{
 sk <- function(medias,s2,dfr,prob)
 # esta funcao faz os calculos de bo e da probabilidade
 {
	bo <- 0
	si2 <- s2
	defr <- dfr
	parou <- 1
	np <- length(medias) - 1 #numero de particoes

	for (i in 1:np) #ve qual e o maior B0
	{
	g1 <- medias[1:i]
	g2 <- medias[(i+1):length(medias)]
	B0 <- sum(g1)^2/length(g1) + sum(g2)^2/length(g2) - (sum(g1) + sum(g2))^2/length(c(g1,g2))
	if (B0 > bo)
		{
		 bo <- B0
		 parou <- i
		}
	}
	
	g1 <- medias[1:parou]
	g2 <- medias[(parou+1):length(medias)]
	teste <- c(g1,g2)

	sigm2 <- (sum(teste^2) - sum(teste)^2/length(teste) + defr*si2)/(length(teste) + defr)

	lamb <- pi*bo/(2*sigm2*(pi-2))
	
	v0 <- length(teste)/(pi-2)
	
	p <- pchisq(lamb,v0,lower.tail = FALSE)
	
		if (p < prob) {
			for (i in 1:length(g1)){
			cat(names(g1[i]),"\n",file="skresult",append=TRUE)
			}
		cat("*","\n",file="skresult",append=TRUE)
		}

	if (length(g1)>1) #para calcular os demais grupos
	{
	sk(g1,s2,dfr,prob)
	}
	if (length(g2)>1)
	{
	sk(g2,s2,dfr,prob)
	}
}

medias <- sort(tapply(y,trt,mean),decreasing=TRUE)
dfr <- DFerror

	rep <- tapply(y,trt,length)
#	erro <- anova$fitted.values-anova$model[[1]]
#	s0 <- (erro^2)/anova$df.residual
#	s1 <- tapply(s0,anova$model[[vari]],sum)
#	s2 <- sum(s1/rep)
  s0 <- MSerror <-SSerror/DFerror
	s2 <- s0/rep[1]

cat("\nScott-Knott test\n------------------------------------------------------------------------\n")
#"Confidence Level: ",conf.level,"\n",
#"Independent variable: ", which,"\n","\n")

prob <- alpha

sk(medias,s2,dfr,prob)

f <- names(medias)
names(medias) <- 1:length(medias)		#monta a tabela de resultado
resultado <- data.frame("r"=0,"f"=f,"m"=medias)
                                              #muda o valor do r da tabela, se o arquivo tiver vazio a funcao para
if (file.exists("skresult") == FALSE) {stop} else{
	xx <- read.table("skresult")
	file.remove("skresult")
	x <- xx[[1]]
	x <- as.vector(x)
	z <- 1

	for (j in 1:length(x)){
		if (x[j] == "*")	{z <- z+1}
		for (i in 1:length(resultado$f)){
		if (resultado$f[i]==x[j]){
		resultado$r[i] <- z;
		}
		
		}
	}

}
##############################################################################################################################
#Isso permite que sejam colocadas quantas letras forem necessarias, e nao apenas o numero de letras existentes no alfabeto (26).
#Ou seja, permite comparar qualquer quantidade de tratamentos (Eric)
letras<-letters
 if(length(resultado$r)>26) {
 l<-floor(length(resultado$r)/26)
  for(i in 1:l) letras<-c(letras,paste(letters,i,sep=''))
                            }
##############################################################################################################################
#coloca as letras
res <- 1
for (i in 1:(length(resultado$r)-1))
	{
		if (resultado$r[i] != resultado$r[i+1]){
			resultado$r[i] <- letters[res]
			res <- res+1
				if (i == (length(resultado$r)-1)){
				resultado$r[i+1] <- letters[res]			
				}	
		}
		else{
			resultado$r[i] <- letters[res]
				if (i == (length(resultado$r)-1)){
				resultado$r[i+1] <- letters[res]			
				}	
		}
	}
names(resultado) <- c("Groups","Treatments","Means")
print(resultado)
cat('------------------------------------------------------------------------\n')
}
