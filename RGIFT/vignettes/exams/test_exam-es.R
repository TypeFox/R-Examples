#
#Genera las preguntas para el examen
#
#Tenemos que ejecutar echo("test_exam-es.R")
#para que no nos meta basurilla

library(RGIFT)

#Redefine function to set default values. This has been done to
#simplify the modifications in the original file
pregunta<-function(qtxt, atxt){GIFTMC(qtxt, atxt, wwrong="-33.333")}


#Set locale to UTF-8
#Sys.setlocale(locale="es_ES.UTF-8")

#Fichero en el que escribimos el test
#
#

#Fichero en el que guardamos las preguntas del test en formato GIFT
sink("test_examen.txt")

#
#Práctica 1
#
cat("\n\n$CATEGORY: Ex_Pr1\n\n")

for(i in 1:5)#5 vectores de datos distintos
{

	set.seed(1*i)
	x<-round(rnorm(10,10,1), 2)


     preg<-paste("Sea el vector de datos\\n\\n", deparse(x, width.cutoff=100L),
   "\\n\\n", sep="")

	#Media
	preg1<-paste(preg, "Su media es:")
	sol1<-as.character(mean(x)*c(1,.9, 1.1))
	pregunta(preg1, sol1)

	#Desviación típica
	preg2<-paste(preg, "Su desviación típica es:")
	sol2<-as.character(sd(x)*c(1,.5, 1.5))
	pregunta(preg2, sol2)

	#Varianza
	preg2b<-paste(preg, "Su varianza es:")
	sol2b<-as.character(var(x)*c(1, .5, 1.5))
	pregunta(preg2b, sol2b)

	#Cuantil 0'40
	preg3<-paste(preg, "Su cuantil del 40% es:")
	sol3<-as.character(quantile(x, .4)*c(1,.75,1.25))
	pregunta(preg3, sol3)


	#Cuantil 0'90
	preg4<-paste(preg, "Su cuantil del 90% es:")
	sol4<-as.character(quantile(x, .9)*c(1,.75,1.25))
	pregunta(preg4, sol4)

}


#Algunas preguntas de instrucciones
pregunta("Si quiero saber la longitud de unos datos en R, utilizo la función:",
   c("length()", "long()", "Ninguna de las dos")
)

pregunta("Para realizar un diagrama de barras utilizo la función",
   c("barplot", "plot", "Ninguna de las dos")
)

pregunta("Para realizar un histograma utilizo la función",
   c("Ninguna de las dos","histo", "h")
)

pregunta("Para calcular la raíz cuadrada de 13 utilizo",
   c("sqrt(13)", "sd(13)", "str(13)")
)

pregunta("Para calcular el factorial de 5 utilizo",
   c("factorial(5)", "facto(5)", "f(5)")
)

pregunta("Para generar 10 valores entre 0 y 5, usaré",
   c("seq(0,5,length\\=10)", "0:5", "sequ(0,5)" )
)

pregunta("Para realizar un diagrama de sectores, utilizo la función:",
   c("Ninguna de las dos", "sector()", "tarta()")
)

pregunta("Un diagrama de sectores:",
  c("Es para una variable basada en categorías", 
     "Es para una variable continua", "Ninguna de las dos")
)


#
#Práctica 2
#
cat("\n\n$CATEGORY: Ex_Pr2\n\n")

for(i in 1:5)#5 vectores de datos distintos
{

        set.seed(10*i)
        x<-round(100*rnorm(8,10,1))/100


        preg<-paste("Sea el vector de datos\\n", deparse(x, width.cutoff=100L),
   "\\n", sep="")

        #
        preg1<-paste(preg, "¿Cuánto vale la suma de los elementos menores que la media?")
        sol1<-as.character( sum(x[x<mean(x)])*c(1,.75,.5))
        pregunta(preg1, sol1)

	#
        preg2<-paste(preg, "¿Cuánto vale la suma de los elementos mayores que la media?")
        sol2<-as.character( sum(x[x>mean(x)])*c(1,1.75,1.5))
        pregunta(preg2, sol2)

        #
        preg3<-paste(preg, "¿Cuánto vale la suma de los elementos menores que la mediana?")
        sol3<-as.character( sum(x[x<median(x)])*c(1,.75,.5))
        pregunta(preg3, sol3)

	#
        preg4<-paste(preg, "¿Cuánto vale la suma de los elementos mayores que la mediana?")
        sol4<-as.character( sum(x[x>median(x)])*c(1,.75,1.5))
        pregunta(preg4, sol4)

}


pregunta("Para calcular el determinante de una matriz X usamos:",
  c("det(X)", "deter(X)", "determinate(X)")
)

pregunta("Para calcular la inversa de una matriz X usaremos",
   c("solve(X)", "invers(X)", "1/X")
)

pregunta("Para multiplicar de forma matricial  dos matrices X e Y usaremos:",
   c("X%*%Y", "X*Y", "XY")
)

pregunta("Para multiplicar dos matrices X e Y elemento a elemento usaremos:",
   c("X*Y", "X%*%Y", "X.Y")
)

pregunta("Para calcular una muestra de tamaño 10 de un  vector X sin reposición,
usaremos:",
   c("sample(X,10)", "samples(X,10,TRUE)", "Ninguna de las dos")
)

pregunta("Para calcular una muestra de tamaño 3 de un  vector X con reposición,
usaremos:",
   c("sample(X,3,rep\\=TRUE)", "sample(X,3,FALSE)", "Ninguna de las dos")
)


pregunta("¿Qué instrucción usarías para simular 100 veces el número de figuras obtenidas al extraer 10 cartas de una baraja de 40 cartas con reemplazamiento?",
   c("rbinom(100,10,12/40)", "rbinom(100,10,10/40)", "Ninguna de las dos")
)


#Algunas preguntas sobre distribuciones de probabilidad
nombres<-c("Binomial", "Normal", "Exponencial", "Uniforme") 
distr<-c("binom", "norm", "exp", "unif")
params<-list(c(10, .27), c(4, 2), c(5), c(0, 15))

ndistr<-length(nombres)


for(i in 1:ndistr)
{
	preg<-paste("Sea una variable ", nombres[i], " de parámetros ")
	preg<-paste(preg, paste(params[[i]], collapse=" y "), ". ", sep="")


	#Probabilidad acumulada
	preg1<-paste(preg, "Para calcular la probabilidad de que el valor de la variable aleatoria sea menor o igual que 3, usaremos:", sep="")
	sol1<-paste("p",distr[i], "(3,", paste(params[[i]], collapse= ", "), ")", sep="") 
	sol1<-c(sol1, paste("1-", sol1,sep=""), "Ninguna de las dos")
	
	pregunta(preg1, sol1)

	#1-Probabilidad acumulada
	preg2<-paste(preg, "Para calcular la probabilidad de que el valor de la variable aleatoria sea mayor que 3, usaremos:", sep="")
	sol2<-paste("p",distr[i], "(3,", paste(params[[i]], collapse= ", "), ")", sep="") 
	sol2<-c(paste("1-", sol2,sep=""), sol2, "Ninguna de las dos")
	
	pregunta(preg2, sol2)

	#Cálculo de cuantiles
	qnt<-sample(1:100,1)
	preg3<-paste(preg, "Para calcular el percentil del ",qnt, "% usaremos:", sep="")
	sol3a<-paste("q",distr[i], "(", qnt/100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	sol3b<-paste("r",distr[i], "(", qnt/100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	sol3c<-paste("p",distr[i], "(", qnt/100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	
	pregunta(preg3, c(sol3a, sol3b, sol3c))

	#Simulación
	preg4<-paste(preg, "Para simular 100 observaciones de esa variable aleatoria usaremos:", sep="")
	sol4a<-paste("r",distr[i], "(", 100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	sol4b<-paste("q",distr[i], "(", 100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	sol4c<-paste("p",distr[i], "(", 100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	
	pregunta(preg4, c(sol4a, sol4b, sol4c))
}


#NOTA: POnemos algo más de pintar curvas sobre histogramas??


#
#Práctica 3
#

cat("\n\n$CATEGORY: Ex_Pr3\n\n")


pregunta("Si tengo un factor X y quiero saber el número de datos de cada nivel del factor",
   c("table(X)", "levels(X)", "Ninguna de las anteriores")
)

pregunta("Para crear un factor compuesto por las letras A, B y C donde cada letra se repite 10, 20 y 50 veces respectivamente:",
  c("factor(rep(c(\"A\",\"B\",\"C\"),c(10,20,30)))", 
    "factor(rep(c(\"A\",\"B\",\"C\"),10,20,50))", 
    "Ninguna de las otras")
)

pregunta("Para convertir un vector X en un factor:",
   c("as.factor(X)", "fact(X)", "Ninguna es correcta")
)

pregunta("Para obtener un resumen numérico de un data frame X:",
   c("summary(X)", "Summary(X)", "Ninguna de las otras dos")
)

pregunta("Para conocer los nombres de las variables del data frame X:",
   c("names(X)", "variables(X)", "Ninguna de las dos")
)

#Preguntas del fichero CochesBig.sav
library(foreign)
cochesbig<-read.spss("CochesBig.sav", to.data.frame = TRUE)

variables<-names(cochesbig)
marcas<-levels(cochesbig$marca)

for(v in variables[1:3])
for(m in marcas[1:5])
{
	#Media
	preg1<-paste("Del fichero de CochesBig.sav. La media de ", v, " de la marca ", sep="")
	preg1<-paste(preg1, m, " es", sep="")
	sol1<-as.character(round(mean(cochesbig[cochesbig$marca==m,v])*c(1, .6, .5), 4))
	pregunta(preg1, sol1)


	#Mediana
	preg2<-paste("Del fichero de CochesBig.sav. La mediana de ", v, " de la marca ", sep="")
	preg2<-paste(preg2, m, " es", sep="")
	sol2<-as.character(round(quantile(cochesbig[cochesbig$marca==m, v],.5)*c(1, .6, .5), 4))
	pregunta(preg2, sol2)

	#Cuantil del 0.6
	preg3<-paste("Del fichero de CochesBig.sav. Valor de ", v, " que deja por debajo el 60% de los coches de la marca ", sep="")
	preg3<-paste(preg3, m, " es", sep="")
	sol3<-as.character(round(quantile(cochesbig[cochesbig$marca==m, v],.60)*c(1, .6, .5), 4))
	pregunta(preg3, sol3)

}

#NOTA: Podemos poner alguna pregunta más sobre como comparar modelos, etc.



#
#Práctica 4
#


cat("\n\n$CATEGORY: Ex_Pr4\n\n")

pregunta("La media de 10000 variables aleatorias independientes e igualmente distribuidas, se distribuye según:",
   c("Una normal por el teorema central del límite",
   "Depende de cómo sea la distribución de las variables aleatorias",
   "Ninguna de las otras dos respuestas es correcta")
)

pregunta("La media de 10000 variables aleatorias independientes e igualmente distribuidas, se distribuye según:",
   c("Ninguna de las dos",
   "Para la mayoría de casos según una normal",
   "Depende de cómo sea la distribución de las variables aleatorias")
)


for(m in seq(10, 25, by=5))
{
for(ss in seq(30, 50, by=10))
{
	preg1<-paste("Si tengo una población de media ", m, " y calculo 1000 intervalos de confianza al 95% con muestras de tamaño ", ss, " de esa población:", sep="")
	sol1<-c("Ninguna de las dos", paste("No se puede saber aproximadamente el número de intervalos que contendrán al ", m, "", sep=""), paste("Más del 95% de los intervalos contendrán el valor ", m, sep=""))
	pregunta(preg1, sol1)


	for(dt in seq(20, 50, by=10))
	{
preg2<-paste("Si tengo una variable aleatoria con media ", m, " y desviación típica ", dt, ". ¿Cuál es la desviación típica de la variable media de las muestras de tamaño ", ss, "?", sep="")
sol2<-c(paste(dt, "/sqrt(", ss,")", sep=""),
   paste(dt, ", la misma que la de la variable original", sep=""),
   "Ninguna de las dos")
pregunta(preg2, sol2)

preg3<-paste("Si tengo una variable aleatoria con media ", m, " y desviación típica ", dt, ". ¿Cuál es la media de la variable media de las muestras de tamaño ", ss, "?", sep="")
sol3<-c(paste(m, ", la misma que la de la variable original", sep=""),
paste(dt, "/sqrt(", ss,")", sep=""),
   "Ninguna de las dos")
pregunta(preg3, sol3)


	}

}
}



#Intervalos de confianza  y contrastes con los datos de cochesbig
lm0<-as.list(apply(as.matrix(cochesbig[,1:3]), 2, mean))


for( conf in seq(80, 99, by=5))
{
for(v in variables[1:3])
{
	m0<-round(lm0[[v]]*runif(1, .75,1.25))

	preg1<-paste("Utilizando el archivo CochesBig.sav. La media de ", v,
  " puede ser ", m0, " con una confianza del ", conf,"%", sep="")

	sol1<-c(paste("Sí, ya que ese valor está dentro del intervalo de confianza al ",conf,"%", sep=""),
   paste("No, porque la media vale ",round(mean(cochesbig[,v]), 2), sep=""),
   "Ninguna de las dos")

	tt<-t.test(cochesbig[,v], mu=m0, conf.level=conf/100)

	if(tt$p.value>(1-conf/100))
	{
		pregunta(preg1, sol1)
	}else
	{
		pregunta(preg1, sol1[c(3,1,2)])

	}

	preg2<-paste("Utilizando el archivo CochesBig.sav. Un intervalo de confianza al ", conf,   "% de confianza para la variable ", v, "es:", sep="")
 
	sol2<-c(paste("(", paste(round(as.numeric(tt$conf.int), 4), collapse=", "), ")", sep=""),
   paste("(", paste(5+round(.9*as.numeric(tt$conf.int), 4), collapse=", "), ")", sep=""),
   "Ninguna de las dos preguntas es correcta")

}
}

#AÑADIR: contraste de hipótesis proporciones,



#Contrastes para proporciones

for( conf in seq(80, 99, by=5))
{
for(n in seq(180, 200, by=5))
{
for(x in seq(5, 10, by=5))
{
	p0<-round((x/n)*runif(1,.75,5), 2)

	preg1<-paste("En una fábrica se inspeccionaron ", n, 
   " piezas, de las que ", x, " resultaron ser defectuosas. ",
   "Con un ", conf, "% de confianza, ¿hay evidencia de que la",
   " proporción de piezas defectuosas sea igual a ", p0,"?", sep="")

	ttp<-prop.test(x,n, p=p0, conf.level=conf/100)

	sol1<-c(paste("Sí, porque ese valor está en el intervalo de confianza",sep=""),
   paste("No porque la proporción observada es ", round(x/n, 4), sep=""),
   "Ninguna de las dos")

	if( ttp$p.value>(1-conf/100) )
	{
		pregunta(preg1, sol1)
	}else
	{
		pregunta(preg1, sol1[c(3,2,1)])

	}


	preg2<-paste("En una fábrica se inspeccionaron ", n,
   " piezas, de las que ", x, " resultaron ser defectuosas. ",
"Un intervalo de confianza al ", conf,   "% de confianza para la proporción de piezas defectuosas es", sep="")

        sol2<-c(paste("(", paste(round(as.numeric(ttp$conf.int), 4), collapse=", "), ")", sep=""),
   paste("(", paste(.15+round(1.1*as.numeric(ttp$conf.int), 4), collapse=", "), ")", sep=""),
   "Ninguna de las dos preguntas es correcta")
	pregunta(preg2, sol2)


}
}
}


#NOTA: Podemos poner  que calculen alguna probabilidad sobre 
#la variable promedio de algunos datos, o contrastes de dos poblaciones



#
#Práctica 5
#


cat("\n\n$CATEGORY: Ex_Pr5\n\n")


#Contraste de la media de dos poblaciones
for(conf in seq(80, 99, by=5))
for(m1 in seq(100, 110, by=5))
for(m2 in seq(90, 100, by=5))
for(sigma in seq(5, 5, by=1))
{
	set.seed(conf*m1*m2*sigma)
	hombres<-round(rnorm(10, m1, sigma), 2)
	set.seed(conf*m1*m2*sigma+1)
	mujeres<-round(rnorm(10, m2, sigma), 2)

	preg<-paste("Con los siguientes datos:\\n\\n", 
    "Coeficientes de inteligencia chicas: ",
    deparse(mujeres, width.cutoff=100L), "\\n\\n",
    "Coeficientes de inteligencia chicos: ",
    deparse(hombres, width.cutoff=100L), "\\n\\n", sep="")


   #Contraste
    preg1<-paste(preg, 
   "¿Podemos considerar las medias de chicos y chicas iguales? ",
    "Considera un nivel de significación de ", 1-conf/100, ".", sep="")


   sol1<-c("Sí", "No", "Faltan datos")

   tt<-t.test(mujeres, hombres, conf.level=conf/100)

   if(tt$p.value>(1-conf/100))
   {
	pregunta(preg1, sol1)
   }else
   {
	pregunta(preg1, sol1[c(2,1,3)])
   }



   #p-valor
   preg2<-paste(preg, "El p-valor del contraste de igualdad de medias es:", sep="")
   sol2<-c(as.character(tt$p.value), as.character(1-tt$p.value), "Ninguno de los dos")
   pregunta(preg2, sol2)

}

#Contraste de proporciones
for(n1 in seq(190, 200, by=10))
for(x1 in seq(5, 20, by=5))
for(n2 in seq(90, 100, by=10))
for(x2 in seq(5, 10, by=3))
for( conf in seq(80, 99, by=5))
{
	preg<-paste("Para contrastar dos proporciones tomamos dos muestras ",
    "procedentes de dos poblaciones distintas de tamaños ", n1, " y ", n2,
    ", obteniendo los siguientes resultados:\\n\\n Éxitos en la primera muestra ",
    x1, "\\n\\nÉxitos en la segunda muestra ", x2, ".\\n\\n", sep="")


   #Contraste de medias
   preg1<-paste(preg, 
  "¿Podemos aceptar la igualdad de proporciones con un nivel de confianza del ",
  conf, "%?", sep="") 

	sol1<-c("Sí", "No", "Faltan datos")

	ttp<-prop.test(c(x1,x2), c(n1, n2), conf.level=conf/100)

	if(ttp$p.value>(1-conf/100))
	{
		pregunta(preg1, sol1)
	}else
	{
		pregunta(preg1, sol1[c(2,1,3)])
	}


   #p-valor
   preg2<-paste(preg, "El p-valor del contraste de igualdad de proporciones es:", sep="")
   sol2<-c(as.character(ttp$p.value), as.character(1-ttp$p.value), "Ninguno de los dos")
   pregunta(preg2, sol2)

}


#
#Práctica 6
#

#ANOVA
cat("\n\n$CATEGORY: Ex_Pr6_AOV\n\n")

for(media1 in seq(95, 100, by=5))
for(media2 in seq(95, 100, by=5))
for(media3 in seq(95, 100, by=5))
for(sigma in seq(1, 1, by=1))
{
	set.seed(media1*media2*media3*sigma)
	x1<-round(rnorm(3, media1, sigma), 2)
	set.seed(media1*media2*media3*sigma+1)
	x2<-round(rnorm(3, media2, sigma), 2)
	set.seed(media1*media2*media3*sigma+2)
	x3<-round(rnorm(3, media3, sigma), 2)

catalizador<-data.frame(catalizador=rep(c("C1","C2","C3"), each=3), 
   produccion=c(x1,x2,x3))

aovcat<-aov(produccion~catalizador, data=catalizador)
saovcat<-summary(aovcat)

preg<-"Un experimento se realizó para medir la producción proporcionada por cada uno de tres catalizadores en cierta reacción. El experimento se repitió tres veces. Las producciones en gramos son:\\n\\n"


	preg<-paste(preg, "Catalizador 1:", paste(x1, collapse=","), "\\n\\n")
	preg<-paste(preg, "Catalizador 2:", paste(x2, collapse=","), "\\n\\n")
	preg<-paste(preg, "Catalizador 3:", paste(x3, collapse=","), "\\n\\n")


	#CM ENTRE
	preg1<-paste(preg, "¿Cuál es el cuadrado medio de la variabilidad entre grupos?", sep="")
	sol1<-c(as.character(round(saovcat[[1]]$"Mean Sq", 4)), "Ninguno de los dos")
	pregunta(preg1, sol1)

	#CM DENTRO
	preg2<-paste(preg, "¿Cuál es el cuadrado medio de la variabilidad dentro de los grupos (o residual)?", sep="")
	sol2<-c(as.character(round(saovcat[[1]]$"Mean Sq"[2:1], 4)), "Ninguno de los dos")
	pregunta(preg2, sol2)

	#Estadístico de contraste
	preg3<-paste(preg, "¿Cuál es el valor del estadístico de contraste?", sep="")
	sol3<-round(c(saovcat[[1]]$"F value"[1], saovcat[[1]]$"Pr(>F)"[1], 
   saovcat[[1]]$"Df"[1]), 4)
	sol3<-as.character(sol3)
	pregunta(preg3, sol3)

	#p-valor
	preg4<-paste(preg, "¿Cuál es el p-valor del contraste?", sep="")
	sol4<-round(c(saovcat[[1]]$"Pr(>F)"[1], saovcat[[1]]$"F value"[1], 
   saovcat[[1]]$"Df"[1]), 4)
	sol4<-as.character(sol4)
	pregunta(preg4, sol4)


}



#Regresión lineal
cat("\n\n$CATEGORY: Ex_Pr6_RL\n\n")

lpred<-list(pvp=10000, emision=200, potencia=150, costecar=800, nplaza=3)

for(v in c("pvp", "emision", "potencia", "costecar", "nplaza"))
{
preg<-paste("Sea el fichero de datos cochesBig.sav. Considera el modelo lineal en el que usamos la variable ", v," para explicar el consumo. ", sep="")

	f1<-as.formula(paste("consumo", v, sep="~"))
	f2<-as.formula(paste(v, "consumo", sep="~"))

	m1<-lm(f1, data=cochesbig)
	m2<-lm(f2, data=cochesbig)

	#Interceptación
	preg1<-paste(preg, "¿Cuál es el valor del término independiente de ese modelo?", sep="")
sol1<-as.character(round(c(coef(m1)[[1]], coef(m1)[[2]], coef(m2)[[1]]), 4))
	pregunta(preg1, sol1)



	#Coeficiente
	preg2<-paste(preg, "¿Cuál es el valor del coeficiente de la variable explicativa de ese modelo?", sep="")
sol2<-as.character(round(c(coef(m1)[[2]], coef(m1)[[1]], coef(m2)[[2]]), 4))
	pregunta(preg2, sol2)


	#Predicción
	preg3<-paste(preg, "¿Cuál es el valor de consumo que me predice el modelo para un coche de ", v, " de ", lpred[v][[1]],"?", sep="")
	sol3<-as.character(round(c(predict(m1, as.data.frame(lpred[v]))[[1]],
   predict(m1, .5*as.data.frame(lpred[v]))[[1]]),4))
	sol3<-c(sol3, "Ninguno de los dos")
	pregunta(preg3, sol3)
}

#FIN
sink()


#Convertimos a UTF-8 (SÓLO funcionará en Linux) el fichero
system("iconv -f ISO-8859-1 -t UTF8 test_examen.txt>kk.txt")
system("cp kk.txt test_examen.txt")
system("rm kk.txt")




