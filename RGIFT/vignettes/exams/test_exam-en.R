#
#Produce questions for a exam
#
#Run this sript using echo("test_exam-en.R") to avoid
#garbage in the output  file
#

library(RGIFT)



#Redefine function to set default values. This has been done to
#simplify the modifications in the original file
pregunta<-function(qtxt, atxt){GIFTMC(qtxt, atxt, wwrong="-33.333")}


#Set file where the questions are saved
sink("test_exam.txt")

#
#Computer Lab 1: Intro to R
#
#cat("\n\n$CATEGORY: Ex_Pr1\n\n")
GIFTcategory("Ex_Pr1")

for(i in 1:5)#5 different data sets 
{

	set.seed(1*i)
	x<-round(rnorm(10,10,1), 2)


     preg<-paste("Consider the following vector:\\n\\n", deparse(x, width.cutoff=100L),
   "\\n\\n", sep="")

	#Mean
	preg1<-paste(preg, "Its mean is:")
	sol1<-as.character(mean(x)*c(1,.9, 1.1))
	pregunta(preg1, sol1)

	#Standard deviation
	preg2<-paste(preg, "Its standard deviation is:")
	sol2<-as.character(sd(x)*c(1,.5, 1.5))
	pregunta(preg2, sol2)

	#Variance
	preg2b<-paste(preg, "Its variance is:")
	sol2b<-as.character(var(x)*c(1, .5, 1.5))
	pregunta(preg2b, sol2b)

	#0'40 Quantile
	preg3<-paste(preg, "Its 0.40 quantile is:")
	sol3<-as.character(quantile(x, .4)*c(1,.75,1.25))
	pregunta(preg3, sol3)


	#0'90 Quantile
	preg4<-paste(preg, "Its 0.90 quantile is:")
	sol4<-as.character(quantile(x, .9)*c(1,.75,1.25))
	pregunta(preg4, sol4)

}


#Some questions about specific commands
pregunta("In order to know the length of a vector in R, you will use:",
   c("length()", "long()", "None of the above")
)

pregunta("To make a bar plot you will use:",
   c("barplot()", "plot()", "None of the above")
)

pregunta("To make a histogram you will use:",
   c("None of the above","histo()", "h()")
)

pregunta("The command to compute the square root of 13 is:",
   c("sqrt(13)", "sd(13)", "str(13)")
)

pregunta("To compute factorial of 5, you will use:",
   c("factorial(5)", "facto(5)", "f(5)")
)

pregunta("To make a sequence of 10 values between 0 and 10, you will write:",
   c("seq(0,5,length\\=10)", "0:5", "sequ(0,5)" )
)

pregunta("To make a pie chart, what function will you use?:",
   c("None of the above", "sector()", "tarta()")
)

pregunta("A pie chart is done for ",
  c("a category-based variable", 
     "a continuous variable", "None of the above")
)


#
#Computer Lab 2: Simulation with R (I)
#
#cat("\n\n$CATEGORY: Ex_Pr2\n\n")
GIFTcategory("Ex_Pr2")


for(i in 1:5)#5 different vectors of data
{

        set.seed(10*i)
        x<-round(100*rnorm(8,10,1))/100


        preg<-paste("Consider the following vector:\\n", deparse(x, width.cutoff=100L),
   "\\n", sep="")

        #
        preg1<-paste(preg, "How much is the sum of all the values that are lower than the mean?")
        sol1<-as.character( sum(x[x<mean(x)])*c(1,.75,.5))
        pregunta(preg1, sol1)

	#
        preg2<-paste(preg, "How much is the sum of all the values that are greater than the mean?")
        sol2<-as.character( sum(x[x>mean(x)])*c(1,1.75,1.5))
        pregunta(preg2, sol2)

        #
        preg3<-paste(preg, "How much is the sum of all the values that are lower than the median?")
        sol3<-as.character( sum(x[x<median(x)])*c(1,.75,.5))
        pregunta(preg3, sol3)

	#
        preg4<-paste(preg, "How much is the sum of all the values that are greater than the median?")
        sol4<-as.character( sum(x[x>median(x)])*c(1,.75,1.5))
        pregunta(preg4, sol4)

}


pregunta("In order to compute a matrix determinat you will use:",
  c("det(X)", "deter(X)", "determinate(X)")
)

pregunta("To compute a matrix inverse you will use:",
   c("solve(X)", "invers(X)", "1/X")
)

pregunta("To perform a matrix multiplication on two matrices, you will use:",
   c("X%*%Y", "X*Y", "XY")
)

pregunta("To multiply two matrices element-wise, you will use:",
   c("X*Y", "X%*%Y", "X.Y")
)

pregunta("To obtain a sample of size 10 without replacement from a vector X,
you will use:",
   c("sample(X,10)", "samples(X,10,TRUE)", "Ninguna de las dos")
)

pregunta("To obtain a sample of size 3 with replacement from a vector X,
you will use:",
   c("sample(X,3,rep\\=TRUE)", "sample(X,3,FALSE)", "Ninguna de las dos")
)


pregunta("What command would you use to sample 100 times the number of 
jack, queens and kings in a deck of 52 cards when 10 cards are taken
with replacement?",
   c("rbinom(100,10,12/40)", "rbinom(100,10,10/40)", "Ninguna de las dos")
)


#Some questions about probability distributions
nombres<-c("Binomial", "Normal", "Exponential", "Uniform") 
distr<-c("binom", "norm", "exp", "unif")
params<-list(c(10, .27), c(4, 2), c(5), c(0, 15))

ndistr<-length(nombres)


for(i in 1:ndistr)
{
	preg<-paste("Let X be a ", nombres[i], " random variable with parameters ")
	preg<-paste(preg, paste(params[[i]], collapse=" y "), ". ", sep="")


	#Probabilidad acumulada
	preg1<-paste(preg, "In order to compute the probability of X being lower or equal than 3, you will use:", sep="")
	sol1<-paste("p",distr[i], "(3,", paste(params[[i]], collapse= ", "), ")", sep="") 
	sol1<-c(sol1, paste("1-", sol1,sep=""), "None of the above")
	
	pregunta(preg1, sol1)

	#1-Probabilidad acumulada
	 preg2<-paste(preg, "In order to compute the probability of X being higher than than 3, you will use:", sep="")
	sol2<-paste("p",distr[i], "(3,", paste(params[[i]], collapse= ", "), ")", sep="") 
	sol2<-c(paste("1-", sol2,sep=""), sol2, "None of the above")
	
	pregunta(preg2, sol2)

	#Quantiles
	qnt<-sample(1:100,1)
	preg3<-paste(preg, "To compute the ",qnt, "% percentile, you will use:", sep="")
	sol3a<-paste("q",distr[i], "(", qnt/100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	sol3b<-paste("r",distr[i], "(", qnt/100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	sol3c<-paste("p",distr[i], "(", qnt/100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	
	pregunta(preg3, c(sol3a, sol3b, sol3c))

	#Simulation
	preg4<-paste(preg, "To sample 100 values from than random variable, you will use:", sep="")
	sol4a<-paste("r",distr[i], "(", 100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	sol4b<-paste("q",distr[i], "(", 100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	sol4c<-paste("p",distr[i], "(", 100,",",paste(params[[i]], collapse= ", "), ")", sep="") 
	
	pregunta(preg4, c(sol4a, sol4b, sol4c))
}



#
#Computer Lab 3:  Simulation with R (II)
#

#cat("\n\n$CATEGORY: Ex_Pr3\n\n")
GIFTcategory("Ex_Pr3")


pregunta("Let X be a vector of type factor. In order to know the number of observations in each level, you will use:",
   c("table(X)", "levels(X)", "None of the above")
)

pregunta("In order to create a factor made of letters A, B and C where each letter repeats 10, 20 and 50 times, respectively, you will use:",
  c("factor(rep(c(\"A\",\"B\",\"C\"),c(10,20,30)))", 
    "factor(rep(c(\"A\",\"B\",\"C\"),10,20,50))", 
    "None of the above")
)

pregunta("To convert a vector X into a factor:",
   c("as.factor(X)", "fact(X)", "None of the above")
)

pregunta("To obtain a numeric summary of a data.frame X, you will use:",
   c("summary(X)", "Summary(X)", "None of the above")
)

pregunta("To show the names of the variable in data.frame X, you will use:",
   c("names(X)", "variables(X)", "None of the above")
)

#Questions about data set CochesBig.sav (in Spanish)
library(foreign)
cochesbig<-read.spss("CochesBig.sav", to.data.frame = TRUE)

variables<-names(cochesbig)#List of all variable names
marcas<-levels(cochesbig$marca)#Manufaturers/brands

for(v in variables[1:3])
for(m in marcas[1:5])
{
	#Media
	preg1<-paste("Using the CochesBig.sav data set. The mean of variable ", v, " for manufacturer ", sep="")
	preg1<-paste(preg1, m, " is", sep="")
	sol1<-as.character(round(mean(cochesbig[cochesbig$marca==m,v])*c(1, .6, .5), 4))
	pregunta(preg1, sol1)


	#Mediana
	preg2<-paste("Using the CochesBig.sav data set. The mean of variable ", v, " for manufacturer ", sep="")
	preg2<-paste(preg2, m, " is", sep="")
	sol2<-as.character(round(quantile(cochesbig[cochesbig$marca==m, v],.5)*c(1, .6, .5), 4))
	pregunta(preg2, sol2)

	#Cuantil del 0.6
	preg3<-paste("Using the CochesBig.sav data set. The 60% percentile of variable ", v, " for manufacturer ", sep="")
	preg3<-paste(preg3, m, " is", sep="")
	sol3<-as.character(round(quantile(cochesbig[cochesbig$marca==m, v],.60)*c(1, .6, .5), 4))
	pregunta(preg3, sol3)

}


#
#Computer Lab 4: Central Limit Theorem
#


#cat("\n\n$CATEGORY: Ex_Pr4\n\n")
GIFTcategory("Ex_Pr4")

pregunta("The distribution of the mean of 10000 independent and equally distributed random variables...",
   c("... is Normal, according to the Central Limit Theorem.",
   "... depends on the distribution of the random variables.",
   "None of the above.")
)

pregunta("The distribution of the mean of 10000 independent and equally distributed random variables...",
   c("None of the above",
   "... is Normal in most cases",
   "... depends on the distribution of the random variables.")
)


for(m in seq(10, 25, by=5))
{
for(ss in seq(30, 50, by=10))
{
	preg1<-paste("Given a Normal population with mean ", m, ", if I compute 1000 95% confidence intervals with samples of size ", ss, " from that population:", sep="")
	sol1<-c("None of the above", paste("It is not possible to know how many intervals will contain value ", m, "", sep=""), paste("More than 95% of the intervals will contain value ", m, sep=""))
	pregunta(preg1, sol1)


	for(dt in seq(20, 50, by=10))
	{
preg2<-paste("Let X_i be a Normal random variable with mean ", m, " and standard deviation ", dt, ". What is the standard deviation of the mean of a sample of size ", ss, "?", sep="")
sol2<-c(paste(dt, "/sqrt(", ss,")", sep=""),
   paste(dt, ", same as X_i", sep=""),
   "None of the above")
pregunta(preg2, sol2)

preg3<-paste("Let X_i be a Normal random variable with mean ", m, " and standard deviation ", dt, ". What is the mean of the average of a sample of size ", ss, "?", sep="")
sol3<-c(paste(m, ", same as X_i", sep=""),
paste(dt, "/sqrt(", ss,")", sep=""),
   "None of the above")
pregunta(preg3, sol3)


	}

}
}



#COnfidence intervals and tests based on  cochesbig
lm0<-as.list(apply(as.matrix(cochesbig[,1:3]), 2, mean))


for( conf in seq(80, 99, by=5))
{
for(v in variables[1:3])
{
	m0<-round(lm0[[v]]*runif(1, .75,1.25))

	preg1<-paste("Using the CochesBig.sav data set. The mean of ", v,
  " can be ", m0, " with a ", conf,"% confidence.", sep="")

	sol1<-c(paste("Yes, becuase that value is inside the ",conf,"% confidence interval", sep=""),
   paste("No, because the mean is ",round(mean(cochesbig[,v]), 2), sep=""),
   "None of the above")

	tt<-t.test(cochesbig[,v], mu=m0, conf.level=conf/100)

	if(tt$p.value>(1-conf/100))
	{
		pregunta(preg1, sol1)
	}else
	{
		pregunta(preg1, sol1[c(3,1,2)])

	}

	preg2<-paste("Using the CochesBig.sav data set. A ", conf,   "% confidence interval for variable ", v, "is:", sep="")
 
	sol2<-c(paste("(", paste(round(as.numeric(tt$conf.int), 4), collapse=", "), ")", sep=""),
   paste("(", paste(5+round(.9*as.numeric(tt$conf.int), 4), collapse=", "), ")", sep=""),
   "None of the above")

}
}



#Tests of proportions

for( conf in seq(80, 99, by=5))
{
for(n in seq(180, 200, by=5))
{
for(x in seq(5, 10, by=5))
{
	p0<-round((x/n)*runif(1,.75,5), 2)

	preg1<-paste("In a factory, ", n, 
   " pieces were inspected and ", x, " of them were found to be defective. ",
   " With a confidence of ", conf, "%, is there any evidence that",
   " the proportion of defective pieces is equal to ", p0,"?", sep="")

	ttp<-prop.test(x,n, p=p0, conf.level=conf/100)

	sol1<-c(paste("Yes, because that value is inside the confidence interval.",sep=""),
   paste("No because the observed proportion is ", round(x/n, 4), sep=""),
   "None of the above")

	if( ttp$p.value>(1-conf/100) )
	{
		pregunta(preg1, sol1)
	}else
	{
		pregunta(preg1, sol1[c(3,2,1)])

	}


	preg2<-paste("In a factory ", n,
   " pieces were inspected and ", x, " of them were found to be defective. ",
"A ", conf,   "% confidence interval for the proportion of defective pieces is", sep="")

        sol2<-c(paste("(", paste(round(as.numeric(ttp$conf.int), 4), collapse=", "), ")", sep=""),
   paste("(", paste(.15+round(1.1*as.numeric(ttp$conf.int), 4), collapse=", "), ")", sep=""),
   "None of the above")
	pregunta(preg2, sol2)


}
}
}



#
#Computer Lab 5:
#


#cat("\n\n$CATEGORY: Ex_Pr5\n\n")
GIFTcategory("Ex_Pr5")

#Test on the means of two populations
for(conf in seq(80, 99, by=5))
for(m1 in seq(100, 110, by=5))
for(m2 in seq(90, 100, by=5))
for(sigma in seq(5, 5, by=1))
{
	set.seed(conf*m1*m2*sigma)
	hombres<-round(rnorm(10, m1, sigma), 2)
	set.seed(conf*m1*m2*sigma+1)
	mujeres<-round(rnorm(10, m2, sigma), 2)

	preg<-paste("Consider the following data:\\n\\n", 
    "Girls' IQ: ",
    deparse(mujeres, width.cutoff=100L), "\\n\\n",
    "Boys' IQ: ",
    deparse(hombres, width.cutoff=100L), "\\n\\n", sep="")


   #Contraste
    preg1<-paste(preg, 
   "Could we say that the mean IQ for boys and girls is the same? ",
    "Consider a significance level of ", 1-conf/100, ".", sep="")


   sol1<-c("Yes", "No", "We need more data")

   tt<-t.test(mujeres, hombres, conf.level=conf/100)

   if(tt$p.value>(1-conf/100))
   {
	pregunta(preg1, sol1)
   }else
   {
	pregunta(preg1, sol1[c(2,1,3)])
   }



   #p-valor
   preg2<-paste(preg, "The p-value associated to the test is:", sep="")
   sol2<-c(as.character(tt$p.value), as.character(1-tt$p.value), "None of the above")
   pregunta(preg2, sol2)

}

#Test for rates
for(n1 in seq(190, 200, by=10))
for(x1 in seq(5, 20, by=5))
for(n2 in seq(90, 100, by=10))
for(x2 in seq(5, 10, by=3))
for( conf in seq(80, 99, by=5))
{
	preg<-paste("In order to test two rates  two samples ",
    "are taken from two different populations. Their sample sizes are ", n1, " and ", n2,
    ", and the following values were obtained:\\n\\n Number of successes in the first sample: ",
    x1, "\\n\\nNumber of sucesses in the second smaple: ", x2, ".\\n\\n", sep="")


   #Contraste de medias
   preg1<-paste(preg, 
  "Could we say that these two rates are equal with a confidence level of  ",
  conf, "%?", sep="") 

	sol1<-c("Yes", "No", "More data are needed.")

	ttp<-prop.test(c(x1,x2), c(n1, n2), conf.level=conf/100)

	if(ttp$p.value>(1-conf/100))
	{
		pregunta(preg1, sol1)
	}else
	{
		pregunta(preg1, sol1[c(2,1,3)])
	}


   #p-valor
   preg2<-paste(preg, "The p-value associated with the test is:", sep="")
   sol2<-c(as.character(ttp$p.value), as.character(1-ttp$p.value), "None of the above")
   pregunta(preg2, sol2)

}


#
#Computer Lab 6: Linear regression and ANOVA
#

#ANOVA
#cat("\n\n$CATEGORY: Ex_Pr6_AOV\n\n")
GIFTcategory("Ex_PR6_AOV")


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

preg<-"An experiment was made to measure the production related to three different catalyzers used in an industrial process. The experiment was repeated three times for each catalyzer. Productions in grams were:\\n\\n"


	preg<-paste(preg, "Catalyzer 1:", paste(x1, collapse=","), "\\n\\n")
	preg<-paste(preg, "Catalyzer 2:", paste(x2, collapse=","), "\\n\\n")
	preg<-paste(preg, "Catalyzer 3:", paste(x3, collapse=","), "\\n\\n")


	#CM ENTRE
	preg1<-paste(preg, "The Mean Square of the variability between groups is:", sep="")
	sol1<-c(as.character(round(saovcat[[1]]$"Mean Sq", 4)), "None of the above")
	pregunta(preg1, sol1)

	#CM DENTRO
	preg2<-paste(preg, "The Mean Square of the variability within groups is:", sep="")
	sol2<-c(as.character(round(saovcat[[1]]$"Mean Sq"[2:1], 4)), "None of the above")
	pregunta(preg2, sol2)

	#Estadístico de contraste
	preg3<-paste(preg, "The value of the test statistic is:", sep="")
	sol3<-round(c(saovcat[[1]]$"F value"[1], saovcat[[1]]$"Pr(>F)"[1], 
   saovcat[[1]]$"Df"[1]), 4)
	sol3<-as.character(sol3)
	pregunta(preg3, sol3)

	#p-valor
	preg4<-paste(preg, "The p-value of the test is:", sep="")
	sol4<-round(c(saovcat[[1]]$"Pr(>F)"[1], saovcat[[1]]$"F value"[1], 
   saovcat[[1]]$"Df"[1]), 4)
	sol4<-as.character(sol4)
	pregunta(preg4, sol4)
}



#Linear regression
#cat("\n\n$CATEGORY: Ex_Pr6_RL\n\n")
GIFTcategory("Ex_Pr6_LR")


lpred<-list(pvp=10000, emision=200, potencia=150, costecar=800, nplaza=3)

for(v in c("pvp", "emision", "potencia", "costecar", "nplaza"))
{
preg<-paste("Using the cochesBig.sav data set, consider the linear model where variable  ", v," is used to explain 'consumo'. ", sep="")

	f1<-as.formula(paste("consumo", v, sep="~"))
	f2<-as.formula(paste(v, "consumo", sep="~"))

	m1<-lm(f1, data=cochesbig)
	m2<-lm(f2, data=cochesbig)

	#Intercept
	preg1<-paste(preg, "What is the value of the intercept of the model?", sep="")
sol1<-as.character(round(c(coef(m1)[[1]], coef(m1)[[2]], coef(m2)[[1]]), 4))
	pregunta(preg1, sol1)

	#Coefficient 
	preg2<-paste(preg, "What's the value of the regression coefficient of the predictive variable of the model?", sep="")
sol2<-as.character(round(c(coef(m1)[[2]], coef(m1)[[1]], coef(m2)[[2]]), 4))
	pregunta(preg2, sol2)


	#Prediction
	preg3<-paste(preg, "What is the value of 'consumo' for a car with a value of ",  
        v, " of ", lpred[v][[1]],"?", sep="")
	sol3<-as.character(round(c(predict(m1, as.data.frame(lpred[v]))[[1]],
   predict(m1, .5*as.data.frame(lpred[v]))[[1]]),4))
	sol3<-c(sol3, "None of the above")
	pregunta(preg3, sol3)
}

#THE END
sink()


#Convert file to UTF-8 (this will only work on LInux)
#system("iconv -f ISO-8859-1 -t UTF8 test_examen.txt>kk.txt")
#system("cp kk.txt test_examen.txt")
#system("rm kk.txt")

