# Fichier pour tester l'enchainement des fonctions sur des exemples simples
library("planor")
toto <- planor.factors(LETTERS[1:5],c(2,6,12,3,2))
tata <- planor.model(y~A+B+C+D+E)
titi <- planor.designkey(factors=toto, model=tata, nunits=8*9, base=~A+B)
tyty <- planor.design(pick(titi,c(1,1)))
titi <- planor.designkey(factors=toto, model=tata, nunits=8*9, base=~A+B+E)

 F2 <- planor.factors( factors=c(LETTERS[1:4], "bloc"), nlevels=c(6,6,4,2,6) )
 planor.harmonize(factors=F2[,1:5],
                  model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D,
                  base=~A+B+D)
