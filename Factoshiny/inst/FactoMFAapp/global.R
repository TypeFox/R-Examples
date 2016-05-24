#global AFM
quanti=names(which(sapply(x,is.numeric)))
quali=names(which(!(sapply(x,is.numeric))))
VariableChoices=quanti
nom=rownames(x)
num=c(1:length(nom))
QualiChoice=quali
IdChoices=c(1:length(VariableChoices))
Idqualisup=c(1:length(QualiChoice))

