`adonis.pertables` <-
function(formula=X~., data, permutations = 5, method = "bray"){

require(vegan)
        
fmla<-fmla0<- formula
fmla[[2]]<-substitute(X)
Xsim<-eval(fmla0[[2]])[[2]] # selecciona el segundo elemento de un objeto pertables, que es la lista de tablas permutadas
Xraw<-eval(fmla0[[2]])[[3]] # selecciona el tercer elemento, que es la raw data table
    
# calculo de adonis sobre los raw data
X<-Xraw
adonis.raw<-adonis(formula=fmla, data=data, permutations = permutations, method = method)
namestab<-rownames(adonis.raw$aov.tab)

# c<U+00E1>lculo de adonis sobre los datos simulados
X<-Xsim
results<- lapply(X, function(X)adonis(formula=fmla, data=data, permutations=permutations,method=method))

F<- sapply(results, function(x) x$aov.tab$F.Model)
rownames(F) <- namestab

R2<- sapply(results, function(x) x$aov.tab$R2)
rownames(R2) <- namestab

p<- sapply(results, function(x) x$aov.tab$Pr)
rownames(p) <- namestab
R2.quant<- apply(R2, 1, quantile,c(0,0.005,0.025,0.5,0.975,0.995,1),na.rm=TRUE)[,1:(length(namestab)-2)] #quitamos las columnas de total y residuales
p.quant<- apply(p, 1, quantile,c(0,0.005,0.025,0.5,0.975,0.995,1),na.rm=TRUE)[,1:(length(namestab)-2)] #quitamos las columnas de total y residuales



F.raw<-adonis.raw$aov.tab$F.Model
Fls <- cbind(F.raw, F)
Fls <- Fls[-c((length(F.raw)-1):length(F.raw)), ]
ptax<- apply(Fls, 1, function(x) (rank(x)/(length(x)))[1])
ptax<-ifelse(ptax<=0.5,ptax,1-ptax)
adonis.raw$aov.tab$Prtax <- c(ptax, NA, NA)
names(adonis.raw$aov.tab)[7] <-"Pr(tax)"
adonis.raw$call <- match.call()
adonis.output<- list(raw=adonis.raw, simulation=list(F=F, R2 = R2, pvalue = p, R2.quant = R2.quant, p.quant = p.quant))
class(adonis.output) <- c("adonis.pertables",class(adonis.output))
return(adonis.output)
}

