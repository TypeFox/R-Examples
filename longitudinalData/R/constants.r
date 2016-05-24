cat("\n### Constantes ###\n")

MAX_CLUSTERS <- 26
CLUSTER_NAMES <- paste("c",1:MAX_CLUSTERS,sep="")


CRITERION_NAMES <- c("Calinski.Harabatz","Calinski.Harabatz2","Calinski.Harabatz3","Ray.Turi","Davies.Bouldin","BIC","BIC2","AIC","AICc","AICc2",#"entropie","ICL","ICL2",
                     "postProbaGlobal","random")

#CRITERION_MIN_OR_MAX<- c(Calinski.Harabatz=1,Kryszczuk.Calinski=1,Genolini.Calinski=1,Ray.Turi=-1,Davies.Bouldin=-1,random=1)

DISTANCE_METHODS <- c("manhattan", "euclidean", "minkowski", "maximum", "canberra", "binary")

CHOICE_STYLE <- list(
    typeTraj=c("l","l","n"),
    colTraj=c("clusters","black","black"),
    typeMean=c("b","b","b","b","l","l","n"),
    colMean=c("clusters","black","clusters","black","clusters","black","black"),
    pchMean=c("letters","letters","symbols","symbols","letters","letters","letters")
)

