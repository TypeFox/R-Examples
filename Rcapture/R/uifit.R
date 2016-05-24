uifit <- function(x.closedp)
{

    ############################################################################################################################
    # Validation de l'argument fourni en entrée
    if(!any(class(x.closedp)=="closedp.t")) stop("'x.closedp' must be an object produced with 'closedp' or 'closedp.t")
    ############################################################################################################################
    # Ma fonction fonctionne correctement car les éléments glm et parameter de l'objet de type closedp
    # contiennent des sorties pour les mêmes modèles, et ce, dans le même ordre.

        t <- x.closedp$t

        ifirstcap <- NULL
        for (i in 1:t) { ifirstcap <- c(ifirstcap,rep(i,2^(t-i))) }

        # Identification des modèles qui ont été ajustés
        lmn<-rownames(x.closedp$results)
        nm<-length(lmn)

        # Initialisation du tableau principal de sorties
        tableau <- matrix(nrow=t+5,ncol=nm+1)
        dimnames(tableau) <- list(paste("u",1:(t+5),sep = ""),c("observed",lmn))
        stat <- matrix(nrow=nm,ncol=1)
        dimnames(stat) <- list(lmn,c("Chi-suare value"))

        # Valeurs observées
        desc<- descriptive(x.closedp$X,x.closedp$dfreq)
        tableau[,1]<-c(desc$base.freq[,2],rep(NA,5))

    # Boucle qui traite tous les modèles
    for (j in 1:nm)
    {
        glmo <- x.closedp$glm[[j]]
        N <- x.closedp$parameters[[j]][1,1]
        if (lmn[j]=="M0")
        {
            p <- exp(glmo$coef[2])/(1+exp(glmo$coef[2]))
            tableau[,j+1] <- N*p*(1-p)^(0:(t+4))
        } else
        if (lmn[j]=="Mb")
        {
            p <- 1-exp(glmo$coef[2])/(1+exp(glmo$coef[3]))
            tableau[,j+1] <- N*p*(1-p)^(0:(t+4))
        } else
        if (lmn[j]=="Mh Poisson2")
        {
            EprobaP_general <- function(x,beta,tau,a,t,k){
            (exp(beta)*(1+exp(beta)*a^x)^(t-k))*(a*tau)^x/(factorial(x)*sum(na.rm=TRUE,choose(t,0:t)*exp(beta*(0:t)+tau*a^(0:t))))
            }
            value_Eproba <- rep(0,t+5)
            for (i in 1:(t+5))
            {
                EprobaP <- function(x){ EprobaP_general(x,glmo$coef[2],glmo$coef[3],2,t,i)}
                value_Eproba[i] <- sum(na.rm=TRUE,EprobaP(0:100))
            }
            tableau[,j+1] <- N*value_Eproba
        } else
        if (lmn[j]=="Mh Darroch")
        {
            if(glmo$coef[3]>0)
            {
                EprobaD_general <- function(x,beta,tau,t,k){
                exp(-(x^2)/(2*tau))*((1+exp(beta+x))^(t-k))*exp(beta+x)/
                (sqrt(2*pi*tau)*sum(na.rm=TRUE,choose(t,0:t)*exp(beta*(0:t)+tau*((0:t)^2)/2)))
                }
                value_Eproba <- rep(0,t+5)
                for (i in 1:(t+5))
                {
                    EprobaD <- function(x){EprobaD_general(x,glmo$coef[2],glmo$coef[3],t,i)}
                    value_Eproba[i] <- integrate(EprobaD,-100,100)$value
                }
                tableau[,j+1] <- N*value_Eproba
            } else {
                for ( i in 1:t ) { tableau[i,j+1] <- sum(na.rm=TRUE,glmo$fitted.values[ifirstcap==i]) }
            }
        } else {
            for ( i in 1:t ) { tableau[i,j+1] <- sum(na.rm=TRUE,glmo$fitted.values[ifirstcap==i]) }
        }


        # Stat d'ajustement du chi-deux
        stat[j,1] <- sum(na.rm=TRUE,((tableau[1:t,1]-tableau[1:t,j+1])^2)/tableau[1:t,j+1])

    }

    # Statistiques sur le jour de la première capture
    Mean <- colSums((1:t)*tableau[1:t,])/colSums(tableau[1:t,])
    Variance <- colSums(((1:t)^2)*tableau[1:t,])/colSums(tableau[1:t,]) - Mean^2
    firstcapt <- cbind(Mean,Variance)

    # output
    list(predicted=tableau,fit.stat=stat,day.first.capt=firstcapt)
}
