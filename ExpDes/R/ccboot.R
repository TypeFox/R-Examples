ccboot <-
function(y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE,
    main = NULL, B = 1000){
    trans.fator = function(matriz){
        trt <- matriz[,1]
        y <- matriz[,2]
        list(trt=trt,y=y)
    }
    anova.medias=function(dados){
        medias <- tapply(dados$y,list(dados$trt),mean)
    }
    qihb = function(dados,r,k,B){
        for(b in 1:B)
        {
            dadosb <- reamostra.dados(dados,r,k)
            anavab <- anova.dados(dadosb)
            mediasb <- anova.medias(dadosb)
            med.ordb <- sort(mediasb)
            if(b == 1){
                qihb <- qih.f(med.ordb,r,k,anavab)}
            else {
                qihb <- c(qihb,qih.f(med.ordb,r,k,anavab))}
            }
        qihb
    }
    dih.f = function(dad,r,k){
        D <- matrix(0,k,k)
        for (i in 1:(k-1))
        {
            seq <- (i+1):k
            for (j in seq)
            {
                D[i,j] <- dad[j]-dad[i]
            }
        }
        D    
    }
    hclstr = function(D,k,valores,qihb){
        v <- matrix(1:k,k,1)
        res <- matrix(1:(k-1),k-1,6)
        ct <- 1
        fmin  <- function(D,k)
        {
            min <- D[1,2];II<-1;JJ<-2
            for (ii in 1:(k-1))
            for (jj in (ii+1):k)
            {
                if (min>D[ii,jj])
                {
                    min <- D[ii,jj];II<-ii;JJ<-jj
                }
            }
            list(min=min,II=II,JJ=JJ)
        }
        arma.res = function(ct,min,II,JJ)
        {
            if ((v[II] > 0) & (v[JJ] > 0))
            {
                res[ct,2] <- min(v[II],v[JJ])
                res[ct,3] <- max(v[II],v[JJ])
                res[ct,4] <- min
            }
            if ((v[II] < 0) & (v[JJ] > 0))
            {
                L <- res[abs(v[II]),2]
                U <- res[abs(v[II]),3]
                min1 <- min(L,U);max1<-max(L,U) 
                res[ct,2] <- min(min1,v[JJ])
                res[ct,3] <- max(max1,v[JJ])
                res[ct,4] <- min
            }
            if ((v[II] > 0) & (v[JJ] < 0))
            {
                L <- res[abs(v[JJ]),2]
                U <- res[abs(v[JJ]),3]
                min1 <- min(L,U);max1<-max(L,U) 
                res[ct,2] <- min(min1,v[II])
                res[ct,3] <- max(max1,v[II])
                res[ct,4] <- min
            }
            if ((v[II] < 0) & (v[JJ] < 0))
            {
                L <- res[abs(v[II]),2]
                U <- res[abs(v[II]),3]
                min1 <- min(L,U);max1<-max(L,U)   
                L <- res[abs(v[JJ]),2]
                U <- res[abs(v[JJ]),3]
                min2 <- min(L,U);max2<-max(L,U) 
                res[ct,2] <- min(min1,min2)
                res[ct,3] <- max(max1,max2)
                res[ct,4] <- min
            }
        R <- res[ct,4]
        Q <- R/valores[3]
        P <- 1 - ptukey(Q,k,valores[2])
        res[ct,5] <- P
        valorpb <- (length(qihb[Q<qihb]))/B
        res[ct,6] <- valorpb
        res
    } 
    monta.res = function(D,k,II,JJ)
    {
        D1 <- matrix(0,k-1,k-1);iii<-2;jjj<-iii+1
        for (ii in 1:(k-1))
        for (jj in (ii+1):k)
        {
            if ((abs(ii-II) > 1e-6) & (abs(jj - JJ)>1e-6) & 
            (abs(ii-JJ) > 1e-6) & (abs(jj - II)>1e-6))
            {
                D1[iii,jjj] <- D[ii,jj]
                jjj <- jjj +1
                if (jjj >= k) 
                {
                    iii <- iii + 1; jjj <- iii +1  
                }  
            }
        }
        jj <- 2;
        for (ii in 1:k)
        {
            if ((abs(ii-II) > 1e-6) & (abs(ii-JJ) > 1e-6))
            {
                if ((ii>II) & (ii>JJ)) aux <- max(D[II,ii],D[JJ,ii])
                if ((ii>II) & (ii<JJ)) aux <- max(D[II,ii],D[jj,JJ]) 
                if ((ii<II) & (ii>JJ)) aux<-max(D[ii,II],D[JJ,ii])
                if ((ii<II) & (ii<JJ)) aux<-max(D[ii,II],D[ii,JJ])  
                D1[1,jj] <- aux
                jj <- jj + 1
            }          
        }
    D1          
    }
    monta.vet = function(v,II,JJ)
    {
        v1<-v[v != v[II]]
        v1<-v1[v1 != v[JJ]]
        v1<-c(-ct,v1)
        v1
    }
    Dct <- D;kt <- k 
    for (ct in 1:(k-1))
    {
        estmin <- fmin(Dct,kt)
        res <- arma.res(ct,estmin$min,estmin$II,estmin$JJ)
        Dct <- monta.res(Dct,kt,estmin$II,estmin$JJ)
        v <- monta.vet(v,estmin$II,estmin$JJ)
        kt <- kt - 1
    }
    res 
    }
    agrupar = function(cluster,k){
        grupos <- matrix(1:k,k,4)
        grupos[,1] <- seq(1:k)
        grupos[,2] <- grupos[,4] <- 0
        grupos[,3] <- 1
        i <- 1
        while((cluster[i,6] > alpha) & (i < k)){
            for(j in cluster[i,2]:cluster[i,3]){
            grupos[j,2] <- i
            grupos[j,4] <- 1
            }
            if(i != (k-1)) i <- i + 1
            else break
        }
        aux <- 1
        for(q in 1:k){
            if(q == 1)
            {
                grupos[q,3] <- aux
            }
            else{
                if((grupos[q,2] == grupos[q-1,2]) & (grupos[q,4] == 1))
                {
                    grupos[q,3] <- aux
                }
                if((grupos[q,2] == grupos[q-1,2]) & (grupos[q,4] == 0))
                {
                    aux <- aux + 1
                    grupos[q,3] <- aux
                }
                if((grupos[q,2] != grupos[q-1,2]) & (grupos[q,4] == 1))
                {
                    aux <- aux + 1
                    grupos[q,3] <- aux
                }
                if((grupos[q,2] != grupos[q-1,2]) & (grupos[q,4] == 0))
                {
                    aux <- aux + 1
                    grupos[q,3] <- aux
                }
            }
        }
        grupos
    } 
    desmascara.grupo = function(grupo,k,nomes){
        grupo[1:k,1] <- nomes[grupo[,1]]
        grupo
    }
    reamostra.dados = function(dados,r,k){
        x <- matrix(0,r*k,2)
        x[,1] <-  dados$trt
        numdados<-r*k
        unif <- trunc(runif(numdados)*numdados)+1
        x[,2] <- dados$y[unif]   
        list(trt=x[,1],y=x[,2])
    }
    anova.dados=function(dados){  
        dados$trt<-as.factor(dados$trt)  
        anovadados<-lm(dados$y ~ dados$trt)
        anovadad<-anova(anovadados)
        anovadad
    }
    qih.f = function(dad,r,k,anava){
        q <- (dad[k]-dad[1])/(anava$"Mean Sq"[2]/r)**0.5
    }
    k <- length(tapply(y, trt, length))
    r <- length(trt)/k
    mediast <- sort(tapply(y, trt, mean))  #ordem crescente
    nomes <- rownames(mediast)
    names(mediast) <- 1:length(mediast)
    dados <- cbind(trt,y)
    dados <- trans.fator(dados)
    qihb1 <- qihb(dados,r,k,B)
    D <- dih.f(mediast,r,k)
    for(i in 1:(k-1))
    {
        for(j in (i+1):k)
	  {
	      D[j,i] <- D[i,j]
	  }
    }
    S2fraiz <- cbind(SSerror/DFerror,DFerror,sqrt(SSerror/DFerror/r))
    valores <- as.vector(S2fraiz)
    names(valores) <- c("S2","f","raiz")
    cluster <- hclstr(D,k,valores,qihb1) 
    grupos <- agrupar(cluster,k)
    grupo.desm <- desmascara.grupo(grupos,k,nomes)
    agrupamento <- matrix(0,k,2)
    agrupamento[,1] <- grupo.desm[,1]
    agrupamento[,2] <- grupo.desm[,3]
    j = k
    group <- matrix(0,k,2)
    for(i in 1:k){
        group[i,1] <- agrupamento[j,1]
        group[i,2] <- agrupamento[j,2]
        j = j - 1
    }
    mediasd <- sort(mediast, decreasing = TRUE)
    names(mediasd) <- names(mediast)
    resultado <- data.frame(r = 0, f = group[,1], m = mediasd)
    res <- 1
    for (i in 1:(k-1)) {
        if (group[i,2] != group[i + 1,2]) {
            resultado$r[i] <- letters[res]
            res <- res + 1
            if (i == (length(resultado$r) - 1)) {
                resultado$r[i + 1] <- letters[res]
            }
        }
        else {
            resultado$r[i] <- letters[res]
            if (i == (length(resultado$r) - 1)) {
                resultado$r[i + 1] <- letters[res]
            }
        }
    }
    cat('\nBootstrap multiple comparison test\n------------------------------------------------------------------------\n')
    names(resultado) <- c("Groups", "Treatments", "Means")
    print(resultado)    
    cat("------------------------------------------------------------------------\n Note: Bootstrap multiple comparison method might present different results\n in each run since it is based on simulation.\n------------------------------------------------------------------------\n")
}
