############################## Magic function to make the matrix start with constrain ###############################################
constrainAB <- function(cA)
{
 # The matrices cA contain the constraint values of the mixing parameters. 
 # For fixed value constraint use the values at respective location in the matrix.
 # For aki, all row sums must be less or equal 1. 
 # For equality value constraint use whole numbers starting form 2. 
 # Same numbers at different locations of the matrix show equal parameters. 
 # USE "NA" TO FILL UP THE REST OF THE MATRICES. 

 #----------------------------------------------------------------------------
 #checking cA
 a  <- a1 <- cA

 I <- nrow(cA)
 K <- ncol(cA) 
 if(any(apply(cA,
              1,
              function(x) any(sum(x[x < 1],
                                  na.rm = TRUE) > 1))) == TRUE) {

  stop("All row sums of the fixed constraints mixing parameter matrix must be 
       less or equal 1")

 }

 if(any(apply(cA,
              1, 
              function(x) length(x[!is.na(x) & x < 1])) == (K-1))) {

  stop("Any row of the fixed constraints mixing parameter matrix with K-1 
       elements must be completed so that it will add up to one")

 }

 colnames(a) <- 1:K

 aux1 <- matrix(a[which(apply(a,
                              1,
                              function(x) x%*%rep(1,K)) != 'NA'),],
                ncol = K) # Extraindo da matrix as linhas completas (sem NA). 

 rownames(aux1) <- which(apply(a,
                               1,
                               function(x) x%*%rep(1,K)) != 'NA')

 if(all(aux1 < 1)){

  aux1 <- matrix(logical(0),
                 ncol = K)

 }


 aux2 <- lapply(apply(a,
                      1,
                      function(x) x[x < 1] ),
                function(x)x[!is.na(x)])# Extraindo da matrix as linhas completas com somente fixed constraint (caso tenha é claro!).

 names(aux2) <- 1:I

 if(dim(aux1)[1] > 1){

  aux3 <- combn(1:dim(aux1)[1],
                2) # Objeto auxiliar para computar se as linhas completas são idênticas.

  aux4 <- sapply(1:dim(aux3)[2],
                 function(i) setequal(aux1[aux3[1,i],],
                                      aux1[aux3[2,i],])) # As linhas completas possuem os mesmos equalitys?

  aux5 <- sapply(1:dim(aux3)[2],
                 function(i) identical(ftable(aux1[aux3[1,i],]),
                                       ftable(aux1[aux3[2,i],]))) # As linhas completas possuem as mesmas frequências?

  if(dim(aux1)[1] > 1 & aux4 == TRUE & aux5 == FALSE & length(unlist(lapply(aux2,
                                                                            function(x) x[!is.na(x)]))) == 0) {

   stop("Your matrix has more than one column with same equality constrain, therefore, this equalitys need to be with same proportions for both sum 1.")

  }

  if(any(apply(aux1,
               1,
               function(x)any(x<1)) == TRUE)){
   aux51 <- apply(aux1,
                  1,
                  function(x) x[any(x <1)])

   aux52 <- aux51[sapply(aux51,
                         function(x) length(x) != 0)]

   aux53 <- unlist(lapply(aux52,
                          function(x) x[x>1]))

   aux54 <- unique(aux53)

   if(length(aux54) == 1){
    stop('Any row of your matrix with fixed and equality constrain can be completed by yourself')
   }  
  }
 }


 if(any(table(matrix(a[a>1],
                     nrow=1)) == 1)){

  stop("The equalitys must be frequencies greater or equal two.")

 } 

 #---------------------------- BUILDING cA -------------------------------# 

 maxca <- max(cA, 
              na.rm=TRUE)

 #    --------          -----------         ---------       --------   ---- #
 # If TRUE there are equality parameters find the indices #
 #    --------          -----------         ---------       --------   ---- # 
 if(maxca > 1) { 

  al <- list()

  for(i in 2:max(a, 
                 na.rm=TRUE)){

   al[[i]] <- which(a==i, 
                    arr.ind=TRUE)

   al[[i]] <- al[[i]][order(al[[i]][,1],
                            al[[i]][,2]),]
  }

  al <- al[!sapply(al,
                   is.null)]
  al <- lapply(al, 
               function(x) matrix(x, 
                                  ncol = 2))
  names(al) <- 2:max(a,
                     na.rm=T)
  # ------- --------- --------- ---------- ---------- --------- ------- - #
  # A partir de agora, iremos buscar na matriz de equality constraints se #
  # tem alguma linha toda preenchida,  pois existe  uma  necessidade   de #
  # escrever um programa específico para estes casos.                     #
  # -------- -------- --------- -------- -------- -------- --------- ---  #

  # Este objeto checa na matriz se tem alguma linha toda preenchida, ou seja, sem NA.
  if(dim(aux1)[1] < 1){

   al <- list()

   for(i in 2:max(a, 
                  na.rm=TRUE)){

    al[[i]] <- which(a==i, 
                     arr.ind=TRUE)

    al[[i]] <- al[[i]][order(al[[i]][,1],
                             al[[i]][,2]),]

    a[al[[i]]] <- runif(1)/K
   }
  } else {

   # Vamos começar o árduo trabalho no caso em que temos uma ou mais 
   # linhas todas preenchidas.

   posequalaux <- apply(aux1,
                        1,
                        function(x) x[x > 1]) # Tirando os fixed!

   if(is.list(posequalaux)){

    posequal <-  posequalaux[sapply(posequalaux, function(x) length(x)!=0)]  

   } else 

    if(is.vector(posequalaux)) {

     posequal <- list(posequalaux)
     names(posequal) <- rownames(aux1)  
    } else {

     posequal <- list()

     for(i in 1:dim(posequalaux)[2]){

      posequal[[i]] <- posequalaux[,i]

     }  
     names(posequal) <- rownames(aux1)    
    }

   # ------ ------ ------- ------- -------- -------- -------- ------ --#
   # Quando os equalitys não são os mesmos, precisamos saber quais são #
   # estes equality. Neste caso,  precisamos    saber  primeiro, quais #
   # aperecem com menor frequência, pois estes, determinarão os valores#
   # dos demais. Este objeto retorna isto.                             #
   # --------- ----- ------ ------ ------ ----- ------ ------ ----- -- #

   auxEE <- lapply(posequal, 
                   function(x) names(table(x))[table(x) != min(table(x))]) # Este objeto checa se os equalitys estão na mesma proporção.

   aux6 <- list()

   for(i in 1:length(auxEE)){
    aux6[[i]] <-  ifelse(length(auxEE[[i]]) == 0,
                         names(table(posequal[[i]]))[1], # Equalitys na mesma proporção
                         names(table(posequal[[i]]))[table(posequal[[i]]) == min(table(posequal[[i]]))]) # Os equalitys com menores frequências em cada linha.

   }

   # Este objeto retorna a frequência absoluta de cada equality constrain.
   aux7  <-  lapply(posequal, 
                    table)

   # Este objeto retorna a frequência absoluta apenas dos equality com menores frequências. 
   aux8 <-  mapply(function(x,y) x[y], 
                   aux7, 
                   aux6,
                   SIMPLIFY = FALSE)

   names(aux6) <- names(aux8)
   # ------ ------ ----- ------ ----- ------ ----- -------- ------ --#
   # Agora a coisa começa a ficar um pouco mais complicado. O valor  #
   # aleatório dos equalitys na primeira linha identificada completa #
   # (sem NA) irá determinar os valores das demais linhas caso tenha #
   # equalitys iguais.
   # ----- ------- ------ ------ ------ ----- ------ ---- ------ -- -#

   posequaldifaux <- sapply(lapply(posequal[-1], 
                                   function(x) x%in%posequal[[1]]), 
                            function(x) all(x == FALSE)) # Este objeto retorna os valores aleatórios do(s) equality com menor(es) frequência(s) apenas das linhas completas que os equalitys são todos diferentes, ou seja, não tem nenhum equality em comum.  

   posequaldif <- names(posequaldifaux[posequaldifaux == TRUE])

   namesvaluesfr <- c(names(posequal[1]),
                      posequaldif)

   valuesfraux <- list(list())

   for(i in namesvaluesfr){

    for(j in 1:length(aux6[[i]])){

     valuesfraux[[i]][[j]] <- rep(runif(1)/K,
                                  aux8[[i]][j])

    }

    names(valuesfraux[[i]]) <- aux6[[i]]

   }

   # Este objeto retorna os equality com as demais frequências das linhas completas com equalitys totalmente diferentes.

   aux9 <- list()

   for(i in namesvaluesfr){

    aux9[[i]] <- names(table(posequal[[i]]))[table(posequal[[i]]) != min(table(posequal[[i]]))] 

    ifelse(length(aux9[[i]]) == 0,
           aux9[[i]]  <- names(table(posequal[[i]]))[-1],
           aux9[[i]])   
   }

   # Escolhemos um deles para gerar um valor aleatório, pois o outro deverá ter um valor no qual a soma deva ser 1. Por comodidade escolhemos o primeiro no vetor.
   aux10  <- lapply(aux9, 
                    function(x) ifelse(length(x)==1, 
                                       x, 
                                       x[-length(x)]))

   # Este objeto retorna a frequência absoluta do equality constrain.
   aux11 <- list()

   for(i in namesvaluesfr) {

    aux11[[i]] <- aux7[[i]][aux10[[i]]]

   }

   # Agora, o último equality será determinado de modo que a soma nas linhas deva ser 1. Este valor só é necessário se aux9 > 1.
   if(all(unlist(lapply(aux9, 
                        function(x) any(length(x) > 1))) == TRUE)){

    # Este objeto retornb o(s) valor(es) dele(s).
    values2aux <- list(list())

    for(i in namesvaluesfr){

     for(j in 1:length(aux10[[i]])){

      values2aux[[i]][[j]] <- rep(runif(1)/K,
                                  aux11[[i]][j])

     }

     names(values2aux[[i]]) <- aux10[[i]]  

    }

    values3aux <- list()

    for(i in namesvaluesfr){

     values3aux[[i]] <- rep((1-sum(c(valuesfraux[[i]],
                                     values2aux[[i]][[1]])))/as.numeric(aux7[[i]][aux9[[i]][length(aux9[[i]])]]), 
                            as.numeric(aux7[[i]][aux9[[i]][length(aux9[[i]])]]))
    }

    values3 <- lapply(values3aux, 
                      function(x) unique(x))

    for(i in namesvaluesfr) {

     names(values3[[i]]) <- aux9[[i]][length(aux9[[i]])]

    }

    values2 <- lapply(values2aux, 
                      function(x) unlist(lapply(x, 
                                                function(x) unique(x))))

    values2 <- values2[!sapply(values2,
                               is.null)]

    valuesfr <- lapply(valuesfraux, 
                       function(x) unlist(x))

    valuesfr <- valuesfr[!sapply(valuesfr,
                                 is.null)] 

    valuesf <- mapply(function(x,y,z) c(x,y,z),
                      valuesfr,
                      values2,
                      values3,
                      SIMPLIFY=F)

    valuesf <- lapply(valuesf, 
                      function(x) x[order(names(x))])

    # Iremos bgora substituir os equalitys pelos valores encontrados.
    equali <- lapply(namesvaluesfr, 
                     function(x) sort(unique(posequal[[x]])))

    names(equali) <- namesvaluesfr

    for(i in namesvaluesfr){

     for(j in 1:length(equali[[i]])){

      a[al[[as.character(equali[[i]][j])]]] <- valuesf[[i]][j]

     }
    }
   } else {

    valuesfr <- lapply(valuesfraux, 
                       function(x) unlist(x))

    valuesfr <- valuesfr[!sapply(valuesfr,
                                 is.null)] 

    aux22 <- all(is.na(unique(unlist(aux2))) == TRUE) # Não tem nenhum fixed constrain?

    values2 <- list()

    for(i in 1:length(names(aux9))){

     values2[[i]] <- (1-sum(valuesfr[[i]],
                            ifelse(aux22,
                                   0,
                                   sum(aux2[[namesvaluesfr]]))))/aux7[[names(aux9)[[i]]]][unlist(aux9)[i]] 
    }

    values2 <- unique(unlist(values2))

    names(values2) <- unlist(aux10)

    valuesfr <- unique(unlist(valuesfr))

    names(valuesfr) <- sapply(aux6,function(x)x)

    valuesf <- c(valuesfr, 
                 values2)

    valuesf <- valuesf[order(names(valuesf))]

    # Iremos agora suastituir os equalitys pelos valores encontrados.
    equali <- unlist(lapply(namesvaluesfr, 
                            function(x) sort(unique(posequal[[x]]))))

    for(j in 1:length(equali)){

     a[al[[as.character(equali[j])]]] <- valuesf[j]

    }
   }  

   # DEVEMOS CHECAR PRIMEIRO SE ALGUMA LINHA COMPLETA FICOU COM EQUALITY SEM VALOR.

   aux14 <- a[,as.numeric(colnames(aux1)[colnames(aux1) != namesvaluesfr])]

   aux15 <- aux14[aux14 > 1]

   aux16 <- unique(aux15)

   namesaux15 <- names(table(aux15))

   if(length(aux16) > 0) {

    if(length(aux16) > 1){

     valuesau <- list()

     for(i in as.numeric(namesaux15)[-length(aux16)]){

      valuesau[[i]] <- rep(runif(1)/K,
                           table(aux15)[names(table(aux15)[as.character(i)])])

     }
     valuesau <- unlist(valuesau)

     valuesauf <- (1 - sum(aux14[aux14<1],
                           valuesau))/table(aux15)[names(table(aux15)[length(aux16)])] 
    } else {

     valuesauf <- (1 - sum(aux14[aux14<1]))/table(aux15)[names(table(aux15)[length(aux16)])]

     a[a == as.numeric(namesaux15)] <- valuesauf

    }
   }
  }
  # DEVEMOS CHECAR SE FICOU ALGUM EQUALITY SEM VALOR. 
  if(any(a > 1,
         na.rm = TRUE)){
   namesauxx <- sort(unique(a[a > 1][!is.na(a[a > 1])]))

   alauxx <- list()
   for(i in namesauxx){
    alauxx[[i]] <- which(a==i, 
                         arr.ind=TRUE)
    alauxx[[i]] <- alauxx[[i]][order(alauxx[[i]][,1],
                                     alauxx[[i]][,2]),]
    a[alauxx[[i]]] <- runif(1)/K
   }   
  } 

  awe <- which(is.na(a), 
               arr.ind=TRUE)

  awe <- awe[order(awe[,1],
                   awe[,2]), ] #positions of NA positions

  a[awe] <- runif(nrow(awe))/K #filling up NA positions with random values    
  y <- 0

  for(i in 1:nrow(a)) y[i] <- sum(awe[,1]==i)

  for(i in 1:nrow(a)){ 

   for(j in 0:K) if(y[i]==j){

    sba <- sum(a[i,][is.na(a1[i,])])

    sba.1 <- 1 - sum(a[i,][!is.na(a1[i,])])

    a[i,][is.na(a1[i,])] <- a[i,][is.na(a1[i,])]/sba*sba.1 

   }
  }
  #adjusted random values to satisfy the constraints sum = 1 
 } else { #only fixed parameters
  #creating random generated values for alpha(k|i)

  awf <- which(is.na(a), 
               arr.ind=TRUE)

  awf <- awf[order(awf[,1],
                   awf[,2]), ] #indices of cA = NA

  a[awf] <- runif(nrow(awf)) #random numbes at places with NA

  y <- 0

  for(i in 1:nrow(a)) y[i] <- sum(awf[,1]==i)

  for(i in 1:nrow(a)) y[i] <- sum(awf[,1]==i)

  for(i in 1:nrow(a)){ 

   for(j in 0:K){ 

    if(y[i]==j){

     sba <- sum(a[i,][is.na(a1[i,])])

     sba.1 <- 1-sum(a[i,][!is.na(a1[i,])])

     a[i,][is.na(a1[i,])] <- a[i,][is.na(a1[i,])]/sba*sba.1 

    }
   }
  }  
 }
 invisible(a)
}
