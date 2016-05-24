reliefcat <-
function (data, nosample, threshold, vnom,repet) 
{#require(StatMatch)
    p = dim(data)[2]
    class=as.numeric(factor(data[,p]))
    data <- cbind(data[,-p],class)
    f = p - 1
    acum <- rep(0, f)
    ngroups = length(unique(data[, p]))
    pesototal = rep(0, f)
    inst <- length(data[, 1])
    priors <- table(data[, p])/inst
    for (rep in 1:repet) {
        pesos <- rep(0, f)
        pesosmiss <- matrix(0, ngroups-1, f)
        indices <- sample(inst, nosample, replace = FALSE)
        muestra <- data[indices, ]
        for (i in 1:nosample) {
            datatemp <- data[-indices[i], ]
            data1=split.data.frame(datatemp[,1:f],datatemp[,p])
            indg <- muestra[i,p]
#Finding the gower distance to the nearhit
near=acugow(muestra[i,-p],data1[[indg]],vnom)
#print(class(nearhit)
ind1=order(near$dist)[1]
#print(ind1)
#Finding the nearhit weigths
pesonear=near$matdist[ind1,]
#Finding the nearmisses and their weights
for(j in 1:(ngroups-1)) {
            for (kk in (1:ngroups)[-indg]) {
                  #nearmiss <- near2(muestra[i,-p], data1[[kk]], vnom)
                    b=acugow(muestra[i,-p],data1[[kk]],vnom)
                    ind2=order(b$dist)[1]
 #                   print(ind2)
                    nearmiss=data1[[kk]][ind2,]
                         pesotempo=b$matdist[ind2,]*priors[kk]/(1-priors[indg]) 
 #               print(pesotempo)
                  }
           pesosmiss[j,]=pesotempo
       }
pesos=pesos+colSums(pesosmiss)/nosample-pesonear/nosample
}
        o1 <- order(-pesos)
        o2 <- pesos[o1]
        o3 <- o1[o2 > threshold]
        pesototal = pesototal + pesos
        acum[o3] = acum[o3] + 1
#cat("pesototal")
#print(pesototal)
    }
    pesotota = as.matrix(pesototal)
#print(pesotota)
    of1 <- order(-pesotota)
    of2 <- pesotota[of1]/repet
    acum = as.matrix(acum)
    tabla = cbind(1:f, acum, pesotota/repet)
    colnames(tabla) = c("feature", "frequency", "weight")
    tabla = tabla[order(-tabla[, 3]), ]
    cat("Features appearing in at least half of repetitions ordered by their average relevance weight: \n")
    print(tabla[tabla[, 2] >= repet/2, ])
    plot(of2, ylab = "weights")
    text(1:f, of2, tabla[, 1], pos = 4)
    relevant1 = which(acum >= repet/2)
    relevant2 = which(pesotota/repet > threshold)
    relevant = c(relevant1, relevant2)
    relevant = relevant[duplicated(relevant)]
    cat("selected features", "\n")
    relevant = tabla[1:length(relevant), 1]
    return(relevant)
}
