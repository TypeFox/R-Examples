reliefcont <-
function (data, nosample, threshold,repet) 
{
 #Data containing  continuous features must be normalized first using mmnorm
    p = dim(data)[2]
    class=as.numeric(factor(data[,p]))
    data <- cbind(data[,-p],class)
    data=as.matrix(data)
    f = p - 1
    acum <- rep(0, f)
    features <- seq(f)
    ngroups = length(unique(data[, p]))
    pesototal = rep(0, f)
    inst <- length(data[, 1])
    priors <- table(data[, p])/inst
    dh <- rep(0, f)
    for (j in 1:f) {
        dh[j] <- diff(range(data[, j]))
    }
    for (rep in 1:repet) {
        nearhit <- matrix(0, nosample, f)
        pesos <- rep(0, f)
        tempo <- matrix(0, ngroups, f)
        indices <- sample(inst, nosample, replace = FALSE)
        muestra <- data[indices, ]
        #class=as.numeric(factor(data[,p]))
        for (i in 1:nosample) {
            datatemp <- data[-indices[i], ]
            #data1 = split.data.frame(datatemp[, 1:f], class[-indices[i]])
            data1=split.data.frame(datatemp[,1:f],datatemp[,p])
            indg <- muestra[i,p]
#           print(indg)
            nearhit[i, ] <- near3(muestra[i,-p], data1[[indg]])
#           print(nearhit[i,])
            for (kk in 1:ngroups) {
                if (kk != indg) {
                  nearmiss <- near3(muestra[i,-p], data1[[kk]])
                  tempo[kk, ] <- (muestra[i,-p] - as.vector(nearmiss))
                }
                for (ii in 1:f) {
                  tempo[kk, ii] <- (1/nosample) * abs(tempo[kk, ii])/dh[ii]
                }
            }
            pesomiss <- rep(0, f)
            for (jj in 1:f) {
                for (kk in 1:ngroups) {
                  if (kk != indg) {
                    pesomiss[jj] <- pesomiss[jj] + priors[kk] * 
                      tempo[kk, jj]
                  }
                }
                pesomiss[jj] <- pesomiss[jj]/(1 - priors[indg])
            }
            for (j in 1:f) {
                diff <- -(1/nosample) * abs(muestra[i,j] - nearhit[i, j])/dh[j] + pesomiss[j]
                pesos[j] <- pesos[j] + diff
            }
        }
        o1 <- order(-pesos)
        o2 <- pesos[o1]
        o3 <- o1[o2 > threshold]
        pesototal = pesototal + pesos
        acum[o3] = acum[o3] + 1
    }
    pesotota = as.matrix(pesototal)
    of1 <- order(-pesotota)
    of2 <- pesotota[of1]/repet
    acum = as.matrix(acum)
    tabla = cbind(1:f, acum, pesotota/repet)
    colnames(tabla) = c("feature", "frequency", "weight")
    tabla = tabla[order(-tabla[, 3]), ]
    cat("Features appearing in at least half of repetitions ordered by their average relevance weight: \n")
    print(tabla[tabla[, 2] >= repet/2, ])
    plot(of2, ylab = "weights")
    text(1:f, of2, tabla[, 1], cex = 0.7, pos = 4)
    relevant1 = which(acum >= repet/2)
    relevant2 = which(pesotota/repet > threshold)
    relevant = c(relevant1, relevant2)
    relevant = relevant[duplicated(relevant)]
    cat("selected features", "\n")
    relevant = tabla[1:length(relevant), 1]
    return(relevant)
}
