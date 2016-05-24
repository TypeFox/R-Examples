gao_cs <-
function (formula, data, alpha=0.05, silent = FALSE){

#-------Specify the one way factorial model -----------------------------------#
    dat <- model.frame(formula, data)
    if (ncol(dat) != 2) {
        stop("Specify one response and only one class variable in the formula")
    }
    if (is.numeric(dat[, 1]) == FALSE) {
        stop("Response variable must be numeric")
    }
    res<- dat[, 1]
N <- length(res)
    group <- as.factor(dat[, 2])
    fl <- levels(group)
    a <- nlevels(group)
    nc <- a * (a - 1)/2


#------Compute the joint ranks for the first step -----------------------------#

    response <-1/N*(rank(res)-1/2)
    
#-------Compute the rank means and sort them in the next step------------------#

    samples <- split(response, group)
    n <- sapply(samples, "length")
    mm <- sapply(samples, "mean")
    vv <- sapply(samples, "var")

#------------sort the data according to the order of means---------------------#

    order.h1 <- data.frame(Sample = fl, Size = n, Effect = mm, Variance = vv)
    ordered <- cbind(Order = 1:a, order.h1[order(order.h1$Effect, decreasing = FALSE),] )
    rownames(ordered) <- 1:a
    ordered.names <- 1:a
    
    grp.neusort<-list()
    res.neusort<-list()
    for (jj in 1:a){
    sample.name <- ordered$Sample[jj]
    res.neusort[[jj]] <- res[group==sample.name]
    grp.neusort[[jj]] <-rep(jj,length(res.neusort[[jj]]))}
    resp.sort<-unlist(res.neusort)
    grp.sort <-unlist(grp.neusort)

 #-----------Prepare an overview table for the output---------------------------#

    i <- 1:(a - 1)
    h1 <- list()
    for (s in 1:(a - 1)) {
        h1[[s]] <- i[1:s]
    }
    vi <- unlist(h1)
    j <- a:2
    h2 <- list()
    for (s in 1:(a - 1)) {
        h2[[s]] <- j[s:1]
    }
    vj <- unlist(h2)
    h3 <- list()
    h4 <- list()
    for (s in 1:(a - 1)) {
        h3[[s]] <- rep(j[s], s)
        h4[[s]] <- rep(i[s], s)
    }
    Nmean <- unlist(h3)
    Step <- unlist(h4)


    
#----------Compute the rank means according to the order defined previous------#
    rel.diffs <- c()
variance.diffs <- c()
df.diffs <- c()

    for (l in 1:nc){
    ii<-vi[l]
    jj<-vj[l]
    i<-ordered.names[ii]
    j<-ordered.names[jj]


    h5<-i:j

    h6<-list()
    h66<-list()
    for (s in 1:length(h5)) {

    h6[[s]] <- resp.sort[grp.sort==h5[s]]
h66[[s]] <-rep(h5[s],length(h6[[s]]))}

uh6 <- unlist(h6)


#--------Rerank the data for the computation of the statistics-----------------#
    reranks <- 1/length(uh6)*(rank(uh6)-1/2)
    regroup <- unlist(h66)

    resamples <- split(reranks,regroup)

    h7 <- sapply(resamples,"mean")


    h77 <- sapply(resamples,"var")
    h777 <- sapply(resamples,"length")
lh5 <- length(h5)
  rel.diffs[l]<- abs(h7[lh5] - h7[1])
    v.diff1 <- h77[lh5]/h777[lh5]
v.diff2 <- h77[1]/h777[1]
 variance.diffs[l] <-v.diff1 + v.diff2
df.diffs[l] <-(v.diff1+v.diff2)^2 /(v.diff1^2*1/(h777[lh5]-1)+ v.diff2^2*1/ (h777[1]-1))
}
Statistics <- round(sqrt(2)*rel.diffs / sqrt(variance.diffs),4)

#--------------------Compute the pvalues and Quantiles-------------------------#

    pvalues <- round(ptukey(Statistics, Nmean, df.diffs, lower.tail = FALSE),4)
    alpha.level <- round(1 - (1 - alpha)^(Nmean/a),4)
    level1 <- (Nmean == a)
    level2 <- (Nmean == a - 1)
    level3 <- level1 + level2
    alpha.level[level3 == 1] <- alpha
    quantiles <- round(qtukey(1 - alpha.level, Nmean, df.diffs),4)

#-------------Attention: quantiles are not monotone!!!-------------------------#

    for (h in 1:(nc - 1)) {
        if (quantiles[h + 1] >= quantiles[h]) {
            quantiles[h + 1] <- quantiles[h]
        }
    }
    Rejected1 <- (Statistics > quantiles)

#-----------Logical section for the subsets S being not rejected before--------#
    
    for (s in 1:nc) {
        if (Rejected1[s] == FALSE) {
            Under1 <- (vj[s] >= vj)
            Under2 <- (vi[s] <= vi)
            Under3 <- Under1 * Under2
            Under4 <- which(Under3 == 1)
            Rejected1[Under4] <- FALSE
        }
    }
    Out1 <- (Statistics>quantiles)
    Out2 <- (Rejected1 == FALSE)
    Out3 <- Out1 * Out2
    Out4 <- (Out3 == 1)
    
    
    pvalues[Out4] <- paste(">", alpha.level[Out4])
    quantiles[Out4] <- paste(">", Statistics[Out4])

#------Compute Bonferroni and Holm adjusted p-values for the single tests------#

Tests <-sapply(1:nc, function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        ( ordered$Effect[j] - ordered$Effect[i])/sqrt(ordered$Variance[i]*1/ordered$Size[i] +
            ordered$Variance[j]*1/ordered$Size[j])
    })
Effects <- sapply(1:nc, function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        (ordered$Effect[j] - ordered$Effect[i])
    })
df.single <- sapply(1:nc, function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        v.single1 <- ordered$Variance[i] /ordered$Size[i]
        v.single2 <- ordered$Variance[j] /ordered$Size[j]
        (v.single1 + v.single2)^2 / (v.single1^2*1/(ordered$Size[i]-1)+ v.single2^2*1/ (ordered$Size[j]-1))
    })


p.tapp1 <- pt(Tests, df = df.single)
p.tapp <- apply(cbind(2*p.tapp1, 2-2*p.tapp1),1,"min")
p.bonf <- p.adjust(p.tapp, method = "bonferroni")
p.holm <- p.adjust(p.tapp, method = "holm")



#------------Prepare a nice overview for the output----------------------------#

names.ordered <- sapply(1:nc, function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        paste(ordered$Sample[j], "-", ordered$Sample[i], sep = "")
    })

#---------------------------Table output---------------------------------------#

Single.Analysis <- data.frame(Comp = names.ordered,
Effect = round(Effects,4),
Statistic = round(Tests,4),
DF = round( df.single,4),
P.RAW = round(p.tapp,4),
p.BONF = round( p.bonf,4),
p.HOLM = round(p.holm,4))

    MCP <- list(Info=ordered,
    Single.Analysis = Single.Analysis,
CS.Analysis=data.frame(
Comp = names.ordered,
Effect = round(rel.diffs,4),
    Statistic = Statistics,
DF = round(df.diffs,4),
Quantiles =quantiles,
Adj.P = pvalues,
    Alpha = alpha.level,
Rejected = Rejected1, Layer = Step))

#----------------------Information about the procedure-------------------------#
    if (!silent) {
    cat("\n","#----Gao et al's (2008) modification of Campbell and Skillings (1985) (CS) stepwise multiple comparison procedure  \n",
        "#---- This function uses joint ranks of the data. Attention: In the CS algorithm, the samples are jointly reranked! \n",
        "#----Reference: Gao, X. et al. (2008). Nonparametric Multiple Comparison Procedures for Unbalanced One-Way Factorial Designs. JSPI 138, 2574 - 2591. \n")}

return(MCP)
}
