computeOrdering <- function(data, varname, test) {
  if (test == "FCH") {
  ord <- fold.change(data,varname)
  for (i in 1:ncol(ord))  {
    ord[,i]<- rownames(ord)[order(ord[,i], decreasing = TRUE)]
  }
  }
  if (test == "T"){
  ord<-meta.test(data, varname)$test
  for (i in 1:ncol(ord))  {
    ord[,i]<- rownames(ord)[order(ord[,i], decreasing = TRUE)]	
  }
  }
  if (test == "") {}
  rownames(ord) <- NULL
  return(ord)
}

flip <- function(order) {
order<-order[nrow(order):1,]
return(order)
}

commonGenes<-function(ord,n)
{
r<-apply(ord[,-1], 2,function(x) match(ord[,1],x))
r<-apply(r,1,max)
or<-pmax(r[1:n],1:n)
tmp <- table(or)
x <- integer(n)
x[as.integer(names(tmp))] <- tmp
x <- cumsum(x)
return(x[1:n])
}





prelimScore<-function(ordering, alpha, min.weight = 1e-05, two.sided = TRUE) {
n<- -log(min.weight)/alpha
comm.dir  <- commonGenes(ordering,n)
if (two.sided) {
comm.flip <- commonGenes(flip(ordering),n)
cg <- comm.dir + comm.flip
} else {cg <- comm.dir}

w <- exp(-alpha*c(1:length(cg)))
pS<- sum(w * cg)
#pS<- weighted.mean(cg, w)
return(pS)
}





###
preparePermutations<-function(id, B, sample.ratio = 0.8)
{
    n1 <- floor(sum(id) * sample.ratio)
    n2 <- floor(sum(1 - id) * sample.ratio)
    per<- function ( v, m)
    {
    x <- sort(unique(v))
    M <- matrix(NA, m, length(v))
    #M[1, ] <- v #prvy riadok povodne - potrebujem??
    for (i in 2:m) {M[i, ] <- sample(v)}
    return(M[-1,])
    }
        yperm <- per(id, B + 1)
        ysubs <- matrix(nrow = B, ncol = (n1 + n2))
        for (i in 1:B) {
            x <- sample(1:(sum(id)), n1)
            y <- sample((sum(id) + 1):(length(id)), n2)
            ysubs[i, ] <- c(x, y)
        }

    return(list(yperm = yperm, ysubs = ysubs))
}

RandomScore<-function(data, varname, B, alpha, test, which=c("random", "empirical", "subsample"), two.sided = TRUE){
N=length(clinical(data))
n<-nrow(GEDM(data)[[1]])
random = NULL
empirical.ci = NULL
subsample = NULL
classlab<-list()
for (i in 1:N) classlab[[i]]<-as.numeric(clinical(data)[[i]][,varname])-1
classlab<-sapply(classlab, function(x) 1-x)
p<-sapply(classlab, preparePermutations, B)


if ("random" %in% which) {
pp<-p[1,]
Score<- function(j, data, varname, pp, test, two.sided) {
    prData<-prepareData(j, data, varname, pp, 1)
    ordering <- computeOrdering(prData, varname, test)
    comm <- commonGenes( ordering,n)
    if (two.sided) {
      comm2<- commonGenes( flip(ordering),n)
      cg <- comm + comm2	} else cg <- comm
    SC<-sapply(as.list(alpha), function(x,n,cg) {sum(exp(-x*c(1:n))*cg)}, n,cg)
  return(SC)
  }
random<-sapply(1:B,Score, data, varname, pp, test, two.sided)

}

if ("empirical" %in% which) {
pp<-p[1,]
Empirical<-function(j, data, varname, p, test) {
prData<-prepareData(j, data, varname, p, 1)
ordering <- computeOrdering(prData, varname, test)
cg<-commonGenes(ordering, n)
cg2<-commonGenes(flip(ordering),n)
res<-list(top=cg, bottom=cg2)
return(res)
}
res<-sapply(1:B,Empirical, data, varname, pp, test)
top<- t(matrix(unlist(res["top",]), nrow = B, byrow = TRUE))
bottom<- t(matrix(unlist(res["bottom", ]), nrow = B, byrow = TRUE))
top.ci<- t(apply(top, 1, quantile, probs = c(0.025, 0.5, 0.975)))
bottom.ci<- t(apply(bottom, 1, quantile, probs = c(0.025, 0.5, 0.975)))
empirical.ci <- list(top = top.ci, bottom = bottom.ci)
}
if ("subsample" %in% which) {
pp<-p[2,]
Subsample<- function(j, data, varname, pp, test, two.sided) {
    preData<-prepareData(j, data, varname, pp, 2)
    ordering <- computeOrdering(preData, varname, test)
    comm <- commonGenes( ordering,n)
    if (two.sided) {
      comm2<- commonGenes( flip(ordering),n)
      cg <- comm + comm2	} else cg <- comm
    subSC<-sapply(as.list(alpha), function(x,n,cg) {sum(exp(-x*c(1:n))*cg)}, n,cg)
  return(subSC)
  }
subsample<-sapply(1:B,Subsample, data, varname, pp, test, two.sided)
}


return(res=list(random = random, empirical.ci =  empirical.ci, subsample = subsample))
}



prepareData<-function(j, data, varname , p, type) {
if (type == 1) {
cl<-clinical(data)
for (i in 1:length(cl)) {
cl[[i]][,varname]<-p[[i]][j,]+1
}
perData<- new ("MetaArray", GEDM = GEDM(data), clinical = cl, datanames = datanames(data))
}
if (type == 2) {
dumExpr<-list()
dumClin<-list()
for (i in 1:length(GEDM(data))) {
dumExpr[[i]]<-GEDM(data)[[i]][,p[[i]][j,]]
dum<- as.data.frame(clinical(data)[[i]][p[[i]][j,],])
if (dim(dum)[2] == 1) {colnames(dum) <- varname
    rownames(dum) <- colnames(dumExpr[[i]]) }
dumClin[[i]]<- dum
}
perData<- new ("MetaArray", GEDM = dumExpr, clinical = dumClin, datanames = datanames(data))
}
return(perData)
}



##
# Vypocet alfa tak, ze pri exp vahach na hodnote <n> klesne vaha pod <min.weight>
# Hodnoti sa maximalne 2500 genov
computeAlpha <- function(n = NULL, min.weight = 1e-05, ngenes) {
if (is.null(n)) {
  n <- c(100, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 2500)
  alpha <- -log(min.weight)/n
  } else {
    alpha <- -log(min.weight)/n
    }
    select <- n <= ngenes
    alpha <- alpha[select]
    return(alpha)
}


# Vypocet pauc podla OrderedList

# Vyber alfa podla pAUC
selectAlpha <- function (alpha, subsample, random){
pAUC <- numeric(length(alpha))
B <- dim(subsample)[2]
y <- c(rep(0, B), rep(1, B))
pauc <- function(x, A, B) {
        n <- length(A)
        o <- order(A, decreasing = TRUE)
        r <- c(0, cumsum(B[o[-n]]))
        t <- (0:(n - 1)) - r
        r <- r/sum(B)
        t <- t/sum(1 - B)
        roc <- numeric(length(x))
        for (i in 1:length(x)) {
            z <- which(t <= x[i])
            z <- z[length(z)]
            roc[i] <- r[z]
        }
        return(roc)
    }  
for (i in 1:length(alpha)) {
        X <- c(random[i,], subsample[i,])
        pAUC[i] <- integrate(pauc, 0, 0.1, A = X, B = y, stop.on.error = FALSE)$value
    }
    x=list(alpha = alpha[which.max(pAUC)], pAUC=pAUC)
    return(x)
}

 

# Vyznamnost skore
sigScore <- function (ranking, alpha, B, min.weight = 1e-05, two.sided = TRUE) {
s <- prelimScore(ranking, alpha, min.weight, two.sided ) 
# random rankings and their score
rs <-  numeric(B)
rrank <- ranking
for (i in 1:B) {
  
  rrank<-apply(rrank, 2, function(x) sample(x))
  rs[i]<- prelimScore(rrank, alpha, min.weight, two.sided)
  }
# empirical probability
return (sum(rs>=s)/B)
}

# Geny
selectGenes<-function(ordering, alpha, percent, min.weight = 1e-05, two.sided = TRUE) {
n<- -log(min.weight)/alpha
comm.dir  <- commonGenes(ordering,n)
if (two.sided) {
comm.flip <- commonGenes(flip(ordering),n)
cg <- comm.dir + comm.flip} else cg<- comm.dir
w <- exp(-alpha*c(1:length(cg)))
pS <- w * cg
y <- min(which(cumsum(pS) >= percent * sum(pS) - (1e-05))) 
y1 <- colIntersect(ordering[1:y,])
y2 <- colIntersect(flip(ordering)[1:y,])
genes<- sort(c(y1,y2))
return(genes)
}


colIntersect <- function(x)
{
  N<-ncol(x)
  dum<-intersect(x[,1], x[,2])
  if (N >= 3) {
  for (i in 3:N) {
    dum<-intersect(dum,x[,i])
  }            
  }
  return(dum)
}


      
    


# Wraper function
performSOGL <- function(data, varname, test, B, which=c("score", "empirical"), min.weight = 1e-05, two.sided = TRUE, percent = 0.95){
cat("Processing data...")
if (!all(sapply(1:(length(GEDM(data))-1), function(x) all(rownames(GEDM(data)[[x]])==rownames(GEDM(data)[[x+1]]))))) stop("The gene expression data matrices have not equal rownames")
all.genes<-rownames(GEDM(data)[[1]])
ordering<- computeOrdering(data, varname, test)
A<-computeAlpha(ngenes=nrow(GEDM(data)[[1]]))
cat("Tuning alpha..")
if  ("empirical" %in% which) sampl <- c("random", "empirical", "subsample") else sampl <- c("random", "subsample")
sampling<-RandomScore(data, varname, B, A, test, which=sampl)
a<-selectAlpha(A, sampling$subsample, sampling$random)
score<-prelimScore(ordering, a$alpha)
n<- -log(min.weight)/a$alpha
comm.dir  <- commonGenes(ordering,n)
if (two.sided) {
comm.flip <- commonGenes(flip(ordering),n)
cg <- list(top = comm.dir, bottom = comm.flip) } else cg <- list(top = comm.dir, bottom = NULL)
cat("Significance and genes...")
sig<-sigScore(ordering, a$alpha, B)
genes<-selectGenes(ordering, a$alpha, percent)

alph.pos<-match(a$alpha, A)
res <- list(ordering = ordering, alpha.selected = a$alpha, alpha.considered = A, pAUC=a$pAUC, 
	random = sampling$random[alph.pos,], subsample = sampling$subsample[alph.pos,], emp.ci = sampling$empirical.ci, common.genes = cg, score = score, significance = sig, genes = genes, all.genes=all.genes)
class(res)<-"SOGLresult"
return(res)
}




