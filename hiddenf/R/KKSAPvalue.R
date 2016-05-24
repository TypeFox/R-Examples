KKSAPvalue <-
function(hfobj)
{
# a <- nlevels(hfobj$tall$rows)
r <- nlevels(hfobj$tall$rows)
if(r<4){print("KKSA needs at least 4 levels of row factor");return(list(pvalue=NA))}
rows <- hfobj$tall$rows
cols <- hfobj$tall$cols
y <- hfobj$tall$y
rkconfig.mtx <- rkconfig.fcn(rows)
rows <- as.factor(rows)
# a <- length(table(rows))
# if(a<4)stop("KKSA only applicable if a > 3")
cols <- as.factor(cols)
# b <- length(table(cols))
c <- length(table(cols))
# cc <- 2^(a-1)-1-a # no singletons
cc <- 2^(r-1)-1-r # no singletons
fvalues <- rep(NA,cc)
pvalues <- rep(NA,cc)
for(i in 1:cc)
{
i1 <- (rkconfig.mtx[,i]==1)
i2 <- !i1
y1.tmpout <- lm(y[i1] ~ cols[i1] + rows[i1])
ms1 <- anova(y1.tmpout)[3,3]
df1 <- (sum(i1)/c-1)*(c-1)
y2.tmpout <- lm(y[i2] ~ cols[i2] + rows[i2])
ms2 <- anova(y2.tmpout)[3,3]
df2 <- (sum(i2)/c-1)*(c-1)
fstat <- ms1/ms2
finv <- ms2/ms1
if(fstat>finv){
fmax <- fstat
# dfn <- df0
# dfd <- df1
}
else{
fmax <- finv
# dfn <- df1
# dfd <- df0
}
# print(c(i,fstat,finv,fmax))
fvalues[i] <- fmax
#pvalues[i] <- 1-pf(fmax,dfn,dfd)
#pvalues[i] <- min(1,2*(1-pf(fmax,dfn,dfd))) HERE
#pvalues[i] <- 1-pf(fmax,df0,df1)+pf(1/fmax,df0,df1)
pvalues[i] <- 1-pf(fmax,df1,df2)+pf(1/fmax,df1,df2)
}
pvalue <- min(pvalues)*cc
pvalue <- min(1,pvalue)
config <- which.min(pvalues)

config.vector <- rkconfig.mtx[,config]
grp.vector <- 1+config.vector[1+c*(0:(r-1))]
tall <- list(y=y,rows=rows,cols=cols)
# kksa.out <- list(pvalues=pvalues,pvalue=pvalue,tall=tall,fvalues=fvalues)
#kksa.out <- list(pvalues=pvalues,pvalue=pvalue,config=config,config.vector=config.vector,tall=tall)
#KKSA.out <- list(fmax <- fvalues[config],pvalue=pvalue,grp.vector=grp.vector,tall=tall,NumDf=dfn,DenomDf=dfd)
#KKSA.out <- list(fmax = fvalues[config],pvalue=pvalue,grp.vector=grp.vector,tall=tall,NumDf=df1,DenomDf=df2)
KKSA.out <- list(fmax = fvalues[config],pvalue=pvalue,grp.vector=grp.vector,NumDf=df1,DenomDf=df2)
return(KKSA.out)
}
