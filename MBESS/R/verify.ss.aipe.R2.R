verify.ss.aipe.R2 <- function(Population.R2=NULL, conf.level=.95, width=NULL, Random.Predictors=TRUE, 
which.width="Full", p=NULL, n=NULL, degree.of.certainty=NULL, g=500, G=10000, print.iter=FALSE, ...)
{
# Based on method called for Sample size.

print("An internal Monte Carlo simulation has begun and the results may take some time.")
print("Please be patient, as thousands and thousands of replications are being performed")
print("to ensure the returned value of sample size is the exact value to satisfy the specified goals.")

Summary.of.SS <- ss.aipe.R2.sensitivity(Selected.N=n,
True.R2=Population.R2, Estimated.R2=NULL, w=width, p=p, 
Random.Predictors=Random.Predictors, degree.of.certainty=degree.of.certainty, 
conf.level=conf.level, Generate.Random.Predictors=Random.Predictors, 
rho.yx=.3, rho.xx=.3, G=g, print.iter=print.iter)$Summary

i <- 1

if(is.null(degree.of.certainty))
{

E.Width.Info <- as.matrix(cbind(E.Width.CI=Summary.of.SS$mean.CI.width.R2, N=n))
w.Bigger <- E.Width.Info[1,1] > width
w.Smaller <- E.Width.Info[1,1] < width

while(i != 0)
{
i <- i + 1

if(w.Bigger==TRUE) n <- n + 1
if(w.Smaller==TRUE) n <- n - 1

Summary.of.SS <- ss.aipe.R2.sensitivity(Selected.N=n, True.R2=Population.R2,
Estimated.R2=NULL, w=width, p=p, Random.Predictors=Random.Predictors,
degree.of.certainty=NULL, conf.level=conf.level,
Generate.Random.Predictors=Random.Predictors, 
rho.yx=.3, rho.xx=.3, G=g, print.iter=print.iter)$Summary

E.Width.Info <- rbind(E.Width.Info, c(E.Width.CI=Summary.of.SS$mean.CI.width.R2, N=n))

if(w.Bigger==FALSE) w.Bigger <- E.Width.Info[i,1] > width
if(w.Smaller==FALSE) w.Smaller <- E.Width.Info[i,1] < width

if((w.Bigger==TRUE) & (w.Smaller==TRUE)) break()
}

# Now get exact sample size (using previous results).
n <- E.Width.Info[nrow(E.Width.Info),2]

# Based on method called for Sample size.
Summary.of.SS <- ss.aipe.R2.sensitivity(Selected.N=n,
True.R2=Population.R2, Estimated.R2=NULL, w=width, p=p, 
Random.Predictors=Random.Predictors, degree.of.certainty=degree.of.certainty, 
conf.level=conf.level, Generate.Random.Predictors=Random.Predictors, 
rho.yx=.3, rho.xx=.3, G=G, print.iter=print.iter)$Summary

E.Width.Info <- as.matrix(cbind(E.Width.CI=Summary.of.SS$mean.CI.width.R2, N=n))

i <- 1

w.Bigger <- E.Width.Info[1,1] > width
w.Smaller <- E.Width.Info[1,1] < width

while(i != 0)
{
i <- i + 1

if(w.Bigger==TRUE) n <- n + 1
if(w.Smaller==TRUE) n <- n - 1

Summary.of.SS <- ss.aipe.R2.sensitivity(Selected.N=n, True.R2=Population.R2,
Estimated.R2=NULL, w=width, p=p, Random.Predictors=Random.Predictors,
degree.of.certainty=NULL, conf.level=conf.level,
Generate.Random.Predictors=Random.Predictors, 
rho.yx=.3, rho.xx=.3, G=G, print.iter=print.iter)$Summary

E.Width.Info <- rbind(E.Width.Info, c(E.Width.CI=Summary.of.SS$mean.CI.width.R2, N=n))

if(w.Bigger==FALSE) w.Bigger <- E.Width.Info[i,1] > width
if(w.Smaller==FALSE) w.Smaller <- E.Width.Info[i,1] < width

if((w.Bigger==TRUE) & (w.Smaller==TRUE)) break()
}

Contending <- (E.Width.Info[,1] <= width)
To.Use <- min(E.Width.Info[Contending,2])
Here <- E.Width.Info[which(E.Width.Info[,2]==To.Use),]

Result.Full <- list(Required.Sample.Size=as.numeric(Here[2]), Expected.Width=as.numeric(Here[1]))
Result.Full <- as.numeric(Here[2])

print("Information on the Monte Carlo simulation regarding the expected widths follows:")
print(E.Width.Info)

return(Result.Full)
}
############################################################################
############################################################################
# The following is used if there is a desired degree of certainty specified.
############################################################################
############################################################################

if(!is.null(degree.of.certainty))
{

E.Gamma.Info <- as.matrix(cbind(E.Gamma=Summary.of.SS$Pct.Less.w, N=n))
g.Bigger <- E.Gamma.Info[1,1] > degree.of.certainty
g.Smaller <- E.Gamma.Info[1,1] < degree.of.certainty


while(i != 0)
{
i <- i + 1

if(g.Bigger==TRUE) n <- n - 1
if(g.Smaller==TRUE) n <- n + 1

Summary.of.SS <- ss.aipe.R2.sensitivity(Selected.N=n, True.R2=Population.R2,
Estimated.R2=NULL, w=width, p=p, Random.Predictors=Random.Predictors,
degree.of.certainty=NULL, conf.level=conf.level,
Generate.Random.Predictors=Random.Predictors, 
rho.yx=.3, rho.xx=.3, G=g, print.iter=print.iter)$Summary

E.Gamma.Info <- rbind(E.Gamma.Info, c(E.Gamma=Summary.of.SS$Pct.Less.w, N=n))

if(g.Bigger==FALSE) g.Bigger <- E.Gamma.Info[i,1] > degree.of.certainty
if(g.Smaller==FALSE) g.Smaller <- E.Gamma.Info[i,1] < degree.of.certainty

if((g.Bigger==TRUE) & (g.Smaller==TRUE)) break()
}

# Now get exact sample size (using previous results).
n <- E.Gamma.Info[nrow(E.Gamma.Info),2]

# Based on method called for Sample size.
Summary.of.SS <- ss.aipe.R2.sensitivity(Selected.N=n,
True.R2=Population.R2, Estimated.R2=NULL, w=width, p=p, 
Random.Predictors=Random.Predictors, degree.of.certainty=degree.of.certainty, 
conf.level=conf.level, Generate.Random.Predictors=Random.Predictors, 
rho.yx=.3, rho.xx=.3, G=G, print.iter=print.iter)$Summary

E.Gamma.Info <- as.matrix(cbind(E.Gamma=Summary.of.SS$Pct.Less.w, N=n))

i <- 1

g.Bigger <- E.Gamma.Info[1,1] > degree.of.certainty
g.Smaller <- E.Gamma.Info[1,1] < degree.of.certainty

while(i != 0)
{
i <- i + 1

if(g.Bigger==TRUE) n <- n - 1
if(g.Smaller==TRUE) n <- n + 1

Summary.of.SS <- ss.aipe.R2.sensitivity(Selected.N=n, True.R2=Population.R2,
Estimated.R2=NULL, w=width, p=p, Random.Predictors=Random.Predictors,
degree.of.certainty=NULL, conf.level=conf.level,
Generate.Random.Predictors=Random.Predictors, 
rho.yx=.3, rho.xx=.3, G=G, print.iter=print.iter)$Summary

E.Gamma.Info <- rbind(E.Gamma.Info, c(E.Gamma=Summary.of.SS$Pct.Less.w, N=n))

if(g.Bigger==FALSE) g.Bigger <- E.Gamma.Info[i,1] > degree.of.certainty
if(g.Smaller==FALSE) g.Smaller <- E.Gamma.Info[i,1] < degree.of.certainty

if((g.Bigger==TRUE) & (g.Smaller==TRUE)) break()
}

Contending <- (E.Gamma.Info[,1] >= degree.of.certainty)
To.Use <- min(E.Gamma.Info[Contending,2])
Here <- E.Gamma.Info[which(E.Gamma.Info[,2]==To.Use),]

# Result.Full <- list(Required.Sample.Size=as.numeric(Here[2]), Empirical.Degree.of.Certainty=as.numeric(Here[1]))
print("Information on the Monte Carlo simulation regarding the empirical degree of certainty follows:")
print(E.Gamma.Info)
Result.Full <- as.numeric(Here[2])
return(Result.Full)
}
}
