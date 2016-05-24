`setModel` <-
function(n.components, ploidy.level, random.effect=FALSE,
                     seg.ratios=NULL, ploidy.name=NULL, equal.variances=TRUE,
                     type.parents=c("heterogeneous","homozygous"))
{

  ## mainly used just to generate bugs code without priors
  ## replace STEM on line 1, UNSPECIFIED on line 3 etc
  ## eg by using > gsub("STEM","test",bc) etc etc

  ## ploidy level etc - not really necessary to define it here but it
  ## might as well go here

  type <- match.arg(type.parents)
  E.segRatio <- expected.segRatio(ploidy.level, type.parents=type)

  if (random.effect & equal.variances)
    warning("Both 'random.effect' and 'equal.variances' set - silly?")
  
  ## set up segregation proportions if specified
  if (length(seg.ratios)>1){
    E.segRatio$ratio <- seg.ratios
    if (n.components>length(seg.ratios))
      stop("Number of components must be less that no. of segregation ratios")
    if (length(ploidy.name>1)){
      E.segRatio$ploidy.name <- ploidy.name
    }
  }
  
  E.segRatio <- expected.segRatio(ploidy.level, type.parents=type)

  ## set up bugs code
  
  bc <- c("## File:    STEM.bug",
          "## Generated using makeModel {polySegratioMM} at",
          paste("## Purpose: fit", n.components,
                "component mixture for", E.segRatio$ploidy.name,"with",
                type,"parents"),
          "##          Priors: see below")
  if (random.effect) {
    bc <- c(bc,"##          WITH a random effect to approximate measurement",
            "##              error on logit scale under binomial model")
  } else {
    bc <- c(bc,"##          WITHOUT a random effect")
  }

  if (equal.variances) {
    bc <- c(bc,"##          EQUAL VARIANCES for each component")
  } else {
    bc <- c(bc,"##          SEPARATE VARIANCES for each component")
  }
  ## variable names
  bc <- c(bc,"var",
          "    n[N],         # bionomial n for each marker",
          "    r[N],         # r = no. of ones for each marker",
          "    pbin[N],         # segregation ratio for each marker",
          "    p[N],            # logit seg ratio without measurement error",
          "    T[N],            # true groups (labelled 1,2,3, ...)",
          paste("    mu[",n.components,"],           # means of groups",
                sep=""),
          paste("    theta[",n.components-1,
                "],        # scaled positive shift between groups",sep=""))
  if (equal.variances) {
    bc <- c(bc,     
            "    tau,            # sampling precision",
            "    sigma,          # sampling standard deviation")
  } else {
    bc <- c(bc,     
            paste("    tau[",n.components,
                  "],          # sampling precision",sep=""),
            paste("    sigma[",n.components,
                  "],        # sampling standard deviation",sep=""))
  }
  
  bc <- c(bc,     
          paste("    alpha[",n.components,
                "],        # prior parameters for proportions",sep=""))
  if (!random.effect){
    bc <- c(bc,"    b[N],            # N(0,taub) random effect for each observation",
            "    taub,            # tau of random effect")
  }
  bc <- c(bc,paste("    P[",n.components,
                   "];            # proportion in each group",sep=""))
  ## model and transforms
  bc <- c(bc,"model {",
          "#   T[N] fixed for several values to aid convergence",
          "   for (i in 1:N){")
  if (equal.variances) {  
    bc <- c(bc,"       p[i]~ dnorm(mu[T[i]],tau);")
  } else {
    bc <- c(bc,"       p[i]~ dnorm(mu[T[i]],tau[T[i]]);")
  }
  if (random.effect){
    bc <- c(bc,"       b[i] ~ dnorm(0,taub);",
            "       logit(pbin[i]) <-  p[i] + b[i];")
  } else {
    bc <- c(bc,"       logit(pbin[i]) <-  p[i];")
  }
  bc <- c(bc,"       T[i]  ~ dcat(P[]);",
          "       r[i] ~ dbin(pbin[i],n[i]);",
          "   }")
  if (equal.variances) {
    bc <- c(bc,"sigma        <- 1/sqrt(tau);")
  } else {
    bc <- c(bc, paste("for (j in 1:",n.components,"){",sep=""))
    bc <- c(bc,"  sigma[j]     <- 1/sqrt(tau[j]);","}")
  }

  ## set up list of which variables to monitor

  monitor.var <- c("P","mu","sigma")
  if (random.effect) monitor.var <- c(monitor.var,"b","sigmab")
  monitor.var <- c(monitor.var,"T")
  
  ## return S3 class 

  res <- list(bugs.code=bc, n.components=n.components, monitor.var=monitor.var,
              ploidy.level=ploidy.level, random.effect=random.effect,
              equal.variances=equal.variances,
              E.segRatio=E.segRatio, type.parents=type, call=match.call())
  oldClass(res) <- "modelSegratioMM"
  return(res)
}

