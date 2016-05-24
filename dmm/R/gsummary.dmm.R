gsummary.dmm <-
function(dmmobj,traitset="all",componentset="all",bytrait=T,gls=F,digits=3, ...)
# gsummary.dmm() - make summary tables of genetic parameters for a dmm object
{
   if(traitset[1] == "all"){
    traits <- dimnames(dmmobj$b)[[2]][1:ncol(dmmobj$b)]
  }
  else {
    traits <- traitset
  }
  traitpairs <- permpaste(traits)
  l <- length(traits)
  alltraitpairs <- permpaste(dimnames(dmmobj$b)[[2]])
 
  if(componentset[1] == "all") {
    components <- dimnames(dmmobj$variance.components)[[1]]
  }
  else {
    components <- componentset
  }
  coefs <- dimnames(dmmobj$b)[[1]]

  if(bytrait) {
    ftables <- vector("list",l)  # one table per trait
    count <- 0
    for(j in traits) {
      count <- count + 1
      ci95lo <- dmmobj$fraction[components,j] - 1.96 * dmmobj$fraction.se[components,j]
      ci95hi <- dmmobj$fraction[components,j] + 1.96 * dmmobj$fraction.se[components,j]
      ftable <- data.frame(Trait=j, Estimate=dmmobj$fraction[components,j],
                   StdErr=dmmobj$fraction.se[components,j],CI95lo=ci95lo,CI95hi=ci95hi,
                   row.names=components)
      # save ftables as one element of a list of tables
      ftables[[count]] <- ftable
    }

    rtables <- vector("list",l*l)  # one table per traitpair
    count <- 0
    for(i in traits) {
      for(j in traits) {
        traitpair <- paste(i,":",j,sep="",collapse=NULL)
        ij <- match(traitpair,alltraitpairs)
        count <- count + 1
        ci95lo <- dmmobj$correlation[components,ij] - 1.96 * dmmobj$correlation.se[components,ij]
        ci95hi <- dmmobj$correlation[components,ij] + 1.96 * dmmobj$correlation.se[components,ij]
        rtable <- data.frame(Traitpair=alltraitpairs[ij],
                     Estimate=dmmobj$correlation[components,ij],
                     StdErr=dmmobj$correlation.se[components,ij],
                     CI95lo=ci95lo,CI95hi=ci95hi,row.names=components)
        rtables[[count]] <- rtable
      }
    }
  }
  else {  # not bytrait
    ftables <- vector("list",length(components))  # one table per component
    count <- 0
    for(j in components) {
      count <- count + 1
      ci95lo <- dmmobj$fraction[j,traits] - 1.96 * dmmobj$fraction.se[j,traits]
      ci95hi <- dmmobj$fraction[j,traits] + 1.96 * dmmobj$fraction.se[j,traits]
      ftable <- data.frame(Component=j, Estimate=dmmobj$fraction[j,traits],
                   StdErr=dmmobj$fraction.se[j,traits],CI95lo=ci95lo,CI95hi=ci95hi,
                   row.names=traits)
      # save ftables as one element of a list of tables
      ftables[[count]] <- ftable
    }

    rtables <- vector("list",length(components))  # one table per component
    count <- 0
    for(i in components) {
       count <- count + 1
       ci95lo <- dmmobj$correlation[i,traitpairs] - 1.96 * dmmobj$correlation.se[i,traitpairs]
       ci95hi <- dmmobj$correlation[i,traitpairs] + 1.96 * dmmobj$correlation.se[i,traitpairs]
       rtable <- data.frame(Component=i,
                    Estimate=dmmobj$correlation[i,traitpairs],
                    StdErr=dmmobj$correlation.se[i,traitpairs],
                    CI95lo=ci95lo,CI95hi=ci95hi,row.names=traitpairs)
       rtables[[count]] <- rtable
    }
  }

  ptables <- vector("list",1)
  ci95lo <- as.vector(dmmobj$phenotypic.variance[traits,traits]) - 1.96 * as.vector(dmmobj$phenotypic.variance.se[traits,traits])
  ci95hi <- as.vector(dmmobj$phenotypic.variance[traits,traits]) + 1.96 * as.vector(dmmobj$phenotypic.variance.se[traits,traits])
  ptable <- data.frame(Traitpair=traitpairs,
                       Estimate=as.vector(dmmobj$phenotypic.variance[traits,traits]),
                       StdErr=as.vector(dmmobj$phenotypic.variance.se[traits,traits]),
                       CI95lo=ci95lo,CI95hi=ci95hi)
  ptables[[1]] <- ptable

  retobj <- list(ftables=ftables, rtables=rtables, ptables=ptables, traits=traits, components=components, bytrait=bytrait, gls=gls, digits=digits)

  if(gls) {
  if(bytrait) {
    gftables <- vector("list",l)  # one table per trait
    count <- 0
    for(j in traits) {
      count <- count + 1
      ci95lo <- dmmobj$gls$fraction[components,j] - 1.96 * dmmobj$gls$fraction.se[components,j]
      ci95hi <- dmmobj$gls$fraction[components,j] + 1.96 * dmmobj$gls$fraction.se[components,j]
      ftable <- data.frame(Trait=j, Estimate=dmmobj$gls$fraction[components,j],
                StdErr=dmmobj$gls$fraction.se[components,j],CI95lo=ci95lo,CI95hi=ci95hi,
                row.names=components)
      # save ftables as one element of a list of tables
      gftables[[count]] <- ftable
    }

    grtables <- vector("list",l*l)  # one table per traitpair
    count <- 0
    for(i in traits) {
      for(j in traits) {
        traitpair <- paste(i,":",j,sep="",collapse=NULL)
        ij <- match(traitpair,alltraitpairs)
        count <- count + 1
        ci95lo <- dmmobj$gls$correlation[components,ij] - 1.96 * dmmobj$gls$correlation.se[components,ij]
        ci95hi <- dmmobj$gls$correlation[components,ij] + 1.96 * dmmobj$gls$correlation.se[components,ij]
        rtable <- data.frame(Traitpair=alltraitpairs[ij],
                     Estimate=dmmobj$gls$correlation[components,ij],
                     StdErr=dmmobj$gls$correlation.se[components,ij],
                     CI95lo=ci95lo,CI95hi=ci95hi,row.names=components)
        grtables[[count]] <- rtable
      }
    }
  }
  else {  # not bytrait
    gftables <- vector("list",length(components))  # one table per component
    count <- 0
    for(j in components) {
      count <- count + 1
      ci95lo <- dmmobj$gls$fraction[j,traits] - 1.96 * dmmobj$gls$fraction.se[j,traits]
      ci95hi <- dmmobj$gls$fraction[j,traits] + 1.96 * dmmobj$gls$fraction.se[j,traits]
      ftable <- data.frame(Component=j, Estimate=dmmobj$gls$fraction[j,traits],
                  StdErr=dmmobj$gls$fraction.se[j,traits],CI95lo=ci95lo,CI95hi=ci95hi,
                  row.names=traits)
      # save ftables as one element of a list of tables
      gftables[[count]] <- ftable
    }

    grtables <- vector("list",length(components))  # one table per component
    count <- 0
    for(i in components) {
       count <- count + 1
       ci95lo <- dmmobj$gls$correlation[i,traitpairs] - 1.96 * dmmobj$gls$correlation.se[i,traitpairs]
       ci95hi <- dmmobj$gls$correlation[i,traitpairs] + 1.96 * dmmobj$gls$correlation.se[i,traitpairs]
       rtable <- data.frame(Component=i,
                    Estimate=dmmobj$gls$correlation[i,traitpairs],
                    StdErr=dmmobj$gls$correlation.se[i,traitpairs],
                    CI95lo=ci95lo,CI95hi=ci95hi,row.names=traitpairs)
       grtables[[count]] <- rtable
    }
  }

  gptables <- vector("list",1)
  ci95lo <- as.vector(dmmobj$gls$phenotypic.variance[traits,traits]) - 1.96 * as.vector(dmmobj$gls$phenotypic.variance.se[traits,traits])
  ci95hi <- as.vector(dmmobj$gls$phenotypic.variance[traits,traits]) + 1.96 * as.vector(dmmobj$gls$phenotypic.variance.se[traits,traits])
  ptable <- data.frame(Traitpair=traitpairs,
                       Estimate=as.vector(dmmobj$gls$phenotypic.variance[traits,traits]),
                       StdErr=as.vector(dmmobj$gls$phenotypic.variance.se[traits,traits]),
                       CI95lo=ci95lo,CI95hi=ci95hi)
  gptables[[1]] <- ptable

  retobj <- list(ftables=ftables,rtables=rtables,ptables=ptables,gftables=gftables,grtables=grtables,gptables=gptables,traits=traits, components=components, bytrait=bytrait, gls=gls, digits=digits)
  }

  retobj$call <- match.call()
  class(retobj) <- "gsummary.dmm"
  return(retobj)
}
