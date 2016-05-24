.RMXE.th <- function(th, PFam, modifyfct, loRad = 0, upRad = Inf, z.start = NULL,
             A.start = NULL, upper = NULL, lower = NULL,
             OptOrIter = "iterate", maxiter = 50,
             tol = .Machine$double.eps^0.4, loRad0 = 1e-3, ...){
      PFam <- modifyfct(th,PFam)
      IC <- radiusMinimaxIC(L2Fam=PFam, neighbor= ContNeighborhood(),
                            risk = asMSE(), verbose = FALSE,
                            loRad = loRad, upRad = upRad, z.start = z.start,
                            A.start = A.start, upper = upper, lower = lower,
                            OptOrIter = OptOrIter, maxiter = maxiter,
                            tol = tol, warn = FALSE,
                            loRad0 = loRad0, returnNAifProblem = TRUE)
      if(!is(IC,"IC")) if(is.na(IC)) return(NA)
      txt <- "least favorable radius:"
      wL <- grepl(txt, Infos(IC)[,"message"])
      rad <- as.numeric(gsub(txt, "", Infos(IC)[wL,"message"]))
      return(list(b=clip(IC), a=cent(IC), a.w = cent(weight(IC)),
                           A=stand(IC),  A.w = stand(weight(IC)), rad=rad))
}

.MBRE.th <- function(th, PFam, modifyfct,
             z.start = NULL, A.start = NULL, upper = 1e4,
             lower = 1e-4, OptOrIter = "iterate",
             maxiter = 50, tol = .Machine$double.eps^0.4, ...){
      PFam <- modifyfct(th,PFam)
      RobM <- InfRobModel(center = PFam, neighbor = ContNeighborhood(radius = 6))
      IC <- optIC(model = RobM, risk = asBias(), verbose = FALSE,
             z.start = z.start, A.start = A.start, upper = upper,
             lower = lower, OptOrIter = OptOrIter,
             maxiter = maxiter, tol = tol, warn = TRUE, noLow = FALSE,
             .withEvalAsVar = FALSE, returnNAifProblem = TRUE)
      if(!is(IC,"IC")) if(is.na(IC)) return(NA)
      mA <- max(stand(IC))
      mAw <- max(stand(weight(IC)))
      return(list(b=clip(IC), a=cent(IC), a.w=cent(weight(IC)),
               A=stand(IC)/mA, A.w=stand(weight(IC))/mAw))
}

.OMSE.th <- function(th, PFam, modifyfct, radius = 0.5,
             z.start = NULL, A.start = NULL, upper = 1e4,
             lower = 1e-4, OptOrIter = "iterate",
             maxiter = 50, tol = .Machine$double.eps^0.4, ...){
      PFam <- modifyfct(th,PFam)
      RobM <- InfRobModel(center = PFam,
                          neighbor = ContNeighborhood(radius = radius))
      IC <- optIC(model = RobM, risk = asMSE(), verbose = FALSE,
             z.start = z.start, A.start = A.start, upper = upper,
             lower = lower, OptOrIter = OptOrIter,
             maxiter = maxiter, tol = tol, warn = TRUE, noLow = FALSE,
             .withEvalAsVar = FALSE, returnNAifProblem = TRUE)
      if(!is(IC,"IC")) if(is.na(IC)) return(NA)
      res=list(b=clip(IC), a=cent(IC), a.w = cent(weight(IC)),
                A=stand(IC), A.w = stand(weight(IC)))
      return(res)
}


.getLMGrid <- function(thGrid, PFam, optFct = .RMXE.th, modifyfct, radius = 0.5,
                       GridFileName="LMGrid.Rdata", withPrint = FALSE,
                       upper = 1e4, lower = 1e-4, OptOrIter = "iterate",
                       maxiter = 50, tol = .Machine$double.eps^0.4,
                       loRad = 0, upRad = Inf, loRad0 = 1e-3,
                       loRad.s=0.2, up.Rad.s=1,
                       withStartLM = TRUE
                       ){
   wprint <- function(...){ if (withPrint) print(...)}
   thGrid <- unique(sort(thGrid))
   lG  <- length(thGrid)
   lG2 <- lG%/%2
   olG <- c(lG2:1,(lG2+1):lG)
   thGrid <- thGrid[olG]
   itLM <- 0
   z1 <- z.start <- NULL
   A1 <- A.start <- NULL
   r1l <- r.start.l <- NULL
   r1u <- r.start.u <- NULL
   getLM <- function(th){
               itLM <<- itLM + 1
               if(withPrint) cat("Evaluation Nr.", itLM," at th = ",th,"\n")
               a <- try(
               optFct(th = th, PFam = PFam, modifyfct = modifyfct,
                      z.start = z.start, A.start = A.start,
                      upper = upper, lower = lower, OptOrIter = OptOrIter,
                       maxiter = maxiter, tol = tol,
                       loRad = loRad, upRad = upRad, loRad0 = loRad0,
                       loRad.s = r.start.l, upRad.s = r.start.u),
               silent=TRUE)
               print(a)
               print(A.start)
               print(z.start)
               print(c(r.start.l,r.start.u))
               if(is(a,"try-error")|any(is.na(a))){ a <- rep(NA,13)}else{
                  if(withStartLM){
                     if(itLM==1){
                        z1 <<- a[["a.w"]]
                        A1 <<- a[["A"]]
                        if(!is.null(a$rad)){
                           r1l <<- max(a[["rad"]]/1.3,loRad)
                           r1u <<- min(a[["rad"]]*1.3,upRad)
                        }
                     }
                     z.start <<- a[["a.w"]]
                     A.start <<- a[["A"]]
                     if(!is.null(a$rad)){
                        r.start.l <<- max(a[["rad"]]/1.3,loRad)
                        r.start.u <<- min(a[["rad"]]*1.3,upRad)
                     }
                     if(itLM==lG2){
                        z.start <<- z1
                        A.start <<- A1
                        r.start.l <<- r1l
                        r.start.u <<- r1u
                     }
                     a <- c(a[["b"]],a[["a"]],a[["a.w"]],a[["A"]],a[["A.w"]])
                  }
               }
               print(a)
               return(a)
            }

   distroptions.old <- distroptions()
   distrExOptions.old <- distrExOptions()
   distroptions("withgaps"=FALSE)
   distrExOptions( MCIterations=1e6,
                   GLIntegrateTruncQuantile=.Machine$double.eps,
                   GLIntegrateOrder=1000,
                   ElowerTruncQuantile=1e-7,
                   EupperTruncQuantile=1e-7,
                   ErelativeTolerance = .Machine$double.eps^0.4,
                   m1dfRelativeTolerance = .Machine$double.eps^0.4,
                   m2dfRelativeTolerance = .Machine$double.eps^0.4,
                   nDiscretize = 300, IQR.fac = 20)
   on.exit({do.call(distrExOptions,args=distrExOptions.old)
            do.call(distroptions,args=distroptions.old)
            })
   LMGrid <- t(sapply(thGrid,getLM))

   iNA <- apply(LMGrid,1, function(u) any(is.na(u)))
   LMGrid <- LMGrid[!iNA,,drop=FALSE]
   thGrid <- thGrid[!iNA]
   oG <- order(thGrid)
   thGrid <- thGrid[oG]
   LMGrid <- LMGrid[oG,,drop=FALSE]
   Grid <- cbind(xi=thGrid,LM=LMGrid)

   if(GridFileName!="") save(Grid, file=GridFileName)
   wprint(Grid)
   return(Grid)
}


.saveGridToCSV <- function(Grid, toFileCSV, namPFam, nameInSysdata){
   write.table(format(Grid,digits=21),
               file=toFileCSV, row.names=FALSE, col.names=FALSE)
   toFileTXT <- gsub("(.+\\.)csv$","\\1txt",toFileCSV)
   cat(file=toFileTXT,gsub(" ","",namPFam),"\n",nameInSysdata)
   return(invisible(NULL))
}

.readGridFromCSV <- function(fromFileCSV){
  rg <- read.table(fromFileCSV, colClasses=rep("character",2), sep=" ", header=FALSE)
  nrg <- nrow(rg)
  Grid <- matrix(as.numeric(as.matrix(rg)),nrow=nrg)

  as.matrix(read.csv(fromFileCSV)); dimnames(Grid) <- NULL
  fromFileTXT <- gsub("(.+\\.)csv$","\\1txt",fromFileCSV)
  res2 <- scan(file=fromFileTXT, what=c("character","character"))
  return(list(Grid=Grid, namPFam=res2[1], namInSysdata=res2[2]))
}

.generateInterpGrid <- function(thGrid, PFam, toFileCSV = "temp.csv",
            getFun = .getLMGrid, ..., modifyfct, nameInSysdata,
            GridFileName, withPrint = TRUE){
  if(missing(GridFileName))
     GridFileName <- paste(gsub("^\\.(.+)","\\1",nameInSysdata),".Rdata",sep="")
  Grid <- getFun(thGrid = thGrid, PFam = PFam, ..., modifyfct = modifyfct,
                 withPrint = withPrint, GridFileName = GridFileName)
  .saveGridToCSV(Grid,toFileCSV,name(PFam),nameInSysdata)
  return(invisible(NULL))
}

