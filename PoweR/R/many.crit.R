many.crit <- function(law.index, stat.indices, M = 10 ^ 3, vectn = c(20, 50, 100), levels = c(0.05, 0.1), alter = create.alter(stat.indices), law.pars = NULL, parstats = NULL, model = NULL, Rlaw = NULL, Rstats = NULL, center = FALSE, scale = FALSE) {


  if (any(stat.indices == 0) & is.null(Rstats)) stop("'Rstats' should be a list whose components are R functions.")
  if (any(stat.indices == 0)) {
    if (!is.list(Rstats)) stop("'Rstats' should be a list whose components are R functions.")
    for (i in 1:length(stat.indices)) if ((stat.indices[i] == 0) & !is.function(Rstats[[i]])) stop(paste("The ",i,"th component of 'Rstats' should be an R function",sep=""))
  }
  if (is.null(Rstats)) {
    Rstats <- as.list(1:length(stat.indices))
    for (i in 1:length(stat.indices)) Rstats[[i]] <- list(NULL)
  }


  Rcpp <- any(law.index == 0)
  
  if (is.function(Rlaw) & (law.index != 0)) stop("You should set 'law.index' to 0 when 'Rlaw' is a (random generating) function.")

  if(getRversion() < "3.1.0") dontCheck <- identity
  
  stats.len <- length(stat.indices)

# Management of alter
  if (!is.null(alter)) {
    if (!is.list(alter)) stop("'alter' should be a list")
    if (is.null(names(alter))) stop("'alter' should be a named list")
    if (any(is.na(names(alter)))) stop("'alter' names should all be defined")
    if (length(alter) != length(stat.indices)) stop("'alter' and 'stat.indices' should have the same length")
    for (s in 1:stats.len) {
        if (stat.indices[s] != 0) {
            if (names(alter)[s] != paste("stat",stat.indices[s],sep="")) stop(paste("Name of 'alter'[[",s,"]] should be equal to 'stat",stat.indices[s],sep=""))
            if (!(alter[[s]] %in% 0:4)) stop(paste("'alter'[[",s,"]] should be in  {0,1,2,3,4}.",sep=""))
            Cstat.name <- paste("stat",as.character(stat.indices[s]),sep="")
            alter.true <- .C(dontCheck(Cstat.name),as.double(0.0),1L,0.05,1L,rep(" ",50),1L,0.0,0L,0.0,0.0,0.0,0L,alter=as.integer(alter[s]),0L,rep(0.0,4),0L,PACKAGE="PoweR")$alter
            if (alter[[s]] != alter.true) {
                warning(paste("'alter'[[",s,"]] should be set to ",alter.true,". We have done this for you!"),sep="")
                alter[[s]] <- alter.true
            }
        }
    }
    alter <- unlist(alter)
  } else { # alter is NULL
    alter <- rep(0,stats.len)
    names(alter) <- paste("stat",stat.indices,sep="")
  }

 
# Management of parstats
  nbparstats <- rep(NA, length(stat.indices))
  nbparstats[stat.indices != 0] <- getnbparstats(stat.indices[stat.indices != 0])
  parstatstmp <- c()
  if (!is.null(parstats)) {
   if (!is.list(parstats)) stop("'parstats' should be a list")
   if (is.null(names(parstats))) stop("'parstats' should be a named list")
   if (any(is.na(names(parstats)))) stop("'parstats' names should all be defined")
   if (length(parstats) != length(stat.indices)) stop("'parstats' and 'stat.indices' should have the same length")
   for (s in 1:stats.len) {
       if (stat.indices[s] != 0) {
           if (names(parstats)[s] != paste("stat",stat.indices[s],sep="")) stop(paste("Name of 'parstats'[[",s,"]] should be equal to 'stat",stat.indices[s],sep=""))
           if (!is.na(parstats[[s]]) && (nbparstats[s] == 0)) stop(paste("'parstats[['",s,"]] should be equal to NA",sep=""))
           if ((nbparstats[s] != 0) && (length(parstats[[s]]) != nbparstats[s])) stop(paste("The length of parstats[[",s,"]] should be ",nbparstats[s],sep=""))
           parstatstmp <- c(parstatstmp, parstats[[s]])
       }
   }
} else {
    for (s in 1:stats.len) {
        if (stat.indices[s] != 0) {parstatstmp <- c(parstatstmp, stat.cstr(stat.indices[s])$stat.pars)}
    }
}
#  parstats <- parstatstmp
#  parstats[is.na(parstats)] <- 0
   parstatstmp <- parstatstmp[!is.na(parstatstmp)]
  

 
  nblevels <- length(levels)
  vectn.len <- length(vectn)

######################################
# debut modif ....

  stats.len <- length(stat.indices)
  decision.len <- stats.len * vectn.len * nblevels
  decision <- rep(0L, decision.len)
  critvalL <- rep(0.0, M * vectn.len * stats.len)
  nbparlaws <- length(law.pars)

  if (is.null(law.pars)) {
      tmp2 <- law.cstr(law.index)$law.pars
      law.pars <- c(tmp2, rep(0.0, 4 - length(tmp2)))
  }

  
  thetavec <- 0
  xvec <- 0
  p <- length(thetavec)
  np <- length(xvec)

  if (Rcpp | any(stat.indices == 0)) {
      out <- .Call("powcompfastRcpp", M = as.integer(M), law.index = as.integer(law.index), laws.len = 1L, vectn = as.integer(vectn), vectn.len = as.integer(vectn.len),
                   stat.indices = as.integer(stat.indices), stats.len = as.integer(stats.len), decision = as.integer(decision), decision.len = as.integer(decision.len),
                   levels = as.double(levels), nblevels = as.integer(nblevels), cL = as.double(critvalL), cR = as.double(0.0), usecrit = 0L, alter = as.integer(alter),
                   nbparlaws = as.integer(nbparlaws), parlaws = as.double(law.pars), nbparstats = as.integer(nbparstats), parstats = as.double(parstatstmp),
                   modelnum = 1L, funclist = list(function(){}), as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),
                   as.list(Rlaw), Rstats, as.integer(center), as.integer(scale), compquant = 1L, PACKAGE="PoweR",NAOK=TRUE)$cL
  } else {
        
      out <- .C("powcompfast", M = as.integer(M), law.index = as.integer(law.index), laws.len = 1L, vectn = as.integer(vectn), vectn.len = as.integer(vectn.len),
                stat.indices = as.integer(stat.indices), stats.len = as.integer(stats.len), decision = as.integer(decision), decision.len = as.integer(decision.len),
                levels = as.double(levels), nblevels = as.integer(nblevels), cL = as.double(critvalL), cR = as.double(0.0), usecrit = 0L, alter = as.integer(alter),
                nbparlaws = as.integer(nbparlaws), parlaws = as.double(law.pars), nbparstats = as.integer(nbparstats), parstats = as.double(parstatstmp),
                modelnum = 1L, funclist = list(function(){}), as.double(thetavec), as.double(xvec), as.integer(p), as.integer(np),
                as.integer(center), as.integer(scale), compquant= 1L, PACKAGE = "PoweR", NAOK = TRUE)$cL

  }
#####################################

  
 mylist <- as.list(c())
 for (s in 1:stats.len) {
 
   stat.index <- stat.indices[s]
   statname <- paste("stat",as.character(stat.index),sep="")
    
   mylist[[s]] <- matrix(NaN,nrow=length(vectn)*length(levels),ncol=2,dimnames=list(NULL,c("critL","critR")))
   
   for (l in 1:nblevels) {  
     for (n in 1:vectn.len) {
       if (alter[[s]] == 0) { # two.sided test

           res <- quantile(out[(s - 1) * M * vectn.len + (n - 1) * M + (1:M)], probs = c(levels[l] / 2, 1 - levels[l] / 2))
           mylist[[s]][n + vectn.len * (l - 1), ] <- res
           
             # mylist[[s]][n + vectn.len * (l - 1), ] <- compquant(n=vectn[n],law.index=law.index,stat.index=stat.index,probs=c(levels[l]/2,1-levels[l]/2),M=M,law.pars=law.pars,stat.pars=parstats[[s]],model=model,Rlaw=Rlaw,Rstat=Rstats[[s]],center=center,scale=scale)$quant
       } else if ((alter[[s]] == 1) || (alter[[s]] == 4)) { # less test, or bilateral test that reject H0 only for small values of the test statistic

           res <- quantile(out[(s - 1) * M * vectn.len + (n - 1) * M + (1:M)], probs = levels[l])
          
#         res <- compquant(n=vectn[n],law.index=law.index,stat.index=stat.index,probs=levels[l],M=M,law.pars=law.pars,stat.pars=parstats[[s]],model=model,Rlaw=Rlaw,Rstat=Rstats[[s]],center=center,scale=scale)$quant
         mylist[[s]][n + vectn.len * (l - 1), ] <- c(res,NA)
       } else if ((alter[[s]] == 2) || (alter[[s]] == 3)) { # greater test, or bilateral test that reject H0 only for large values of the test statistic

           res <- quantile(out[(s - 1) * M * vectn.len + (n - 1) * M + (1:M)], probs = 1 - levels[l])
           
#         res <- compquant(n=vectn[n],law.index=law.index,stat.index=stat.index,probs=1-levels[l],M=M,law.pars=law.pars,stat.pars=parstats[[s]],model=model,Rlaw=Rlaw,Rstat=Rstats[[s]],center=center,scale=scale)$quant
         mylist[[s]][n + vectn.len * (l - 1), ] <- c(NA,res)
       }
     }
   }

   # add a new column named "param" to distinguish the critical values from the same statistic with different parameters
   # if parstats is NULL, "param" will appear as NA
   mat <- expand.grid(vectn,levels,paste(parstats[[s]],collapse=" "))
   colnames(mat) <- c("n","level","param")
   mylist[[s]] <- cbind(mat,mylist[[s]])

   names(mylist)[s] <- statname
   
}

 list.names.replicated <- names(mylist)[table(names(mylist)) > 1]
 
 for (name in list.names.replicated) {

   mask <- names(mylist) %in% name

   names(mylist)[mask] <- paste(name,1:sum(mask),sep=".")
      
 }
  
 class(mylist) <- "critvalues"
 
 return(mylist)
 
 
}





























## First solution for creating vector t :
# nb <- c(2,2,3)
# t <- 1:sum(nb)
# toto <- c()
# for (s in 1:length(nb)) {
  
  # toto <- c(toto,rep(s,nb[s]))

# }
# t <- split(t,factor(toto))


## Second solution for creating vector t : 
# apply(cbind(cumsum(x) - x + 1,x),MARGIN=1,FUN=function(x){seq(from=x[1],length=x[2])})
   
# t <- apply(cbind(cumsum(nbparstats2) - nbparstats2 + 1,nbparstats2),MARGIN=1,FUN=function(x){seq(from=x[1],length=x[2])})
 
 
 
