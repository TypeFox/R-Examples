#' Range Cohesion Model for Ordered (and Non-spatial) Data
#' @description rangemod.1d takes observed site by species matrix and returns
#'               expected species richness values of each site
#' @param spmat a site by species matrix or data frame with species in columns
#' @param reps number of replicates
#' @param nb a neighbour object similar to one generated with
#'        '\code{\link[spdep]{poly2nb}}' of '\pkg{\link[spdep]{spdep}}'.
#'        If NULL then a list resembling object of class 'nb' is created.
#'        If NA then result is range scatter.
#' @param var an optional vector containing explanatory variable for
#'        constraining the randomization
#' @param first If TRUE, 'var' is used while choosing the first occurence as
#'        well.if 'var' is null, first is always set 'FALSE'
#' @param degen If true, each randomized site by species matrix is saved and
#'        provided in output
#' @param rsize which rangesizes to use for simulation, can be an integer
#'        vector of same length as number of species(collumns) or either
#'        'observed' or'unif'. See details for explanations
#' @details rangemod.1d impliments simulations used by Rahbeck et.al (2007) to
#'          data which are only in form of a site by species matrix and without
#'          any spatial information. A list similar to an nb object of spdep can
#'          prepared according to order in which the rows (sites) are arranged.
#'          A manually prepared list of neighbours for each site can also be
#'          used.It is important that each site must have at least one neighbour.
#'          'rsize' provides a vector of rangesizes.It can be 'unif' -
#'           ranges are drawn from a uniform distribution,between 1 to
#'          number of sites or 'observed' - rangesize
#'          of each species is exactly the same as in the observed matrix.
#'          Alternatively a it can also be a user specified integer vector, of
#'          same length as number of species.
#' @return If degen is FALSE, a data frame with four colums for mean, SD and
#'        confidence intervals of expected richness
#'
#' \itemize{
#'  \item{"mod.rich"}{ mean richness of each site}
#'  \item{"mod.sd"}{ standard deviation of species richness}
#'  \item{"q2.5"}{ lower limit of the confidence interval}
#'  \item{"q97.5"}{ upper limit of the confidence interval}
#' }
#'        If degen is TRUE, then a list containing above data frame and a list
#'        of all the randomized matrices
#' @references Rahbek, C., Gotelli, N., Colwell, R., Entsminger, G., Rangel, T.
#'              & Graves, G. (2007) Predicting continental-scale patterns of
#'              bird species richness with spatially explicit
#'              models. Proceedings of the Royal Society B: Biological Sciences,
#'              274, 165.
#'
#'              Gotelli, N.J., Anderson, M.J., Arita, H.T., Chao, A., Colwell,
#'              R.K., Connolly, S.R., Currie, D.J., Dunn, R.R., Graves, G.R. &
#'              Green, J.L. (2009) Patterns and causes of species richness:
#'              a general simulation model for macroecology. Ecology Letters,
#'              12, 873-886.
#' @examples
#' tempmat <- matrix(0,nrow=10,ncol=200,dimnames=list(letters[1:10],1:200))
#' tempmat <- as.matrix(apply(tempmat,2,function(x){rbinom(nrow(tempmat),1,
#'                      runif(1,0.1,1))}))
#' rownames(tempmat) <- letters[1:10]
#' temp <- rangemod.1d(tempmat,nb = NULL,var = NULL,rsize = "observed",reps = 5)
#' plot(temp[,1],ylim= c(min(temp[,1] -2),max(temp[,1]+2)),pch = 16,ylab = 'Species Richness')
#' segments(1:10,y0=temp[,1]-temp[,2],y1= temp[,1]+temp[,2])
#' @export
rangemod.1d <- function(spmat,nb = NULL,var = NULL,first = FALSE,degen = FALSE,
                        rsize = c("observed","unif"),reps){
  ####sanity check of arguments####
  if(!is.null(nb)){
    if(!is.na(nb)){
      if(!length(nb) == nrow(spmat)){
        stop("length(of 'nb' should be same as number of sites: ",
             length(nb)," and ", nrow(spmat))
      }
    }
  }

  if(!is.null(var)&& !length(var) == nrow(spmat)){
    stop("'var' should be of same length as number of sites: ",
         length(var)," and ",nrow(spmat),".")
  }

  if(is.vector(rsize,mode = "numeric")&& !length(rsize) == ncol(spmat)){
    stop("rsize should be of same length as number of species: ",
         length(rsize)," and ",ncol(spmat))
  }

  if(is.null(rownames(spmat))){
    warning("No rownames for 'spmat', setting rownames as 1:nrow(spmat)")
    rownames(spmat) <- 1:nrow(spmat)
  }
  ####chunk1- make input and output objects ####
  spmat[spmat>0] <- 1
  spmat <- as.matrix(spmat)
  keep <- which(colSums(spmat) > 0)
  spmat <- spmat[,keep]
  #check rangesizes
  if(is.vector(rsize,mode = "numeric")){
    range.size <- rsize[-remove]
  }else{
    rsize <- match.arg(rsize)
    range.size <- switch(rsize,observed = {colSums(spmat)},
                         unif = {sample(1:nrow(spmat),ncol(spmat),replace = T)})
  }
  ##sanity check for 'first'
  suppressWarnings(if(is.null(var)){
    first <- FALSE
  }else{first})

  mat.temp <- spmat
  mat.out <- matrix(nrow = nrow(spmat),ncol = reps,
                    dimnames = list(rownames(spmat),1:reps))
  degen.mats <- list()

  ####chunk2- get neighbour list####

  suppressWarnings(if(is.null(nb)){
    nblist <- list()
    nblist[[1]] <- 2
    nblist[[nrow(spmat)]] <- nrow(spmat)-1
    for(i in 2:(nrow(spmat)-1)){
      nblist[[i]] <- c(i-1,i+1)
    }
  }else{nblist <- nb})


  ####chunk3- spread each range in matrix####

  for(j in 1:reps){
    mat.temp[which(mat.temp > 0)] <- 0
    for(k in 1:length(range.size)){
      temp.vec1 <- random.range(uid = rownames(spmat),
                                nb=nblist,range.size = range.size[k],
                                var,first)
      mat.temp[which(rownames(mat.temp)%in%as.character(temp.vec1)),k] <- 1
    }
    mat.out[,j] <- rowSums(mat.temp)
    if(degen == TRUE){degen.mats[[j]] <- mat.temp}

  }
  ####output####

  out.df <- data.frame(mod.rich = apply(mat.out,1,mean),
                       mod.sd =  apply(mat.out,1,stats::sd),
                       q2.5 = apply(mat.out,1,stats::quantile,probs = 0.025),
                       q97.5 = apply(mat.out,1,stats::quantile,probs = 0.975))
    if(degen == TRUE){
    outlist <- list(out.df =  out.df,degenerate.matrices = degen.mats)
    outlist
  }else{out.df}

}
