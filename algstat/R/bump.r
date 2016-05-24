#' Convert Dimensions of Approval Data
#'
#' Convert dimensions of approval data, see details for more.
#'
#'In a survey in which k objects are approved from a list of n objects by N voters, the survey responses can be summarized with choose(n, k) frequencies.  bump can summarize these frequencies by computing the number of votes for any lower order grouping primarily using the Tmaker function.  We refer to this as "bumping down".  Given a ith summary (the number of votes for each i-grouping), we can compute the expected (i+1)-group votes by evenly distributing each i-groups votes to the (i+1)-groups containing the i-groups and summing over each i-group's contribution.  This is "bumping up".
#'
#' See examples for the cookie example.  
#'
#'As a simple example of bumping up, suppose we have a survey in which 100 individuals select their favorite 3 items out of 6 items.  The total number of votes cast is therefore 100*3 = 300.  
#'
#' If that is all that is known, then the (one) bumped up expected dataset would simply expect each of the 6 items to be listed equally frequently: a vector of length 6 with each element equal to 300/6 = 50.  We would expect each of the 15 pairings to have 300/choose(6, 2) = 300 / 15 = 20 votes, and each of the 20 triples to have 300/choose(6, 3) = 5 votes.
#'
#' Now suppose we learn that the six objects were voted for 30, 40, 50, 50, 60, and 70 times, respectively, but nothing more.  Obviously, we could then compute 300 votes had been cast (V0 = 100 voters), but "looking up" we could also guess that the pairing [12] was voted for 30/choose(5, 1) + 40/choose(5, 1) = 14 times.  The reasoning is that if object 1 were observed 30 times and 2 40 times, if each were paired evenly with each of the others 1 would contribute 30/choose(5, 1) votes to the pairing and 2 40/choose(5, 1).  The choose(5, 1) corresponds to the number of pairings that 1 (or 2) is present in: [12], [13], [14], [15], [16].  The same idea can be used to estimate the number of votes for each of the choose(6, 2) = 15 pairs.  (See examples.)  This is bumping up; for any level of summary, we can bump it up as far as we like (all the way up to the set itself).
#'
#' Bumping down is easier.  The only thing needed to know is that it follows the order of the subsets function.  For example, in the above voting example, the 15 pairs votes are assumed to be in the order subsets(6, 2), and the result is given in the order of subsets(6, 1).
#'
#' If method = "even", exactly the above is done.  If method = "popular", then when bumping up the number of votes for [12] isn't determined by each of 1 and 2 donating their samples evenly to each of their upstream pairs; but rather 1 and 2 donating to each other (as contributors of the pair [12]) according to how popular it is amongst the alternatives.  In other words, 1 is thought to have (in this case) 30 votes to "give" to either 2, 3, 4, 5, or 6.  If method = "even", it donates 30/5 to each.  If method = "popular", it donates 40/(40+50+50+60+70) of the 30 votes to 2 (as a contributor of [12]), 50/(40+50+50+60+70) of the 30 to 3, and so on.  The expected frequency of [12] is therefore made up as the sum of such contributions from each of the places the contribution might come.  Here, the contributors of [12] are [1] and [2], with contributions 30*40/(40+50+50+60+70) and 40*30/(30+50+50+60+70), for a total of 9.06 expected [12] votes.  The same is done for higher order votes; e.g. [123] takes even or popular contributions from [12], [13], and [14]. 
#' 
#' @param x the summary given
#' @param n the number of objects in the running
#' @param k the number of objects approved
#' @param vin the level of summary given
#' @param vout the level of summary/expectation desired
#' @param method "popular" (default) or "even", see details
#' @return ...
#' @export bump
#' @examples
#'
#' \dontrun{
#' 
#' V0 <- 100 # V0 = number of voters (not votes)
#' bump(V0, 6, 3, 0, 0) # no bump
#' bump(V0, 6, 3, 0, 1) # 1-up
#' bump(V0, 6, 3, 0, 2) # 2-up
#' bump(V0, 6, 3, 0, 3) # 3-up
#' 
#' V1 <- c(30, 40, 50, 50, 60, 70)
#' bump(V1, 6, 3, 1, 0) # bump down
#' bump(V1, 6, 3, 1, 1) # no bump
#' bump(V1, 6, 3, 1, 2) # 1-up
#' bump(V1, 6, 3, 1, 3) # 2-up
#' 
#' cbind(
#'   bump(V1, 6, 3, 1, 2, "popular"),
#'   bump(V1, 6, 3, 1, 2, "even")
#' )
#' 
#' 
#' 
#' 
#' 
#' data(cookie)
#' (out <- spectral(cookie$freq, 6, 3, cookie$cookies))
#' 
#' (V0 <- out$obs$V0)
#' bump(V0, 6, 3, 0, 0)
#' bump(V0, 6, 3, 0, 1)
#' bump(V0, 6, 3, 0, 2)
#' bump(V0, 6, 3, 0, 3)
#' out$fullExp$V0
#' out$decompose(out$effects[,1])
#' 
#' (V1 <- out$obs$V1)
#' bump(V1, 6, 3, 1, 0) # cbind(bump(V1, 6, 3, 1, 0), out$fullExp$V1[[1]])
#' bump(V1, 6, 3, 1, 1) # cbind(bump(V1, 6, 3, 1, 1), out$fullExp$V1[[2]])
#' bump(V1, 6, 3, 1, 2) # cbind(bump(V1, 6, 3, 1, 2), out$fullExp$V1[[3]])
#' bump(V1, 6, 3, 1, 3) # cbind(bump(V1, 6, 3, 1, 3), out$fullExp$V1[[4]])
#' out$fullExp$V1 # the sampler doesn't distribute it's samples up evenly
#' 
#' (V2 <- out$obs$V2)
#' bump(V2, 6, 3, 2, 0) # cbind(bump(V2, 6, 3, 2, 0), out$fullExp$V2[[1]])
#' bump(V2, 6, 3, 2, 1) # cbind(bump(V2, 6, 3, 2, 1), out$fullExp$V2[[2]])
#' bump(V2, 6, 3, 2, 2) # cbind(bump(V2, 6, 3, 2, 2), out$fullExp$V2[[3]])
#' bump(V2, 6, 3, 2, 3) # cbind(bump(V2, 6, 3, 2, 3), out$fullExp$V2[[4]])
#' 
#' (V3 <- out$obs$V3)
#' bump(V3, 6, 3, 3, 0)
#' bump(V3, 6, 3, 3, 1)
#' bump(V3, 6, 3, 3, 2)
#' bump(V3, 6, 3, 3, 3)
#' 
#'
#'
#' }
#'
bump <- function(x, n, k, vin, vout, method = c("popular", "even")){
	
  method <- match.arg(method)
  
  if(length(x) != choose(n, vin)) stop(
    paste0("n = ", n, " and vin = ", vin, 
      " suggests the length of x is ", choose(n, vin), 
      ", not ", length(x), "."),
    call. = FALSE
  )
	
  nvinVotes  <- sum(x)
  nVoters    <- nvinVotes / choose(k, vin)
  nvoutVotes <- choose(k, vout) * nVoters
  
  if(vin == vout) out <- unname(x)

  if(vin < vout){
  	if(length(x) == 1){ 
  	  out <- nvoutVotes * rep(1/choose(n, vout), choose(n, vout))
  	} else {
  	  if(method == "even"){
      	message("Bumping up with even allocation...")  	  	
        upMat <- Emaker(n, vin, vout)
        out <- (upMat / sum(upMat[,1])) %*% x 
        out <- (nvoutVotes/sum(out)) * out        
      } else if(method == "popular"){
      	message("Bumping up with popular allocation...")
      	for(vinEff in vin:(vout-1)){
          upMat <- Emaker(n, vinEff, vinEff+1)
          contribsMat <- 0 * upMat
          for(i in 1:choose(n, vinEff)){
            # find all other same-level terms that contribute to 
            # next-level terms
            v <- upMat[,i]
            ndcsContribsUp <- which(as.logical(v))
          
            # extract those rows (the submatrix)
            upMatRelated <- upMat[ndcsContribsUp,]
          
            # find x-indices that correspond to same-level contributors
            upMatRelated[,i] <- 0L
            ndcsContribsSide <- which(colSums(upMatRelated) > 0)
          
            # allocate x[i]'s contributions
            contribDist <- 
              (x[ndcsContribsSide] / sum(x[ndcsContribsSide])) * x[i]
          
            # put the updated info in the matrix
            contribsMat[ndcsContribsUp,-i][upMat[ndcsContribsUp,-i] == 1] <- 
              contribsMat[ndcsContribsUp,-i][upMat[ndcsContribsUp,-i] == 1] + 
              contribDist
          }
          x <- rowSums(contribsMat)
        }
        out <- x 
        out <- (nvoutVotes/sum(out)) * out
      }
    }
  }

  if(vin > vout){
    message("Bumping down...")  	  	  	
    out <- Tmaker(n, vin, vout) %*% unname(as.numeric(x))
    out <- (nvoutVotes/sum(out)) * out
  }

  as.numeric(out)
}



























vec2subsetNumber <- function(v, n, k){
  sets <- sapply(subsets(n, k), paste, collapse = " ")
  which(paste(v, collapse = " ") == sets)
}

subsetNumber2vec <- function(x, n, k) subsets(n, k)[[x]]

is.subset <- function(x, y) all(x %in% y)

containsQ <- function(y, x) all(x %in% y)




