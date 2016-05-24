
###############################################################################
# Qrefine:                                                                    #
#                                                                             #
# Refine the Q-matrix by minimizing the RSS. We estimate memberships using the#
# non-parametric classification method (weighted hamming).                    #
#                                                                             # 
# Inputs:                                                                     #
# (1) Y: a matrix of binary responses (1=correct, 0=incorrect). Rows          #
#        represent persons and columns represent items.                       #
# (2) Q: the Q-matrix of the test. Rows represent items and columns represent #
#        attributes.                                                          #
# (3) gate: "AND", the examinee needs to possess all related attributes       #
#                  in order to answer an item correctly                       #
#           "OR", the examinee needs to possess only one of the related       #
#                 attributes in order to answer an item correctly             #
# (4) max.ite: the iteration to run until all RSS of all items are stationay. #
#                                                                             #
# Definitions of variables/matrices:                                          #
#  RSS: the vector collecting the item and overall RSSs acorss all iterations #
#  max.rss.item: the item with largest RSS                                    #
#  update.q: the q-vector chosen to replace the old q-vector                  #
#                                                                             #
# Outputs:                                                                    #
# (1) initial.class: initial classification                                   #
# (2) terminal.class: terminal classification                                 #
# (3) modified.Q: the modified Q-matrix                                       #
# (4) modified.entries: the modified q-entries                                #
#                                                                             #
###############################################################################

Qrefine <- function(Y, Q, gate=c("AND", "OR"), max.ite=50)
{
  gate <- match.arg(gate)
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  initial.Q <- Q
  nitem <- dim(Q)[1]
  natt <- dim(Q)[2]
  M <- 2 ^ natt
  index <- NULL
  RSS <- NULL
  
  pattern <- AlphaPermute(natt)
  pattern.Q <- pattern[-1,]
  
  # We may need to run the searching a couple of times until all RSS of all items are stationay.
  
  for (m in 1:max.ite)
  {
    cat(paste("Iteration:", m))
    cat("\n")
    
    max.rss.item <- NULL
    
    for (k in 1:nitem) 
    {
      classification <- AlphaNP(Y, Q, gate, method="Hamming")
      est.class <- classification$est.class
      est.ideal <- classification$est.ideal
      
      # Record the initial classification
      
      if (m==1 & k==1) initial.class <- est.class
      
      # Update Q by minimizing RSS 
      
      diff <- Y - est.ideal
      rss <- apply(diff ^ 2, 2, sum)
      RSS <- rbind(RSS, c(rss, sum(rss)))
      
      # Start the algorithm with the item of largest RSS. Each item is visited only once.
      
      max.rss <- {if (k == 1) which(rss == max(rss)) else which(rss == max(rss[-max.rss.item]))}
      if (length(max.rss) > 1) max.rss <- sample(max.rss, 1) 
      max.rss.item <- c(max.rss.item, max.rss)
      
      # Find the q-vector among the 2^K possible vectors that the corresponding ideal responses yield  
      # minimal RSS with data of item of max RSS. 
      
      update.rss <- NULL
      
      for (i in 1:(M - 1))
      {
        u <- apply(t(pattern[est.class, ]) ^ pattern.Q[i, ], 2, prod)
        temp.rss <- sum((Y[, max.rss] - u) ^ 2)
        update.rss <- c(update.rss, temp.rss)
      }
      
      min.update.rss <- which(update.rss == min(update.rss))
      if (length(min.update.rss) > 1) min.update.rss <- sample(min.update.rss, 1)
      update.q <- pattern.Q[min.update.rss, ]
      Q[max.rss, ] <- update.q
    }
    
    # Stopping criterion: if all of the RSSs between two iterations are identical (which means
    # it's not possible to improve), break the loop.                                        

    if (sum(abs(RSS[((m - 1) * nitem + 1):(m * nitem), nitem + 1] - rep(RSS[((m - 1) * nitem + 1), nitem + 1], nitem))) == 0) break else m <- m + 1
  }

  # Record the terminal classification

  terminal.class <- est.class

  # Report the modified entries

  if (sum((Q - initial.Q) ^ 2) == 0) modified <- "NA" else 
  {
    modified <- which((Q - initial.Q) != 0, arr.ind=TRUE)
    colnames(modified) <- c("Item", "Attribute")
  }
  
  output <- list(modified.Q=Q, modified.entries=modified, initial.class=initial.class, terminal.class=terminal.class, patterns=pattern, initial.Q=initial.Q, Y=Y, gate=gate)
  class(output) <- "Qrefine"
  return(output)
}


