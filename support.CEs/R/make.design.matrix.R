make.design.matrix <- function(choice.experiment.design, optout = TRUE, 
                               categorical.attributes = NULL,
                               continuous.attributes = NULL,
                               unlabeled = TRUE, common = NULL, binary = FALSE) 
{
# Name: make.design.matrix
# Title: Converting a choice experiment design into a design matrix
# Arguments:
#  choice.experiment.design   A data frame containing a choice experiment design created
#                               by the function Lma.design() or rotation.design().
#  optout                     A logical variable describing whether or not the opt-out alternative
#                               is included in the design matrix created by this function.
#  categorical.attributes     A vector containing the names of attributes treated as categorical
#                               independent variables in the analysis.
#  continuous.attributes      A vector containing the names of attributes treated as continuous
#                               independent variables in the analysis.
#  unlabeled                  A logical variable describing the types of a choice experiment
#                               design assigned by the argument choice.experiment.design.
#  common                     A vector containing a fixed combination of attribute-levels 
#                               corresponding to a common base option in each question.
#  binary                     When the function is applied to the conditional logit model, 
#                               the argument is set as FALSE. When the function is applied to
#                               the binary choice model, it is set as TRUE.



# assign design information

  nblocks <- choice.experiment.design$design.information$nblocks
  nquestions <- choice.experiment.design$design.information$nquestions
  nalternatives <- choice.experiment.design$design.information$nalternatives
  nattributes <- choice.experiment.design$design.information$nattributes

# set variables

  variable.names <- NULL
  ced <- choice.experiment.design$alternatives
  
  if ((unlabeled == TRUE) && (is.null(common) == FALSE)) {
    nalternatives <- nalternatives + 1
    common.base <- ced[[1]]
    common.base[, 3] <- nalternatives
    for (i in categorical.attributes) {
      common.base[i] <- factor(common[[i]],
                               levels = levels(common.base[[i]]))
    }
    for (i in continuous.attributes) {
      common.base[i] <- as.numeric(common[[i]])
    }
    ced <- c(ced, alt.common = list(common.base))
  }

  conv.ced <- vector("list", nalternatives)

# create attribute variables

  for (i in 1:nalternatives) {

    # categorical attribute variables
    if (is.null(categorical.attributes) == FALSE) {
      for (j in 1:length(categorical.attributes)) {
        k <- which(names(ced[[i]]) == categorical.attributes[j])
        m <- nlevels(ced[[i]][, k])
        variable.names <- c(variable.names, 
                            levels(ced[[i]][, k])[2:m])
        conv.ced[[i]] <- cbind(conv.ced[[i]],
                               model.matrix(~ ced[[i]][, k] - 1)[, 2:m])
      }
    }

    # continuous attribute variables
    if (is.null(continuous.attributes) == FALSE) {
      for (j in 1:length(continuous.attributes)) {
        k <- which(names(ced[[i]]) == continuous.attributes[j]) 
        variable.names <- c(variable.names, names(ced[[i]])[k])
        conv.ced[[i]] <- cbind(conv.ced[[i]],
                               as.numeric(as.character(ced[[i]][, k])))
      }
    }
  }

# create design matrix

  nvariables <- length(variable.names) / nalternatives

  # multinomial choice 
  if (binary == FALSE) {
    # unlabeled
    if (unlabeled == TRUE) {
      my.design <- conv.ced[[1]]
      for (i in 2:nalternatives) {
        my.design <- rbind(my.design, conv.ced[[i]])
      }
      colnames(my.design) <- variable.names[1: nvariables]
    }
    # labled
    else {
      my.design <- diag.block(conv.ced)
      variable.names <- paste(variable.names,
                              rep(1:nalternatives, each=nvariables),
                              sep = "")
      colnames(my.design) <- variable.names
    }

    # create BLOCK, QES, and ALT variables
    BQS <- ced[[1]][, 1:3]
    for (i in 2:nalternatives) {
        BQS <- rbind(BQS, ced[[i]][, 1:3])
    }

    # ASC for unlabeled design
    if (unlabeled == TRUE) {
      # with opt-out option
      if (optout == TRUE) {
        ASC <- rep(1, nalternatives * nquestions * nblocks)
      }
      # without opt-out option
      else {
        ASC <- c(rep(1, (nalternatives - 1) * nquestions * nblocks),
                 rep(0, nquestions * nblocks))
      }
    }
    # ASC for labeled design
    else {
      # with opt-out option
      if (optout == TRUE) {
        ASC <- kronecker(diag(1, nalternatives), rep(1, nquestions * nblocks))
        if (nalternatives >= 2) {
            colnames(ASC) <- paste("ASC", 1:nalternatives, sep = "")
        }
      }
      # without opt-out option
      else {
        ASC <- rbind(kronecker(diag(1, (nalternatives - 1)),
                               rep(1, nquestions * nblocks)), 
                     matrix(0,
                            nrow = nquestions * nblocks,
                            ncol = (nalternatives - 1)))
        if (nalternatives >= 3) {
          colnames(ASC) <- paste("ASC", 1:(nalternatives - 1), sep = "")
        }
      }
    }

    # add BLOCK, QES, ALT, and ASC variables to design matrix
    my.design <- data.frame(BQS, ASC, my.design) 

    # add rows corresponding to opt-out options to design matrix
    if (optout == TRUE) {
      optout.set <- as.data.frame(matrix(c(rep(c(1:nblocks), each = nquestions),
                                           rep(c(1:nquestions), nblocks),
                                           rep((nalternatives + 1), nquestions * nblocks),
                                           rep(0, nquestions * nblocks * (ncol(my.design) - 3))),
                                         nrow = (nquestions * nblocks),
                                         ncol = ncol(my.design)))
      names(optout.set) <- names(my.design)
      my.design <- rbind(my.design, optout.set)
    }

    # format output
    my.design <- my.design[order(my.design$BLOCK, my.design$QES, my.design$ALT), ]
  }

  # binary choice model
  else {
    # common
    if (is.null(common) == FALSE) {
      my.design <- conv.ced[[1]] - conv.ced[[2]]
      colnames(my.design) <- variable.names[1:nvariables]
    }
    # optout
    else if (optout == TRUE) {
      my.design <- conv.ced[[1]]
      colnames(my.design) <- variable.names[1:nvariables]
    }
    else {
      # forced & Unlabeled
      if (unlabeled == TRUE) {
        my.design <- conv.ced[[1]] - conv.ced[[2]]
        colnames(my.design) <- variable.names[1:nvariables]
      }
      # forced & Labeled
      else {
        my.design <- cbind(conv.ced[[1]], -1 * conv.ced[[2]])
        variable.names <- paste(variable.names, 
                                rep(1:nalternatives,
                                    each = nvariables),
                                sep = "")
        colnames(my.design) <- variable.names
      }
    }
    
    # add BLOCK, QES, ALT, and ASC variables to design matrix
    BQS <- ced[[1]][, 1:3]
    ASC <- rep(1, nquestions * nblocks)
    my.design <- data.frame(BQS, ASC, my.design)
  }

  row.names(my.design) <- 1:nrow(my.design)

  return(my.design)
}

