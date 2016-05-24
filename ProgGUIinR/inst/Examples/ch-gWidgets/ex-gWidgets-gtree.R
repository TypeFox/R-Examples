### R code from vignette source '~/GUI/ch-gWidgets/ex-gWidgets-gtree.Rnw'

###################################################
### code chunk number 1: ex-gWidgets-gtree.Rnw:1-2
###################################################
require(gWidgets)


###################################################
### code chunk number 2: party
###################################################
require(party)
data("GlaucomaM", package = "ipred")      # load data
gt <- ctree(Class ~ ., data = GlaucomaM)  # fit model


###################################################
### code chunk number 3: offspring
###################################################
offspring <- function(key, offspring.data) {
  if(missing(key) || length(key) == 0)  # which party node?
    node <- 1
  else
    node <- as.numeric(key[length(key)]) # key is a vector

  if(nodes(gt, node)[[1]]$terminal)    # return if terminal
    return(data.frame(node = node, hasOffspring = FALSE,
                      description = "terminal",
                      stringsAsFactors = FALSE))

  DF <- data.frame(node = integer(2), hasOffspring = logical(2),
                   description = character(2), 
                   stringsAsFactors = FALSE)
  ## party internals
  children <-  c("left","right")
  ineq <- c(" <= "," >  ")
  varName <- nodes(gt, node)[[1]]$psplit$variableName
  splitPoint <- nodes(gt, node)[[1]]$psplit$splitpoint

  for(i in 1:2) {
    DF[i,1] <- nodes(gt, node)[[1]][[children[i]]][[1]]
    DF[i,2] <- !nodes(gt, DF[i,1])[[1]]$terminal
    DF[i,3] <- paste(varName, splitPoint, sep = ineq[i])
  }
  DF                                    # returns a data frame
}


###################################################
### code chunk number 4: makeGUI
###################################################
window <- gwindow("Example of gtree")
group <- ggroup(cont = window, horizontal = FALSE)
label <- glabel("Click on the tree to investigate the partition", 
            cont = group)
tree <- gtree(offspring, cont = group, expand = TRUE)


###################################################
### code chunk number 5: eventHandler
###################################################
addHandlerDoubleclick(tree, handler = function(h,...) {
  node <- as.numeric(svalue(h$obj))
  if(nodes(gt, node)[[1]]$terminal) {   # if terminal plot
    weights <- as.logical(nodes(gt,node)[[1]]$weights)
    plot(response(gt)[weights, ])
  }})


