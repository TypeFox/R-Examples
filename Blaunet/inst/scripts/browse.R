
brattr <- function(h,...) {
  if ("cov" %in% ls(envir=.GlobalEnv)) {
    nw <- gwindow("Attribute File",width = 800, height = 600)
    group <- ggroup(horizontal = FALSE, cont = nw)
    vars <- gtable(cov, expand = TRUE, cont = group)
  } else gmessage("Sorry! Attribute file is not loaded.", parent = window)
}

bradj <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    nw <- gwindow("Adjacency Matrix",width = 800, height = 600)
    group <- ggroup(horizontal = FALSE, cont = nw)
    ego <- rownames(adj)
    adj1 <- cbind(ego,adj)
    vars <- gtable(adj1, expand = TRUE, cont = group)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}

brel <- function(h,...) {
  if ("el" %in% ls(envir=.GlobalEnv)) {
    nw <- gwindow("Edge List",width = 800, height = 600)
    group <- ggroup(horizontal = FALSE, cont = nw)
    vars <- gtable(el, expand = TRUE, cont = group)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}


