print.partynode <- function(x, data = NULL, names = NULL,
  inner_panel = function(node) "", terminal_panel = function(node) " *",
  prefix = "", first = TRUE, digits = getOption("digits") - 2, ...)
{
  ids <- nodeids(x)
  
  if(first) {
    if(is.null(names)) names <- as.character(ids)
    cat(paste(prefix, "[", names[which(ids == id_node(x))], "] root", sep = ""))

    if(is.terminal(x)) {
      char <- terminal_panel(x)
      if(length(char) > 1L) {
        cat(paste(char[1L], "\n",
          paste(prefix, "    ", char[-1L], sep = "", collapse = "\n"),
          sep = ""), "\n")
      } else {
        cat(char, "\n")
      }
    } else {
      cat("\n")
    }
  }

  if (length(x) > 0) {
    ## add indentation
    nextprefix <- paste(prefix, "|   ", sep = "")
  
    ## split labels
    slabs <- character_split(split_node(x), data = data, digits = digits, ...)
    slabs <- ifelse(substr(slabs$levels, 1, 1) %in% c("<", ">"),
             paste(slabs$name, slabs$levels), 
             paste(slabs$name, "in", slabs$levels))

    ## kid labels
    knodes <- kids_node(x)    
    knam <- sapply(knodes, function(z) names[which(ids == id_node(z))])
    klabs <- sapply(knodes, function(z)
      if(is.terminal(z)) {
        char <- terminal_panel(z)
        if(length(char) > 1L) {
	  paste(char[1L], "\n",
	    paste(nextprefix, "    ", char[-1L], sep = "", collapse = "\n"),
	    sep = "")
	} else {
	  char
	}
      } else {
        paste(inner_panel(z), collapse = "\n")
      })

    ## merge, cat, and call recursively
    labs <- paste("|   ", prefix, "[", knam, "] ", slabs, klabs, "\n", sep = "")          
    for (i in 1:length(x)) {
      cat(labs[i])
      print.partynode(x[i], data = data, names = names[match(nodeids(x[i]), ids)],
        inner_panel = inner_panel, terminal_panel = terminal_panel,
        prefix = nextprefix,  first = FALSE, digits = digits, ...)
    }
  }
}

print.party <- function(x,
  terminal_panel = function(node) formatinfo_node(node, default = "*", prefix = ": "), tp_args = list(),
  inner_panel = function(node) "", ip_args = list(),
  header_panel = function(party) "",
  footer_panel = function(party) "",
  digits = getOption("digits") - 2, ...)
{
  ## header
  cat(paste(header_panel(x), collapse = "\n"))

  ## nodes
  if(inherits(terminal_panel, "grapcon_generator"))
    terminal_panel <- do.call("terminal_panel", c(list(x), as.list(tp_args)))
  if(inherits(inner_panel, "grapcon_generator"))
    inner_panel <- do.call("inner_panel", c(list(x), as.list(ip_args)))
  print(node_party(x), x$data, names = names(x),
    terminal_panel = terminal_panel, inner_panel = inner_panel,
    digits = digits, ...)

  ## footer
  cat(paste(footer_panel(x), collapse = "\n"))
}

print.constparty <- function(x,
  FUN = NULL, digits = getOption("digits") - 4,
  header = NULL, footer = TRUE, ...)
{
  if(is.null(FUN)) return(print(as.simpleparty(x), digits = digits,
    header = header, footer = footer, ...))

  digits <- max(c(0, digits))

  ## FIXME: terms/call/? for "ctree" objects
  if(is.null(header)) header <- !is.null(terms(x))
  header_panel <- if(header) function(party) {
    c("", "Model formula:", deparse(formula(terms(party))), "", "Fitted party:", "")
  } else function(party) ""
  
  footer_panel <- if(footer) function(party) {
    n <- width(party)
    n <- format(c(length(party) - n, n))
    
    c("", paste("Number of inner nodes:   ", n[1]),
      paste("Number of terminal nodes:", n[2]), "")
  } else function (party) ""

  y <- x$fitted[["(response)"]]
  w <- x$fitted[["(weights)"]]
  if(is.null(w)) {
    wdigits <- 0
    wsym <- "n"
  } else {
    if(isTRUE(all.equal(w, round(w)))) {
      wdigits <- 0
      wsym <- "n"
    } else {
      wdigits <- max(c(0, digits - 2))
      wsym <- "w"
    }
  }
  yclass <- class(y)[1]
  if(yclass == "ordered") yclass <- "factor"
  if(!(yclass %in% c("Surv", "factor"))) yclass <- "numeric"
  
  if(is.null(FUN)) FUN <- switch(yclass,
    "numeric" = function(y, w, digits) {
      yhat <- .pred_numeric_response(y, w)
      yerr <- sum(w * (y - yhat)^2)
      digits2 <- max(c(0, digits - 2))
      paste(format(round(yhat, digits = digits), nsmall = digits),
        " (", wsym, " = ", format(round(sum(w), digits = wdigits), nsmall = wdigits), ", err = ",
	format(round(yerr, digits = digits2), nsmall = digits2), ")", sep = "")
    },
    "Surv" = function(y, w, digits) {
      paste(format(round(.pred_Surv_response(y, w), digits = digits), nsmall = digits),
        " (", wsym, " = ", format(round(sum(w), digits = wdigits), nsmall = wdigits), ")", sep = "")
    },
    "factor" = function(y, w, digits) {
      tab <- round(.pred_factor(y, w) * sum(w))
      mc <- round(100 * (1 - max(tab)/sum(w)), digits = max(c(0, digits - 2)))
      paste(names(tab)[which.max(tab)], " (",
        wsym, " = ", format(round(sum(w), digits = wdigits), nsmall = wdigits),
        ", err = ", mc, "%)", sep = "")
    }
  )
  
  node_labs <- nodeapply(x, nodeids(x), function(node) {
    y1 <- node$fitted[["(response)"]]
    w <- node$fitted[["(weights)"]]
    if(is.null(w)) w <- rep.int(1L, NROW(y1))
    FUN(y1, w, digits)
  }, by_node = FALSE)
  node_labs <- paste(":", format(do.call("c", node_labs)))
  
  terminal_panel <- function(node) node_labs[id_node(node)]

  print.party(x, terminal_panel = terminal_panel,
    header_panel = header_panel, footer_panel = footer_panel, ...)
}
