.make_formatinfo_simpleparty <- function(x, digits = getOption("digits") - 4, sep = "")
{
  ## digit processing
  digits <- max(c(0, digits))
  digits2 <- max(c(0, digits - 2))

  ## type of predictions
  y <- node_party(x)$info$prediction
  yclass <- class(y)[1]
  if(yclass == "ordered") yclass <- "factor"
  if(!(yclass %in% c("survfit", "factor"))) yclass <- "numeric"

  ## type of weights
  n <- node_party(x)$info$n
  if(is.null(names(n))) {
    wdigits <- 0
    wsym <- "n"
  } else {
    if(names(n) == "w") {
      wdigits <- max(c(0, digits - 2))
      wsym <- "w"
    } else {
      wdigits <- 0
      wsym <- "n"
    }
  }

  ## compute terminal node labels
  FUN <- function(info) {
    yhat <- info$prediction
    if (yclass == "survfit") {
        yhat <- .median_survival_time(yhat)
        yclass <- "numeric"
    }
    if(yclass == "numeric") yhat <- format(round(yhat, digits = digits), nsmall = digits)
    w <- info$n
    yerr <- if(is.null(info$error)) "" else paste(", err = ",
      format(round(info$error, digits = digits2), nsmall = digits2),
      names(info$error), sep = "")
    rval <- paste(yhat, sep,
      " (", wsym, " = ", format(round(w, digits = wdigits), nsmall = wdigits),
      yerr, ")", sep = "")
    unlist(strsplit(rval, "\n"))
  }
  return(FUN)
}

plot.simpleparty <- function(x, digits = getOption("digits") - 4, tp_args = NULL, ...) {
  if(is.null(tp_args)) tp_args <- list(FUN = .make_formatinfo_simpleparty(x, digits = digits, sep = "\n"))
  plot.party(x, tp_args = tp_args, ...)
}

print.simpleparty <- function(x, digits = getOption("digits") - 4,
  header = NULL, footer = TRUE, ...)
{
  ## header panel
  if(is.null(header)) header <- !is.null(terms(x))
  header_panel <- if(header) function(party) {
    c("", "Model formula:", deparse(formula(terms(party))), "", "Fitted party:", "")
  } else function(party) ""
  
  ## footer panel
  footer_panel <- if(footer) function(party) {
    n <- width(party)
    n <- format(c(length(party) - n, n))
    
    c("", paste("Number of inner nodes:   ", n[1]),
      paste("Number of terminal nodes:", n[2]), "")
  } else function (party) ""

  ## terminal panel
  terminal_panel <- function(node) formatinfo_node(node,
    FUN = .make_formatinfo_simpleparty(x, digits = digits),
    default = "*", prefix = ": ")

  print.party(x, terminal_panel = terminal_panel,
    header_panel = header_panel, footer_panel = footer_panel, ...)
}

predict_party.simpleparty <- function(party, id, newdata = NULL,
    type = c("response", "prob", "node"), ...)
{
  ## get observation names: either node names or
  ## observation names from newdata
  nam <- if(is.null(newdata)) names(party)[id] else rownames(newdata)
  if(length(nam) != length(id)) nam <- NULL

  ## match type
  type <- match.arg(type)

  ## special case: fitted ids
  if(type == "node") return(structure(id, .Names = nam))

  ## predictions
  if(type == "response") {
    FUN <- function(x) x$info$prediction
  } else {
    if(is.null(node_party(party)$info$distribution)) stop("probabilities not available")
    scale <- any(node_party(party)$info$distribution > 1)
    FUN <- function(x) if(scale) prop.table(x$info$distribution) else x$info$distribution
  }
  predict_party.default(party, id, nam, FUN = FUN, ...)
}

as.simpleparty <- function(obj, ...) UseMethod("as.simpleparty")

as.simpleparty.simpleparty <- function(obj, ...) obj

as.simpleparty.party <- function(obj, ...) {
  if (is.simpleparty(obj)) {
      class(obj) <- unique(c("simpleparty", class(obj)))
      return(obj)
  }
  if (is.constparty(obj)) 
      return(as.simpleparty(as.constparty(obj)))
  stop("cannot coerce objects of class ", sQuote(class(obj)), 
       " to class ", sQuote("simpleparty"))
}

as.simpleparty.XMLNode <- function(obj, ...) as.party(obj)

as.simpleparty.constparty <- function(obj, ...) {
  ## extract and delete fitted
  fit <- obj$fitted
  obj$fitted <- NULL

  ## response info
  rtype <- class(fit[["(response)"]])[1]
  if (rtype == "ordered") rtype <- "factor"    
  if (rtype == "integer") rtype <- "numeric"

  ## extract fitted info
  FUN <- function(node, fitted) {
    fitted <- subset(fitted,
      fitted[["(fitted)"]] %in% nodeids(node, terminal = TRUE))

    if (nrow(fitted) == 0)
      return(list(prediction = NA, n = 0,
                  error = NA, distribution = NULL))

    y <- fitted[["(response)"]]
    w <- fitted[["(weights)"]]
    if(is.null(w)) {
      w <- rep(1, nrow(fitted))
      wnam <- "n"
    } else {
      wnam <- if(isTRUE(all.equal(w, round(w)))) "n" else "w" 
    }
    
    ## extract p.value (if any)
    pval <- function(node) {
      p <- info_node(node)
      if(is.list(p)) p$p.value else NULL
    }
    
    switch(rtype,
      "numeric" = {
        yhat <- .pred_numeric_response(y, w)
        list(prediction = yhat, n = structure(sum(w), .Names = wnam),
	  error = sum(w * (y - yhat)^2), distribution = NULL, p.value = pval(node))
      },
      "factor" = {
        yhat <- .pred_factor_response(y, w)
        ytab <- round(.pred_factor(y, w) * sum(w))
        list(prediction = yhat, n = structure(sum(w), .Names = wnam),
	  error = structure(sum(100 * prop.table(ytab)[names(ytab) != yhat]), .Names = "%"),
	  distribution = ytab, p.value = pval(node))
      },
      "Surv" = {
        list(prediction = .pred_Surv(y, w), n = structure(sum(w), .Names = wnam),
	  error = NULL, distribution = NULL, p.value = pval(node)) ## FIXME: change distribution format?
      })
  }

  ## convenience function for computing kid ids
  fit2id <- function(fit, idlist) {
    fit <- factor(fit)
    nlevels <- levels(fit)
    for(i in 1:length(idlist)) nlevels[match(idlist[[i]], levels(fit))] <- i
    levels(fit) <- nlevels
    ret <- factor(as.numeric(as.character(fit)), labels = 1:length(idlist), levels = 1:length(idlist))
    ret
  }

  ## cycle through node
  new_node <- function(onode, fitted) {
    if(is.terminal(onode)) return(partynode(id = onode$id,
      split = NULL, kids = NULL, surrogates = NULL,
      info = FUN(onode, fitted)))
    kids <- kids_node(onode)
    kids_tid <- lapply(kids, nodeids, terminal = TRUE)
    kids_fitted <- base::split.data.frame(fitted, fit2id(fitted[["(fitted)"]], kids_tid), drop = FALSE)
    partynode(id = onode$id, split = onode$split,
      kids = lapply(1:length(kids), function(i) new_node(kids[[i]], kids_fitted[[i]])),
      surrogates = onode$surrogates,
      info = FUN(onode, fitted))
  }
  obj$node <- new_node(node_party(obj), fit)

  class(obj) <- c("simpleparty", "party")
  return(obj)
}

is.simpleparty <- function(party) {

    chkinfo <- function(node)
        all(c("prediction", "n", "error", "distribution") %in% names(info_node(node)))

    all(nodeapply(party, ids = nodeids(party), 
                  FUN = chkinfo, by_node = TRUE))
}

