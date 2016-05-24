##' \code{table} with addmargins and test
##'
##' @param x x
##' @param y y
##' @param useNA useNA
##' @param margin margin
##' @param addmargins addmargins
##' @author David Hajage
##' @keywords internal
n.table <- function (x, y = NULL, useNA = c("no", "ifany", "always"), margin = 0:2, addmargins = FALSE, test = FALSE, test.tabular = test.tabular.auto, show.test = display.test, plim = 4, show.method = TRUE, na = FALSE) {
  if (is.null(y)) {
    n <- table(x, useNA = useNA)
  } else {
    n <- table(x, y, useNA = useNA)
  }
  if (addmargins) {
    if (0 %in% margin)
      margin <- 1:2
    if (is.null(y))
      margin <- 2
    totmargin <- ifelse(margin == 1, 2, margin)
    totmargin <- ifelse(margin == 2, 1, totmargin)
    
    if (length(totmargin) != 0)
      n <- addmargins(n, margin = totmargin, FUN = list(Total = sum), quiet = TRUE)
  }
  if (test)
    attr(n, "test") <- show.test(test.tabular(x, y, na), digits = plim, method = show.method)
  n
}

##' \code{prop.table} with cell proportion, addmargins and propNA
##'
##' @param ... ...
##' @param useNA useNA
##' @param margin margin
##' @param addmargins addmargins
##' @param propNA propNA
##' @author David Hajage
##' @keywords internal
p.table <- function (x, y = NULL, useNA = c("no", "ifany", "always"), margin = 0:2, addmargins = FALSE, propNA = TRUE) {
  if (is.null(y)) {
    n <- table(x, useNA = useNA)
  } else {
    n <- table(x, y, useNA = useNA)
  }
  nn <- n
  if (!propNA) {
    if (length(dim(nn)) == 2) {
      nn[is.na(rownames(nn)), ] <- 0
      nn[, is.na(colnames(nn))] <- 0
    } else {
      nn[is.na(names(nn))] <- 0
    }
  }
  if (is.null(y))
    margin <- margin[margin != 2]
  props <- lapply(margin, function(margin) {
    if (margin != 0) {
      prop <- sweep(nn, margin, margin.table(nn, margin), "/", check.margin = FALSE)
    } else {
      prop <- nn/sum(nn)
    }
    prop
  })
  names(props) <- sapply(as.character(margin), function(x) switch(x, "0" = "cell", "1" = "row", "2" = "col"))

  if (addmargins) {
    if (1 %in% margin) {
      props$row <- as.table(cbind(props$row, Total = margin.table(nn, 1)/sum(margin.table(nn, 1))))
      if (2 %in% margin) {
        props$row <- as.table(rbind(props$row, Total = NA))
      }
    } 
    if (2 %in% margin) {
      props$col <- as.table(rbind(props$col, Total = margin.table(nn, 2)/sum(margin.table(nn, 2))))
      if (1 %in% margin) {
        props$col <- as.table(cbind(props$col, Total = NA))
      }
    }
    if (0 %in% margin) {
      props$cell <- as.table(cbind(props$cell, Total = NA))
      props$cell <- as.table(rbind(props$cell, Total = NA))
    }
  }
  props <- lapply(props, function(prop) {
    if (!propNA) {
      if (length(dim(prop)) == 2) {
        prop[is.na(rownames(prop)), ] <- NA
        prop[, is.na(colnames(prop))] <- NA
      } else {
        prop[is.na(names(prop))] <- NA
      }
    }
    prop
  })
  props
}
