#----------------------------------------------------------------------------------------
# the problem with this method is that the panel.fun do not know the correct 
# subject factor because it had not subdivide in therms of the | in the formula 
# so the we have to include coplotnew() to rectify this 
par.plot<- function( formula = NULL,
                        data = NULL,
                    subjects = NULL,
                       color = TRUE,
                  show.given = TRUE,
                 ... 
                    )
{#1
#----------------------------------------------------------------------
# this defines an ammended version of the coplot with subjects as a new argument
   coplotnew <- function (formula, 
                        data, 
                        subjects, 
                        given.values, 
                        panel = points, 
                        rows, 
                        columns, 
                        show.given = TRUE, 
                        col = par("fg"), 
                        pch = par("pch"), 
                        bar.bg = c(num = gray(0.8), fac = gray(0.95)), 
                        xlab = c(x.name, paste("Given :", a.name)), 
                        ylab = c(y.name, paste("Given :",  b.name)), 
                        subscripts = FALSE, 
                        axlabels = function(f) 
                        abbreviate(levels(f)), 
                        number = 6, 
                        overlap = 0.5, 
                        xlim,  
                        ylim, 
                        ...) 
  {#2
  #--------------------------------------------------------------------
    deparen <- function(expr) 
    {
        while (is.language(expr) && !is.name(expr) && deparse(expr[[1]]) == 
            "(") expr <- expr[[2]]
        expr
    }
  #--------------------------------------------------------------------
    bad.formula <- function() stop("invalid conditioning formula")
    bad.lengths <- function() stop("incompatible variable lengths")
    formula <- deparen(formula)
    if (!inherits(formula, "formula")) 
        bad.formula()
    y <- deparen(formula[[2]])
    rhs <- deparen(formula[[3]])
    if (deparse(rhs[[1]]) != "|") 
        bad.formula()
    x <- deparen(rhs[[2]])
    rhs <- deparen(rhs[[3]])
    if (is.language(rhs) && !is.name(rhs) && (deparse(rhs[[1]]) == 
        "*" || deparse(rhs[[1]]) == "+")) 
        {
        have.b <- TRUE
        a <- deparen(rhs[[2]])
        b <- deparen(rhs[[3]])
        }
    else 
        {
        have.b <- FALSE
        a <- rhs
        }
    if (missing(data)) 
        data <- parent.frame()
    x.name <- deparse(x)
    x <- eval(x, data, parent.frame())
    nobs <- length(x)
    y.name <- deparse(y)
    y <- eval(y, data, parent.frame())
    if (length(y) != nobs) 
        bad.lengths()
    a.name <- deparse(a)
    a <- eval(a, data, parent.frame())
    if (length(a) != nobs) 
        bad.lengths()
    if (is.character(a)) 
        a <- as.factor(a)
    a.is.fac <- is.factor(a)
    if (have.b) 
       {
        b.name <- deparse(b)
        b <- eval(b, data, parent.frame())
        if (length(b) != nobs) 
            bad.lengths()
        if (is.character(b)) 
            b <- as.factor(b)
        b.is.fac <- is.factor(b)
        missingrows <- which(is.na(x) | is.na(y) | is.na(a) | 
            is.na(b))
       }
    else 
       {
        missingrows <- which(is.na(x) | is.na(y) | is.na(a))
        b <- NULL
        b.name <- ""
       }
    number <- as.integer(number)
    if (length(number) == 0 || any(number < 1)) 
        stop("number must be integer >= 1")
    if (any(overlap >= 1)) 
        stop("overlap must be < 1 (and typically >= 0).")
    bad.givens <- function() stop("invalid given.values")
    if (missing(given.values)) 
       {
        a.intervals <- if (a.is.fac) 
            {
            i <- seq(along = a.levels <- levels(a))
            a <- as.numeric(a)
            cbind(i - 0.5, i + 0.5)
            }
        else co.intervals(a, number = number[1], overlap = overlap[1])
        b.intervals <- if (have.b) 
            {
            if (b.is.fac) 
                {
                i <- seq(along = b.levels <- levels(b))
                b <- as.numeric(b)
                cbind(i - 0.5, i + 0.5)
                }
            else 
                {
                if (length(number) == 1) 
                  number <- rep.int(number, 2)
                if (length(overlap) == 1) 
                  overlap <- rep.int(overlap, 2)
                co.intervals(b, number = number[2], overlap = overlap[2])
                }
            }
        }
    else 
        {
        if (!is.list(given.values)) 
            given.values <- list(given.values)
        if (length(given.values) != (if (have.b) 
            2
        else 1)) 
            bad.givens()
        a.intervals <- given.values[[1]]
        if (a.is.fac) 
                {
            a.levels <- levels(a)
            if (is.character(a.intervals)) 
                a.intervals <- match(a.intervals, a.levels)
            a.intervals <- cbind(a.intervals - 0.5, a.intervals + 
                0.5)
            a <- as.numeric(a)
                }
        else if (is.numeric(a)) 
                {
            if (!is.numeric(a.intervals)) 
                bad.givens()
            if (!is.matrix(a.intervals) || ncol(a.intervals) != 
                2) 
                a.intervals <- cbind(a.intervals - 0.5, a.intervals + 
                  0.5)
                }
        if (have.b) 
                {
            b.intervals <- given.values[[2]]
            if (b.is.fac) 
                 {
                b.levels <- levels(b)
                if (is.character(b.intervals)) 
                  b.intervals <- match(b.intervals, b.levels)
                b.intervals <- cbind(b.intervals - 0.5, b.intervals + 
                  0.5)
                b <- as.numeric(b)
                 }
            else if (is.numeric(b)) 
                 {
                if (!is.numeric(b.intervals)) 
                  bad.givens()
                if (!is.matrix(b.intervals) || ncol(b.intervals) != 
                  2) 
                  b.intervals <- cbind(b.intervals - 0.5, b.intervals + 
                    0.5)
                 }
               }
        }
    if (any(is.na(a.intervals)) || (have.b && any(is.na(b.intervals)))) 
        bad.givens()
    if (have.b) 
       {
        rows <- nrow(b.intervals)
        columns <- nrow(a.intervals)
        nplots <- rows * columns
        if (length(show.given) < 2) 
            show.given <- rep.int(show.given, 2)
       }
    else 
       {
        nplots <- nrow(a.intervals)
        if (missing(rows)) 
           {
            if (missing(columns)) 
              {
                rows <- ceiling(round(sqrt(nplots)))
                columns <- ceiling(nplots/rows)
              }
            else rows <- ceiling(nplots/columns)
           }
        else if (missing(columns)) 
            columns <- ceiling(nplots/rows)
        if (rows * columns < nplots) 
            stop("rows * columns too small")
        }
    total.columns <- columns
    total.rows <- rows
    f.col <- f.row <- 1
    if (show.given[1]) 
    {
        total.rows <- rows + 1
        f.row <- rows/total.rows
    }
    if (have.b && show.given[2]) 
    {
        total.columns <- columns + 1
        f.col <- columns/total.columns
    }
    mar <- if (have.b) 
        rep.int(0, 4)
    else c(0.5, 0, 0.5, 0)
    oma <- c(5, 6, 5, 4)
    if (have.b) 
    {
        oma[2] <- 5
        if (!b.is.fac) 
            oma[4] <- 5
    }
    if (a.is.fac && show.given[1]) 
        oma[3] <- oma[3] - 1
    opar <- par(mfrow = c(total.rows, total.columns), oma = oma, 
        mar = mar, xaxs = "r", yaxs = "r", new = FALSE)
    on.exit(par(opar))
    plot.new()
    if (missing(xlim)) 
        xlim <- range(as.numeric(x), finite = TRUE)
    if (missing(ylim)) 
        ylim <- range(as.numeric(y), finite = TRUE)
    pch <- rep(pch, length = nobs)
    col <- rep(col, length = nobs)
    #-----------------------------------
    do.panel <- function(index, subscripts = FALSE, id) 
    {
    #---------
        Paxis <- function(side, x) 
         {
            if (nlevels(x)) 
            {
                lab <- axlabels(x)
                axis(side, labels = lab, at = seq(lab), xpd = NA)
            }
            else axis(side, xpd = NA)
         }
     #--------
        istart <- (total.rows - rows) + 1
        i <- total.rows - ((index - 1)%/%columns)
        j <- (index - 1)%%columns + 1
        par(mfg = c(i, j, total.rows, total.columns))
        plot.new()
        plot.window(xlim, ylim)
        if (any(is.na(id))) 
            id[is.na(id)] <- FALSE
        if (any(id)) 
        {
            grid(lty = "solid")
            if (subscripts) 
                panel(x[id], y[id], subjects[id], subscripts = id, col = col[id], 
                  pch = pch[id], ...)
            else panel(x[id], y[id], subjects[id], col = col[id], pch = pch[id], 
                ...)
        }
        if ((i == total.rows) && (j%%2 == 0)) 
            Paxis(1, x)
        else if ((i == istart || index + columns > nplots) && 
            (j%%2 == 1)) 
            Paxis(3, x)
        if ((j == 1) && ((total.rows - i)%%2 == 0)) 
            Paxis(2, y)
        else if ((j == columns || index == nplots) && ((total.rows - 
            i)%%2 == 1)) 
            Paxis(4, y)
        box()
       } #end of do.panel
    #---------------------------
    if (have.b) 
    {
        count <- 1
        for (i in 1:rows) 
        {
            for (j in 1:columns) 
            {
                id <- ((a.intervals[j, 1] <= a) & (a <= a.intervals[j, 
                  2]) & (b.intervals[i, 1] <= b) & (b <= b.intervals[i, 
                  2]))
                do.panel(count, subscripts, id)
                count <- count + 1
            }
        }
     }
    else 
    {
        for (i in 1:nplots) 
        {
            id <- ((a.intervals[i, 1] <= a) & (a <= a.intervals[i, 
                2]))
            do.panel(i, subscripts, id)
        }
    }
    mtext(xlab[1], side = 1, at = 0.5 * f.col, outer = TRUE, 
        line = 3.5, xpd = NA)
    mtext(ylab[1], side = 2, at = 0.5 * f.row, outer = TRUE, 
        line = 3.5, xpd = NA)
    if (length(xlab) == 1) 
        xlab <- c(xlab, paste("Given :", a.name))
    if (show.given[1]) 
    {
        par(fig = c(0, f.col, f.row, 1), mar = mar + c(3 + (!a.is.fac), 
            0, 0, 0), new = TRUE)
        plot.new()
        nint <- nrow(a.intervals)
        a.range <- range(a.intervals, finite = TRUE)
        plot.window(a.range + c(0.03, -0.03) * diff(a.range), 
            0.5 + c(0, nint))
        rect(a.intervals[, 1], 1:nint - 0.3, a.intervals[, 2], 
            1:nint + 0.3, col = bar.bg[if (a.is.fac) 
                "fac"
            else "num"])
        if (a.is.fac) 
          {
            text(apply(a.intervals, 1, mean), 1:nint, a.levels)
          }
        else 
          {
            axis(3, xpd = NA)
            axis(1, labels = FALSE)
          }
        box()
        mtext(xlab[2], 3, line = 3 - a.is.fac, at = mean(par("usr")[1:2]), 
            xpd = NA)
     }
     else 
     {
        mtext(xlab[2], 3, line = 3.25, outer = TRUE, at = 0.5 * 
            f.col, xpd = NA)
     }
     if (have.b) 
     {
        if (length(ylab) == 1) 
            ylab <- c(ylab, paste("Given :", b.name))
        if (show.given[2]) 
         {
            par(fig = c(f.col, 1, 0, f.row), mar = mar + c(0, 
                3 + (!b.is.fac), 0, 0), new = TRUE)
            plot.new()
            nint <- nrow(b.intervals)
            b.range <- range(b.intervals, finite = TRUE)
            plot.window(0.5 + c(0, nint), b.range + c(0.03, -0.03) * 
                diff(b.range))
            rect(1:nint - 0.3, b.intervals[, 1], 1:nint + 0.3, 
                b.intervals[, 2], col = bar.bg[if (b.is.fac) 
                  "fac"
                else "num"])
            if (b.is.fac) 
             {
                text(1:nint, apply(b.intervals, 1, mean), b.levels, 
                  srt = 90)
             }
            else 
             {
                axis(4, xpd = NA)
                axis(2, labels = FALSE)
             }
            box()
            mtext(ylab[2], 4, line = 3 - b.is.fac, at = mean(par("usr")[3:4]), 
                xpd = NA)
          }
        else 
          {
            mtext(ylab[2], 4, line = 3.25, at = 0.5 * f.row, 
                outer = TRUE, xpd = NA)
          }
      }
    if (length(missingrows) > 0) 
      {
        cat("\nMissing rows:", missingrows, "\n")
        invisible(missingrows)
      }
}
# end of 
# here is the panel function 
#----------------------------------------------------------------------
  panel.fun <- function ( x, y, subjects, 
                           subscripts = id,
                         col = par("col"), 
                          bg = NA, 
                         pch = par("pch"), 
                         cex = par("cex"), 
                           ...) 
   {
     col <- 1
     lll <- nlevels(subjects)
     for (i in 1:lll)
       {
         sy <-  y[subjects==i]
         sx <-  x[subjects==i]
        ssy <- subset(sy,!is.na(sy))
        ssx <- subset(sx,!is.na(sy))
        ssy <- ssy[order(ssx)]
        ssx <- ssx[order(ssx)]
        lines(ssx,ssy,type="o",col=col)
        if (color==TRUE)  col<-col+1
      }
   }
# end of panel.fun
#----------------------------------------------------------------------  
#we need this to interpreted the formula
# rcParseFormula <- function (model) 
#   {
#     parseCond <- function(model) 
#     {
#        model <- eval(parse(text = paste("~", deparse(model))))[[2]]
#        model.vars <- list()
#        while (length(model) == 3 && (model[[1]] == as.name("*") || 
#            model[[1]] == as.name("+"))) {
#            model.vars <- c(model.vars, model[[3]])
#            model <- model[[2]]
#        }
#        rev(c(model.vars, model))
#     } # end of parseCond
#     if (!inherits(model, "formula")) 
#        stop("model must be a formula object")
#     ans <- list(left = NULL, right = NULL, condition = NULL, left.name = character(0), 
#            right.name = character(0))
#         if (length(model) == 3) 
#            {
#             ans$left <- eval(model[[2]])
#             if (inherits(ans$left, "POSIXt")) 
#             ans$left <- as.POSIXct(ans$left)
#             ans$left.name <- deparse(model[[2]])
#             }       
#     model <- model[[length(model)]]
#     if (length(model) == 3 && model[[1]] == as.name("|")) 
#     {
#        model.vars <- parseCond(model[[3]])
#        ans$condition <- vector("list", length(model.vars))
#        names(ans$condition) <- sapply(model.vars, deparse)
#        for (i in seq(along = model.vars)) 
#          {
#            ans$condition[[i]] <- eval(model.vars[[i]],sys.frame(sys.parent(2))) #MSWednesday, April #16, 2003 at 09:58
#            if (inherits(ans$condition[[i]], "POSIXt")) 
#                ans$condition[[i]] <- as.POSIXct(ans$condition[[i]])
#          }
#        model <- model[[2]]
#     }
#     ans$right <- eval(model,sys.frame(sys.parent(2))) #MSWednesday, April 16, 2003 at 09:58
#      if (inherits(ans$right, "POSIXt"))  ans$right <- as.POSIXct(ans$right)
#       ans$right.name <- deparse(model)
#    ans
#   } 
# end of rcParseFormula 
#----------------------------------------------------------------------
# here is the main function    
#----------------------------------------------------------------------
#if (!is.null(data)) {attach(data) ; on.exit(detach(data))} 
if (is.null(data)) stop("The argument data is required")
         subjects <- eval(substitute(subjects), envir=as.environment(data))
            form <- as.formula(formula, env=as.environment(data))
     subjectExist <- length(form[[3]])>1
if (is.null(subjects)) stop("the subjects factor is not set")  
if (!is.factor(subjects)) stop("the subjects argument should be a factor")
  subjects <- factor(as.integer(subjects))
#       form <- rcParseFormula(formula)
if (!subjectExist)
    {    
  id <- NULL # added Thursday, March 27, 2008 at 16:54
          y <- with(data, eval(form[[2]]))
          x <- with(data, eval(form[[3]]))
    # xlabel <- form$right.name 
    # ylabel <- form$left.name
       ymax <- max(y)
       ymin <- min(y)
       xmax <- max(x)
        xmin <- min(x)
         col <- 1
         lll <- nlevels(subjects)
         plot(form, type="n", data=data)
         
         for (i in 1:lll)
           {
            sy <-  y[subjects==i]
            sx <-  x[subjects==i]
           ssy <- subset(sy,!is.na(sy))
           ssx <- subset(sx,!is.na(sy))
           lines(ssx,ssy,type="o",col=col)
           if (color==TRUE)  col<-col+1
           }
    }  
else 
    {  
      coplotnew(formula,
                data = data, 
                subjects =subjects,    
               panel = panel.fun, 
          show.given = show.given,
                  bg = "wheat", 
              bar.bg = c(num="light blue",  fac ="light blue" ),
              ... 
               ) 
    }        
}          
        
          
