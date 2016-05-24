hotelling.stat = function(x, y, shrinkage = FALSE)
{
    ## get the sample sizes for each sample
    nx = nrow(x)
    ny = nrow(y)

    px = ncol(x)
    py = ncol(y)

    if(px!=py)
        stop("Both samples must have the same number of variables (columns)")

    p = px

    if(nx + nx < p - 1 & !shrinkage)
        stop("The sample sizes (nx + ny) must be 1 greater than the number of columns")

    mx = apply(x, 2, mean)
    my = apply(y, 2, mean)

    sx = NULL
    sy = NULL

    if(!shrinkage){
        sx = cov(x)
        sy = cov(y)
    }else{
        sx = cov.shrink(x, verbose = FALSE)
        sy = cov.shrink(y, verbose = FALSE)
    }

    sPooled = ((nx - 1)*sx + (ny - 1)*sy)/(nx + ny - 2)
    sPooledInv = solve(sPooled)

    T2 = t(mx - my) %*% sPooledInv %*% (mx - my)*nx*ny/(nx + ny)
    m = (nx + ny - p - 1)/(p*(nx + ny - 2))

    invisible(list(statistic = as.vector(T2), m = m, df = c(p, nx + ny - p - 1),
                   nx = nx, ny = ny, p = p))
}

hotelling.test = function(x, ...){
    UseMethod("hotelling.test")
}

hotelling.test.default = function(x, y, shrinkage = FALSE, perm = FALSE,
                                  B = 10000, progBar = (perm && TRUE), ...){
    if(!perm){
        stats = hotelling.stat(x, y, shrinkage)
        pVal = with(stats, 1 - pf(m*statistic, df[1], df[2]))
        output = list(stats = stats, pval = pVal)
        class(output) = "hotelling.test"
        invisible(output)
    }else{
        stats = hotelling.stat(x, y, shrinkage)
        res = rep(0, B)

        nx = stats$nx
        ny = stats$ny
        N = nx + ny
        T0 = stats$statistic

        idx = 1:N
        X = rbind(x, y)
        
        onePercent = floor(B / 100)
        pb = NULL
        if(progBar){
          pb = txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
        }
        j = 0
        k = 0

        for(i in 1:B){
            i1 = sample(idx, nx)
            x1 = X[i1,]
            x2 = X[-i1,]

            res[i] = hotelling.stat(x1, x2, shrinkage)$statistic
            j = j + 1
            if(j == onePercent && progBar){
              k = k + 1
              j = 0
              setTxtProgressBar(pb, k)
            }
        }

        pVal = sum(res > T0)/B
        output = list(stats = stats, pval = pVal , results = res)
        class(output) = "hotelling.test"
        invisible(output)
    }
}

hotelling.test.formula = function(x, data = NULL, pair = c(1,2), ...){
    if(missing(x) || class(x) != "formula")
        stop("missing or incorrect formula")

    form = x
    form[[3]] = x[[2]]
    form[[2]] = x[[3]]

    mf = model.frame(form, data)

    group = model.response(mf)
    variables = mf[,-1]

    split.data = split(variables,group)


    x1 = as.matrix(split.data[[pair[1]]])
    x2 = as.matrix(split.data[[pair[2]]])

    hotelling.test(x1, x2, ...)
}

plot.hotelling.test = function(x,...){
    if(is.na(match("results",names(x))))
        stop("Plotting only works if you have used the permutation test")

    dotArgNames = names(list(...))
    
    if("xlim" %in% dotArgNames){
      with(x,{
         h = hist(stats$m*results, main = "Distribution of permuted test stats",
                  xlab = expression(T^2),...);
         T0 = with(stats, m*statistic)
         lines(rep(T0, 2), c(0, max(h$counts)), lwd = 2)
         }
      )
    }else{
      ## Try to do something sensible to make sure the test statistic is visible
      ## in the plotting window, ONLY if the user isn't over-riding xlim
      with(x,{
        r = range(c(stats$m * results, stats$m * stats$statistic));
        pr = pretty(r);
        h = hist(stats$m * results, main = "Distribution of permuted test stats",
             xlab = expression(T^2), xlim = c(min(pr), max(pr)), ...);
        #abline(v = with(stats, m*statistic, lwd = 2))
        T0 = with(stats, m*statistic)
        lines(rep(T0, 2), c(0, max(h$counts)), lwd = 2)
        }
      )
    }
}

print.hotelling.test = function(x, ...){
    if(is.na(match("results",names(x)))){
        with(x,{
            with(stats,{
                cat(paste("Test stat: ", signif(m*statistic, 5), "\n"));
                cat(paste("Numerator df: ", stats$df[1], "\n"))});
            cat(paste("Denominator df: ", x$stats$df[2], "\n"));
            cat(paste("P-value: ", signif(x$pval, 4), "\n"))})
    }else{
        with(x,{
            with(stats,{
                cat(paste("Test stat: ", signif(m*statistic, 5), "\n"));
                cat(paste("Numerator df: ", stats$df[1], "\n"))});
            cat(paste("Denominator df: ", x$stats$df[2], "\n"));
            cat(paste("Permutation P-value: ", signif(x$pval, 4), "\n"));
            cat(paste("Number of permutations :", length(x$results), "\n"))          })
    }
}



hotel.stat = function(x, y, shrinkage = FALSE){
    hotelling.stat(x, y, shrinkage)
}

hotel.test = function(x, ...){
    hotelling.test(x, ...)
}
