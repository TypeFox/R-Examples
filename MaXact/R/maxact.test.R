alterToInt <-
  function(alternative){
    alternative <- charmatch(alternative,
                             c("two.sided","greater", "less"),nomatch=0)
    if(alternative==0)
      stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
    alternative <- alternative - 1
    return(alternative)
  }

`maxact.test` <-
  function(data, max3 = TRUE, exact=TRUE, alternative = "two.sided"){
    alterInt = alterToInt(alternative)
    data.name = paste(deparse(substitute(data)), "(2x3 contingency table)")
    if(!is.matrix(data)){
      stop("'data' must be a matrix")
    } else{
      if (dim(data)[1] != 2 || dim(data)[2] != 3) 
        stop("'data' must have 2 rows and 3 columns")
      if (!is.numeric(data) || any(data < 0) || any(is.na(data))) 
        stop("all entries of 'data' must be nonnegative and finite integers")
      if (!is.integer(data)) {
        datao <- data
        data <- round(data)
        if (any(data > .Machine$integer.max)) 
          stop("'data' has entries too large to be integer")
        if (!identical(TRUE, (ax <- all.equal(datao, data)))) 
          warning("'data' has been rounded to integer: ", 
                  ax)
        storage.mode(data) <- "integer"
      }
    }  

    method = ""
    statistic = double(1)
    names(statistic) = ""
    pval = double(1)
    if(exact){
      if(max3) {
        method = "Exact MAX3 test"
        names(statistic) = "max3"
      }else {
        method = "Exact MAX2 test"
        names(statistic) = "max2"
      }
      cdata = c(t(data))
      tmp = .C(c_maxact_test, as.integer(cdata), as.integer(max3), as.integer(alterInt), statistic=statistic, pval=pval)
      statistic = tmp$statistic
      pval = tmp$pval
    }else{
      cdata = c(t(data))
      statistic = .C(c_maxact, as.integer(cdata), as.integer(max3), as.integer(alterInt), statistic=statistic)$statistic
      if (max3) {
        method = "MAX3 test (normal approximation)"
        names(statistic) = "max3"
        pval = fnorm3(statistic, data, alterInt)
      }else{
        method = "MAX2 test (normal approximation)"
        names(statistic) = "max2"
        pval = fnorm2(statistic, data, alterInt)
      }
    }

    result = list(statistic = statistic, data.name = data.name, 
      p.value = pval, method = method, data = data, alternative=alternative)
    attr(result, "class") <- "htest"
    return(result)
  }

`catt.test` <-
  function(data, theta, exact=TRUE, alternative = "two.sided"){
    alterInt = alterToInt(alternative)
    data.name = paste(deparse(substitute(data)), "(2x3 contingency table)")
    if(!is.matrix(data))
      stop("'data' must be a matrix")
    if(is.matrix(data)){
      if (dim(data)[1] != 2 || dim(data)[2] != 3) 
        stop("'data' must have 2 rows and 3 columns")
      if (!is.numeric(data) || any(data < 0) || any(is.na(data))) 
        stop("all entries of 'data' must be nonnegative and finite integers")
      if (!is.integer(data)) {
        datao <- data
        data <- round(data)
        if (any(data > .Machine$integer.max)) 
          stop("'data' has entries too large to be integer")
        if (!identical(TRUE, (ax <- all.equal(datao, data)))) 
          warning("'data' has been rounded to integer: ", 
                  ax)
        storage.mode(data) <- "integer"
      }
    }  


    method = ""
    statistic = double(1)
    names(statistic) = ""
    pval = double(1)

    method = paste("Cochran-Armitage trend test with theta", deparse(theta))
    cdata <- c(t(data))
    
    if(exact){
      method = paste("Exact", method)
      tmp = .C(c_catt_test, as.integer(cdata), as.double(theta), as.integer(alterInt), statistic=statistic, pval=pval)
      statistic = tmp$statistic
      pval = tmp$pval
    }else{
      method <- paste(method, "(normal approximation)")
      statistic = .C(c_catt, as.integer(cdata), as.double(theta), as.integer(alterInt), statistic=statistic)$statistic
      pval =
        switch(alterInt+1,
               pnorm(abs(statistic), lower.tail=F)*2,
               pnorm(statistic, lower.tail=F),
               pnorm(statistic))
    }

    names(statistic) = paste("CATT(", deparse(theta), ")", sep = "")


    result = list(statistic = statistic, data.name = data.name, 
      p.value = pval, method = method, data = data, alternative=alternative)
    attr(result, "class") <- "htest"
    return(result)
  }

