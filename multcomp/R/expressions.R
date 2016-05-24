
# $Id: expressions.R 435 2016-05-04 01:46:55Z sks $

### determine if an expression `x' can be interpreted as numeric
is_num <- function(x) {
    if (length(x) == 1) return(is.numeric(x))
    if (length(x) == 2) return(is.name(x[[1]]) && is.numeric(x[[2]]))
    return(FALSE)
}

### expressions
as.char <- function(ex) {
    if (length(ex) == 1) return(as.character(ex))
    if (length(ex) == 3 && ex[[1]] == ":")
        return(paste(as.char(ex[[2]]), ":", 
               as.char(ex[[3]]), sep = ""))
    stop("Failed to convert expression ", ex, " to character")
}

### extract coefficients and variable names
coefs <- function(ex) {

    ### `a'
    if (length(ex) == 1 && !is_num(ex))
        return(list(coef = 1, var = as.char(ex)))

    ### `-a'
    if (length(ex) == 2 && (ex[[1]] == "-" && !is_num(ex[[2]])))
        return(list(coef = -1, var = as.char(ex[[2]])))

    if (length(ex) == 3) {

        ### x:y
        if (ex[[1]] == ":")
            return(list(coef = 1, 
                        var = as.char(ex)))

        ### `2 * a'
        if (ex[[1]] == "*" && (is_num(ex[[2]]) && !is_num(ex[[3]])))
            return(list(coef = eval(ex[[2]]), var = as.char(ex[[3]])))

        cf <- coefs(ex[[3]])
        if (ex[[1]] == "-") 
            cf$coef <- cf$coef * (-1)

        return(cf)
    }
    stop("cannot interpret expression ", sQuote(ex), " as linear function")
}

### extract left hand side of an expression
lhs <- function(ex) {

    if (length(ex) != 1)
        stop("expression is not of length 1")

    if (length(ex[[1]]) != 3)
        stop("expression ", sQuote(ex), 
             " does not contain a left and right hand side")

    return(ex[[1]][[2]])
}

### extract right hand side of an expression
rhs <- function(ex) {

    if (length(ex) != 1)
        stop("expression is not of length 1")

    if (length(ex[[1]][[3]]) == 2)
        return(-ex[[1]][[3]][[2]])

    rhs <- ex[[1]][[3]]
    if (!is_num(rhs) || length(rhs) > 1)
        stop("right hand side of expression ", sQuote(ex), 
             " is not a scalar numeric")
    return(rhs)
}

### extract direction of the _alternative_
side <- function(ex) {

    side <- as.char(ex[[1]][[1]])
    if (!(side %in% c("<=", ">=", "==", "=")))
        stop("does not contain ", sQuote("<=, >=, =="))
    alternative <- switch(side, 
        "<=" = "greater",
        ">=" = "less",
        "==" = "two.sided",
        "=" = "two.sided")
    return(alternative)
}


expression2coef <- function(ex, vars, debug = F) {

   ### uses walkCode and makeCodeWalker from codetools

   m.rhs <- rhs(ex)
   m.lhs <- lhs(ex)

   attr( m.lhs, 'coef') <- 1

   if ( debug ) {
        message('expression2coef',': lhs is ', sQuote(paste0(deparse(m.lhs),collapse='')))
        message('expression2coef',': rhs is ', sQuote(paste0(deparse(m.rhs),collapse='')))
   }


   effects <-
   walkCode(m.lhs,
            makeCodeWalker( # dispatch operators
                            handler = function(v, w) {
                                      if ( debug ) w$trace('handler',v,w)

                                      switch( v,
                                             '-'    = w$sub,
                                             '+'    = w$add,
                                             '*'    = w$mul,
                                             '/'    = w$div,
                                             '('    = w$exp,
                                             ':'    = w$ita,
                                 #            w$fatal('handler','operator ', sQuote(v), ' is not supported')
                                              w$eval
                                             )
                            },

                            is.effect<-function(x) {
                                # vars is visible as a parameter to the enclosing function expression2coef
                                as.character(x) %in% vars 
                            },

                            eval = function(v,w) {
                                   if ( debug ) w$trace('eval',v,w)

                                   parms <- c()
                                   for ( e in as.list(v)[-1] )  {
                                         coef <- 1
                                         if ( debug )
                                              message('eval',': walking ', sQuote(e), ' with coef = ', coef)
                                         parms <- c(parms, p <- walkCode(w$setCoef(e,coef),w) )
                                         if (  is.effect( p ) )
                                               w$fatal('eval','within ', sQuote(deparse(v)) , ', the term ', sQuote(p),' ',
                                                       'must not denote an effect. Apart from that, ',
                                                       'the term must evaluate to a real valued constant')
                                   }

                                   cparms <- c()
                                   for ( e in parms ) {
                                         cparms <- c(cparms,parse(text=paste(e,'*',w$getCoef(e))))
                                   }

                                   if ( debug ) {
                                        dumped <- lapply(cparms, function(x,w) paste(x, 'with coef =', w$getCoef(x)),w)
                                        message('eval',': cparms = ', w$enum(dumped))
                                   }

                                   res <- try ( do.call(as.character(v[[1]]), as.list(cparms)), silent=T )

                                   if ( class(res) == 'try-error' )
                                        w$fatal('eval','the evaluation of the expression ', sQuote(deparse(v)),' ',
                                                       'failed with ', dQuote(attr(res,'condition')$message) )

                                   if ( length(res) != 1 || !is.numeric(res) || !is.finite(res)  )
                                        w$fatal('eval','the expression ', sQuote(deparse(v)),' ',
                                                       'did not evaluate to a real valued constant. ',
                                                       'Result is ', sQuote(res) )

                                   res <- w$setCoef(res*w$getCoef(v),1)

                                   if ( debug ) {
                                        dumped <- lapply(res, function(x,w) paste(x, 'with coef =', w$getCoef(x)),w)
                                        message('eval',': res = ', w$enum(dumped))
                                   }
                                   res
                            },

                            # call -- evaluate construct directly ( this should not be reached )
                            call = function(v,w) {
                                 if ( debug ) w$trace('call',v,w)
                                 w$fatal('call',"there is probably a syntax error within subexpression", sQuote(deparse(v)))
                            },

                            # 'a - b' or '-a' -- support for subtraction of effects or constants (but not both)
                            sub = function(v,w) {
                                  if ( debug ) w$trace('sub',v,w)

                                  uminus     <- length(as.list(v)) == 2
                                  minuend    <- as.list(v)[2]
                                  subtrahend <- as.list(v)[3]

                                  if ( debug ) {
                                       message('sub',': minuend is ', sQuote(minuend), ', coef = ', w$getCoef(minuend) )
                                       message('sub',': subtrahend is ', sQuote(subtrahend), ', coef = ', w$getCoef(subtrahend) )
                                       message('sub',': uminus is ', uminus)
                                  }


                                  exp.coef <- ifelse( uminus, -w$getCoef(v), w$getCoef(v) )
                                  res      <- c()

                                  for ( e in minuend ) {
                                        if ( debug )
                                             message('sub',': walking minuend ', sQuote(e), ', coef = ', exp.coef)
                                        res   <- c( res, walkCode(w$setCoef(e, exp.coef ),w))
                                  }

                                  if ( ! uminus ) {
                                       # if uminus becomes true, subtrahend is a list of nulls.
                                       # As a consequence, e would become null and the program would fail

                                       for ( e in subtrahend ) {
                                             if ( debug )
                                                  message('sub',': walking subtrahend ', sQuote(e), ', coef = ', -exp.coef)
                                             res   <- c( res, walkCode(w$setCoef(e, -exp.coef),w))
                                       }
                                  }

                                  sum     <- 0
                                  symbols <- c()
                                  # split result set into constants and symbols
                                  for ( e in res ) {
                                        if ( is.numeric(e) )
                                             sum <- sum + e
                                        else
                                             symbols <- c(symbols, e)
                                  }

                                  # do not allow a reference a single effect to occur multiple times.
                                  # to do: could be folded into a single effect by summing up coeffs
                                  if ( length(dups <- symbols[duplicated(symbols)]) ) {
                                       w$fatal('sub','multiple occurence of ', w$enum(dups), ' ',
                                                     'found within expression ', sQuote(deparse(v)))
                                  }


                                  # fold constants into single number
                                  if ( length(symbols) == 0 ) {
                                       return(w$setCoef(sum,1))
                                  }

                                  if ( sum  ) {
                                       w$fatal('sub','forming a difference between a constant and ',
                                                     'an effect as in ', sQuote(deparse(v)), ' ',
                                                     'is not a sensible operation')
                                  }

                                  symbols
                            },

                            # ': a b ' -- support for interaction of effects as in A:B:C:D
                            ita = function(v,w) {
                                  if ( debug ) w$trace('ita',v,w)

                                  res <- w$setCoef(as.name(deparse(v)), w$getCoef(v) )

                                  if ( debug ) {
                                       dumped <- lapply(res, function(x,w) paste(x, 'with coef =', w$getCoef(x)),w)
                                       message('ita',': res = ', w$enum(dumped))
                                  }

                                  res
                            },

                            # (expression)
                            exp = function(v,w) {
                                  if ( debug ) w$trace('exp',v,w)

                                  res <- c()
                                  for ( e in as.list(v)[-1] ) {
                                        res <- c( res, walkCode(w$setCoef(e,1),w) )
                                  }

                                  symbols <- c()
                                  for ( e in  res ) {
                                        symbols <- c( symbols, w$setCoef(e, w$getCoef(e) * w$getCoef(v) ) )
                                  }

                                  if ( debug ) {
                                       dumped <- lapply(symbols, function(x,w) paste(x, 'with coef =', w$getCoef(x)),w)
                                       message('exp',': res = ', w$enum(dumped))
                                  }
                                  symbols

                            },

                            # '+ a b' or '+ a' -- support for the addition of constants or effects (but not both)
                            add = function(v,w) {
                                  if ( debug ) w$trace('add',v,w)

                                  res <- c()
                                  for ( e in as.list(v)[-1] ) {
                                        res <- c( res, walkCode(w$setCoef(e,1),w) )
                                  }

                                  symbols <- c()
                                  sum     <- 0
                                  for ( e in res ) {
                                        if ( is.numeric(e) )
                                             sum <- sum + e
                                        else
                                             symbols <- c( symbols, e )
                                  }


                                  # fold constants into single number
                                  if ( length(symbols) == 0 ) {
                                       return( w$setCoef( sum * w$getCoef(v), 1 ) )
                                  }


                                  # do not allow to reference a single effect multiple times.
                                  # to do: could be folded into a single effect by summing up coeffs
                                  if ( length(dups <- symbols[duplicated(symbols)]) != 0 ) {
                                       w$fatal('add','multiple occurence of ', w$enum(dups),' ',
                                                     'within subexpression ', sQuote(deparse(v)))
                                  }


                                  if ( sum ) {
                                       w$fatal('add','adding up a constant and an effect ',
                                                     'as in ', sQuote(deparse(v)), ' is not a sensible operation')
                                  }

                                  # associate expression coefficient with all leafs
                                  res <- c()
                                  for ( e in symbols ) {
                                        res <- c( res, w$setCoef(e, w$getCoef(e) * w$getCoef(v) ) )
                                  }

                                  res
                            },

                            # '* a b' -- support multiplication of constants or multiplication of an effect by a constant
                            mul = function(v,w) {
                                  if ( debug ) w$trace('mul',v,w)

                                  # collect all leafs, including constants
                                  res      <- c()
                                  for ( e in as.list(v)[-1] ) {
                                        res <- c(res, walkCode(w$setCoef(e,1),w) )
                                  }

                                  # fold literals
                                  product  <- 1
                                  symbols  <- c()
                                  for ( r in res ) {
                                        if ( is.numeric(r) )
                                             product <- product * r
                                        else
                                             symbols <- c(symbols,r)
                                  }

                                  if ( product  == 0 && length(symbols) ) {
                                       w$fatal('mul','The constant part of the expression ', sQuote(deparse(v)),' ',
                                                     'evaluates to zero. This would zero out the effect(s) ', sQuote(symbols) )

                                  }

                                  # also take the expression coefficient into account
                                  product <- product * w$getCoef(v)

                                  # if only literals, return constant folding result as single number
                                  if ( length(symbols) == 0 ) {
                                       return(w$setCoef(product,1))
                                  }

                                  # prevent multiplication of fixed effects as in 'A * B' but still
                                  # allow the multiplication effects by a constant
                                  if ( length(symbols) > 1 && all(unlist(lapply(res,is.symbol)) ) ) {
                                       w$fatal('mul','the multiplication of effects ', w$enum(symbols),' ',
                                                     'as in ', sQuote(deparse(v)), ' is not a sensible operation')
                                  }

                                  # associate the folded real valued literals as a coefficient with all symbols
                                  res <- c()
                                  for ( s in symbols ) {
                                        res <- c(res, w$setCoef(s, w$getCoef(s) * product ))
                                  }

                                  if ( debug ) {
                                       dumped <- lapply(res, function(x,w) paste(x, 'with coef =', w$getCoef(x)),w)
                                       message('mul',': res = ', w$enum(dumped))
                                  }

                                  res
                            },

                            # '/ a b' --  division of an expression 'a' by a constant expression 'b'
                            div = function(v,w) {
                                  if ( debug ) w$trace('div',v,w)

                                  # const / const is allowed
                                  # fixed / const is allowed:  -> coef(fixed) <- 1/const
                                  # const / fixed is forbidden

                                  # collect all leafs, including constants

                                  if ( (lv<-length(v)) != 3 ) {
                                        w$fatal('div', 'internal error: length of language object ', sQuote(v), ' ',
                                                       'is not 3, but ', lv,'. Please file a bug report')
                                  }

                                  dividend <- c()
                                  for ( e in as.list(v)[2] ) {
                                        dividend <- c(dividend, walkCode(w$setCoef(e,1),w))
                                  }

                                  divisor  <- c()
                                  for ( e in as.list(v)[3] ) {
                                        divisor  <- c(divisor, walkCode(w$setCoef(e,1),w))
                                  }


                                  if ( length(divisor) != 1 ) {
                                       w$fatal('div', "can't divide by ", sQuote(divisor), ' in ', sQuote(deparse(v)))
                                  }

                                  if ( any(unlist(lapply(divisor,is.effect)))  ) {
                                       w$fatal('div', "cant't divide by effect ", sQuote(divisor), ' in ', sQuote(deparse(v)))
                                  }

                                  if ( any(unlist(lapply(divisor,is.symbol)))  ) {
                                       w$fatal('div', "cant't divide by symbol ", sQuote(divisor), ' in ', sQuote(deparse(v)))
                                  }

                                  divisor <- as.numeric(divisor)

                                  if ( !is.finite(divisor) || divisor == 0) {
                                       w$fatal('div', "can't divide by ", sQuote(divisor), ' in ', sQuote(deparse(v)))
                                  }

                                  res <- c()
                                  for ( s in dividend ) {
                                        if ( is.numeric(s) ) {
                                             res <- c(res, s * w$getCoef(v) / divisor )
                                        } else {
                                             res <- c(res, w$setCoef(s, w$getCoef(s) * w$getCoef(v) / divisor ) )
                                        }
                                  }

                                  if ( debug ) {
                                       message("div",": dividend = ", w$enum(dividend))
                                       message("div",": divisor = ",  w$enum(divisor))
                                       message("div",": res = ",      w$enum(res))
                                  }
                                  res
                            },

                            # leaf(e,w) -- gets called with `e` being either an effect name or a literal
                            leaf   = function(e, w) {
                                     if ( debug ) w$trace('leaf',e,w)

                                     #  leafs holding real valued constants tend to lose the coefficient
                                     #  attribute during implicit conversions. Hence, multiply the coefficient
                                     #  with the value and set the coefficient to one.
                                     if ( is.numeric(e) ) {
                                          return(w$setCoef(e * w$getCoef(e),1))
                                     }
                                     e
                            },

                            # return associated coefficient, or 1 if coefficient was not set before
                            getCoef = function(e) {
                                      a <- attr(e,'coef')
                                      ifelse( is.null(a), 1, a )
                            },

                            # set coefficient of `e` to `coef` and return e
                            setCoef = function(e,coef) {
                                      attr(e,'coef') <- coef
                                      e
                            },

                            enum = function(x) {
                                 paste0("'" ,x, "'",collapse=', ')
                            },

                            fatal = function(name,...) {
                                     stop(paste0('expression2coef::walkCode::',name),': ', ...)
                            },

                            trace = function(fn,v,w) {
                                    message(fn,': v = ',       sQuote(v),
                                               ', mode = ',    mode(v),
                                               ', typeof = ',  typeof(v),
                                               ', length = ',  length(v),
                                               ', coef   = ',  w$getCoef(v))
                            }

                            )) # end of walkCode(lhs, makeCodeWalker( ... ) )


   if ( any(idx <- is.numeric(effects) ) ) {
        stop('expression2coef',': The lhs expression ', sQuote(deparse(m.lhs)), ' ',
             'contains a numeric offset term evaluating to ', paste0(effects[idx],collapse=', '), '. ',
             'This is either an internal error or a misspecification from your part. ',
             'If so, please pull these offsets to the right-hand side of the equation')
   }

   effect.names  <- c()
   effect.coeffs <- c()

   # There might be only a single effect as in 'Agriculture = 0'. Thus use
   # c(effects) to prevent the for loop from running into an error condition
   for ( effect in c(effects) ) {
         effect.names  <- c( effect.names, as.character(effect))
         effect.coeffs <- c( effect.coeffs, attr(effect,'coef'))
   }

   list( coef        =  effect.coeffs,
         names       = effect.names,
         m           = m.rhs,
         alternative = side(ex),
         lhs         = deparse( m.lhs, width.cutoff = 500 ) )
}

### interpret character representations of linear functions
chrlinfct2matrix <- function(ex, var) {

    if (!is.character(ex))
        stop("argument ", sQuote(ex), " is not of type character")
        
    if (!is.character(var))
        stop("argument ", sQuote(var), " is not of type character")

    K <- matrix(0, nrow = length(ex), ncol = length(var))
    colnames(K) <- var
    rownames(K) <- seq_along(ex)
    m <- rep(0, length(ex))

    for (i in 1:length(ex)) {

        expr <- parse(text = ex[i])
        if (length(expr[[1]]) != 3)
            stop("argument ", sQuote(ex[i]), 
                 " cannot be interpreted as expression")

        tmp <- expression2coef(expr,vars=var)

        if (!all(tmp$names %in% var))
            stop("variable(s) ", paste(sQuote(tmp$names[!tmp$names %in% var]),collapse=', '), " not found")

        for (n in tmp$names)
            K[i, var == n] <- tmp$coef[tmp$names == n]

        m[i] <- tmp$m

        if (i == 1) {
            alternative <- tmp$alternative
        } else {
            if (tmp$alternative != alternative)
                stop("mix of alternatives currently not implemented")
        }

        rownames(K)[i] <- paste0(tmp$lhs, collapse = "")
    }
    list(K = K, m = m, alternative = alternative)
}
