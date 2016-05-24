# Author: Babak Naimi, naimi.b@gmail.com
# Date :  May 2015
# Version 1.0
# Licence GPL v3

.newFormulaFunction <- function(cls,name,args,getFeature) {
  new('.formulaFunction',cls=cls,name=name,args=args)
}

# adding the classes of formula functions into the corresponding container:
# for each function, a specific class is defined in which the name of feature, 
# and its arguments as well as the function to generate the feature from dataset is specified

.sdmFormulaFuncs <- new('.formulaFunctions')

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.auto',
                                                            representation(x='character',
                                                                           features='character',
                                                                           stat='characterORnull',
                                                                           term='call'
                                                            ))),
                                         name=c('auto','a','Auto'),args=c('x','features','stat')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.hinge',
                                                            representation(x='character',
                                                                           threshold='numeric',
                                                                           increasing='logical',
                                                                           term='call'
                                                            ))),
                                         name=c('hinge','h','H','Hinge','hing','Hing'),args=c('x','threshold','increasing')))


.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.quad',
                                                            representation(x='character',
                                                                           term='call'
                                                            ))),
                                         name=c('quad','q','Q','Quad','quadratic','Quadratic'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.cubic',
                                                            representation(x='character',
                                                                           term='call'
                                                            ))),
                                         name=c('cubic','c','C','Cubic'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.factor',
                                                            representation(x='character',
                                                                           term='call'
                                                            ))),
                                         name=c('factor','f','F','Factor','fact','Fact'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.threshold',
                                                            representation(x='character',
                                                                           threshold='numeric',
                                                                           increasing='logical',
                                                                           term='call'
                                                            ))),
                                         name=c('threshold','th','Th','thereshold','thresh','Thresh'),args=c('x','threshold','increasing')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.poly',
                                                            representation(x='character',
                                                                           degree='numeric',
                                                                           raw='logical',
                                                                           term='call'
                                                            ),
                                                            prototype(
                                                              degree=3,
                                                              raw=TRUE
                                                            ))),
                                         name=c('poly','Po','po','Poly'),args=c('x','degree','raw')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.product',
                                                            representation(x='character',
                                                                           term='call'
                                                            ))),
                                         name=c('product','p','P','Product','prod','Prod'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.func',
                                                            representation(x='call',
                                                                           term='call'
                                                            ))),
                                         name=c('I','i'),args=c('x')))



.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.simple.func',
                                                            representation(term='call'
                                                            ))),
                                         name=c('xxxxxxxx'),args=c('x')))

.sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.log',
                                                            representation(x='character',
                                                                           term='call'
                                                            ))),
                                         name=c('log','log10','exp'),args=c('x')))


.sdmFormulaFuncs$setClasses() # set the classes
# 
# 
# .newFormulaFunction <- function(cls,name,args,getFeature) {
#   new('.formulaFunction',cls=cls,name=name,args=args,getFeature=getFeature)
# }
# 
# # adding the classes of formula functions into the corresponding container:
# # for each function, a specific class is defined in which the name of feature, 
# # and its arguments as well as the function to generate the feature from dataset is specified
# 
# .sdmFormulaFuncs <- new('.formulaFunctions')
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.auto',
#                                                             representation(x='character',
#                                                                            features='character',
#                                                                            stat='characterORnull',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('auto','a','Auto'),args=c('x','features','stat'),getFeature = .getFeature.auto))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.hinge',
#                                                             representation(x='character',
#                                                                            threshold='numeric',
#                                                                            increasing='logical',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('hinge','h','H','Hinge','hing','Hing'),args=c('x','threshold','increasing'),getFeature = .getFeature.hinge))
# 
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.quad',
#                                                             representation(x='character',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('quad','q','Q','Quad','quadratic','Quadratic'),args=c('x'),getFeature = .getFeature.quad))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.cubic',
#                                                             representation(x='character',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('cubic','c','C','Cubic'),args=c('x'),getFeature = .getFeature.cubic))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.factor',
#                                                             representation(x='character',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('factor','f','F','Factor','fact','Fact'),args=c('x'),getFeature = .getFeature.factor))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.threshold',
#                                                             representation(x='character',
#                                                                            threshold='numeric',
#                                                                            increasing='logical',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('threshold','th','Th','thereshold','thresh','Thresh'),args=c('x','threshold','increasing'),getFeature = .getFeature.threshold))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.poly',
#                                                             representation(x='character',
#                                                                            degree='numeric',
#                                                                            raw='logical',
#                                                                            term='call'
#                                                             ),
#                                                             prototype(
#                                                               degree=3,
#                                                               raw=TRUE
#                                                             ))),
#                                          name=c('poly','Po','po','Poly'),args=c('x','degree','raw'),getFeature = .getFeature.poly))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.product',
#                                                             representation(x='character',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('product','p','P','Product','prod','Prod'),args=c('x'),getFeature = .getFeature.product))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.func',
#                                                             representation(x='call',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('I','i'),args=c('x'),getFeature = .getFeature.func))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.simple.func',
#                                                             representation(term='call'
#                                                             ))),
#                                          name=c('xxxxxxxx'),args=c('x'),getFeature = .getFeature.simplefunc))
# 
# .sdmFormulaFuncs$add(.newFormulaFunction(cls=quote(setClass('.log',
#                                                             representation(x='character',
#                                                                            term='call'
#                                                             ))),
#                                          name=c('log','log10','exp'),args=c('x'),getFeature = .getFeature.log))
# 
# .sdmFormulaFuncs$setClasses() # set the classes
########################################################

#------ split a formula (f) based on the sep
.split.formula <- function(f,sep='~') {
  # based on split_formula in the package Formula by Achim Zeileis
  o <- list()
  if (length(f) > 2) {
    while (length(f) > 2 && f[[1]] == sep) {
      o <- c(f[[3]],o)
      f <- f[[2]]
    }
    c(f,o)
  } else  if (length(f) == 2) {
    o[[1]] <- f[[2]]
    o
  } else {
    o[[1]] <- f
    o
  }
}
#-------

.fixFormula <- function(formula) {
  if (length(formula) == 3) {
    nsp <- trim(unlist(strsplit(as.character(formula[2]),'[+|]')))
    nsp <- nsp[nsp != '']
    if ('.' %in% nsp) {
      nsp <- nsp[nsp != '.']
      if (length(nsp) == 0) formula <- as.formula(paste(formula[1],formula[3]),env = parent.frame())
      else formula <- as.formula(paste(paste(nsp,collapse='+'),formula[1],formula[3]),env = parent.frame())
    }
  }
  formula
}
#------- excludes items in y from x
.excludeVector <- function(x,y) {
  w <- unlist(lapply(y,function(z){which(x == z)}))
  if (length(w) > 0) x <- x[-w]
  if (length(x) == 0) x <- NULL
  x
}
#-------

.getRhsFromFormula <- function(f,env=parent.frame()) {
  if (length(f) == 3) as.formula(paste('~',deparse(f[[3]])),env=env)
  else if (length(f) == 2) f
  else as.formula(paste('~',deparse(f)),env=env)
}

#---- Levenshtein distance (similarity of two strings):
.LD <- function(s,t) {
  sl <- unlist(strsplit(s,''))
  tl <- unlist(strsplit(t,''))
  if (s == t) return(0)
  else if (length(sl) == 0) return(length(tl))
  else if (length(tl) == 0) return(length(sl))
  v0 <- 0:length(tl)
  v1 <- rep(NA,length(tl)+1)
  for (i in seq_along(sl)) {
    v1[1] <- i
    for (j in seq_along(tl)) {
      if (sl[i] == tl[j]) cost <- 0
      else cost <- 1
      v1[j+1] <- min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
    }
    for (j in seq_along(v0)) {
      v0[j] <- v1[j]
    }
  }
  return(v1[length(tl)+1])
}
#-----

.LD2 <- function(s,t) {
  sl <- unlist(strsplit(s,''))
  tl <- unlist(strsplit(t,''))
  xx <- c(length(sl),length(tl))
  xx <- min(xx)/max(xx)
  if (s == t) return(0)
  else if (length(sl) == 0) return(length(tl))
  else if (length(tl) == 0) return(length(sl))
  v0 <- 0:length(tl)
  v1 <- rep(NA,length(tl)+1)
  for (i in seq_along(sl)) {
    v1[1] <- i
    for (j in seq_along(tl)) {
      if (sl[i] == tl[j]) cost <- 0
      else cost <- 1
      v1[j+1] <- min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
    }
    for (j in seq_along(v0)) {
      v0[j] <- v1[j]
    }
  }
  return(v1[length(tl)+1]*xx)
}
#------
.pmatch <- function(n,choices) {
  for (i in seq_along(n)) {
    if (n[i] != '') {
      if (!n[i] %in% choices) {
        u <- try(match.arg(tolower(n[i]),tolower(choices)),silent=TRUE)
        if (!inherits(u,"try-error")) {
          n[i] <- choices[which(tolower(choices) == u)]
        } else {
          u <- unlist(strsplit(n[i],''))
          w1 <- which(unlist(lapply(choices,function(x) tolower(strsplit(x,'')[[1]][1]) == tolower(u[1]))))
          w2 <- unlist(lapply(choices,function(x)  length(which(tolower(u) %in% tolower(strsplit(x,'')[[1]])))/length(u)))
          w4 <- unlist(lapply(choices,function(x)  .LD2(n[i],x)))
          w3 <- which(w2 > 0.5)
          if (length(w1) > 0) {
            if (length(w3) > 0) {
              w <- w1[w1 %in% w3]
              if (length(w) > 1) {
                w <- w[which(w2[w] == max(w2))]
                if (length(w) == 1) n[i] <- choices[w]
                else if (length(w1) == 1 && w2[w1] > 0.2) n[i] <- choices[w1]
                else n[i] <- NA
              } else if (length(w) == 1) {
                n[i] <- choices[w]
              } else {
                if (length(w1) == 1 && w2[w1] > 0.2) n[i] <- choices[w1]
                #else stop(paste('no match is found for',n[i]))
                else n[i] <- NA
              }
            } else {
              if (length(w1) == 1 && w2[w1] > 0.2) n[i] <- choices[w1]
              #else stop(paste('no match is found for',n[i]))
              else n[i] <- NA
            }
          } else {
            if (length(which(w2 > 0.7)) > 0) {
              w <- which(w2 > 0.7)[which(w2 > 0.7) %in% which(w4 < 3)]
              if (length(w) == 1) n[i] <- choices[w]
              #else stop(paste('no match is found for',n[i]))
              else n[i] <- NA
            } else n[i] <- NA
          } 
        }
      }
    }
  }
  
  if ('' %in% n) {
    w <- which(n == '')
    for (i in w) if (!choices[i] %in% n) n[i] <- choices[i]
    w <- which(n == '')
    if (length(w) == 1) {
      ww <- which(!choices %in% n)
      if (length(ww) == 1) n[w] <- choices[ww]
    }
  }
  #if (length(unique(n)) < length(n)) stop('repeated arguments!')
  n
}

#################---- detect the terms in the nested formula (model) inside the main formula:
.nested_terms <- function(x,r='.parent',output='.prediction') {
  n <- new('.nestedModel',response=r,output=output)
  if (length(x) > 1) {
    if (x[[1]] == '~') {
      if (length(x) == 3) {
        x <- .fixFormula(x)
        l <- .split.formula(x[[3]],'+')
        if (length(.split.formula(x[[2]],'+')) > 1) stop('nested formula in the rhs, cannot be multi-response!')
        n@response <- as.character(x[[2]])
      } else l <- .split.formula(x[[2]],'+')
    } else l <- .split.formula(x,'+')
  } else l <- list(x)
  n@terms <- lapply(l,.term)
  n
}

#--------- detect the class of the term in the formula
.term <- function(x) {
  if (length(x) == 1) {
    if (x == '.') return(new('.all.vars',names=as.character(x)))
    else return(new('.var',name=as.character(x)))
  } else if (length(x) > 1) {
    if (x[[1]] == 'm') {
      .nested_terms(x[[2]],output='.prediction')
    } else if (x[[1]] == 'r') {
      .nested_terms(x[[2]],output='.residual')
    } else if (x[[1]] == 'select' || x[[1]] == 'se') {
      .select.terms(x)
    } else if (x[[1]] == 'coords' || any(!is.na(pmatch(c("co"),tolower(as.character(x[[1]])))))) {
      .exCoords(x)
    } else if (as.character(x[[1]]) %in% c('g','G','gr','group','Group','GROUP','GR','gro','grop','grp')) {
      .exGroup(x)
    } else if (as.character(x[[1]]) %in% c('time','Time','tim','TIME','Tim')) {
      .exTime(x)
    } else if (as.character(x[[1]]) %in% c('info','Info','inf','INFO')) {
      .exInfo(x)
    } else if (as.character(x[[1]]) %in% c('*',':','product','p','P','Product','prod','Prod','pro','Pro')) {
      .exProduct(x)
    } else .exFunc(x)
  }
}
#-------
.exGroup <- function(x) {
  new('.grouping',group.var=as.character(x[[2]]),term=x)
}
#------
.exTime <- function(x) {
  xx <- as.character(x[[1]])
  s <- list()
  n <- names(x)
  if (!is.null(n)) {
    for (i in 2:length(n)) {
      if (n[i] != '') s[[n[i]]] <- x[[i]]
      else s[[(i-1)]] <- x[[i]]
    }
  } else {
    for (i in 2:length(x)) s[[(i-1)]] <- x[[i]]
  }
  new('.time',name=xx,terms=s,term=x)
}
#--------
.exInfo <- function(x) {
  n <- NULL
  if (length(x[[2]]) > 1 && as.character(x[[2]][[1]]) %in% c('|','+')) {
    if (as.character(x[[2]][[1]]) == '+') n <- as.character(.split.formula(x[[2]],'+'))
    else n <- as.character(.split.formula(x[[2]],'|'))
  } else {
    n <- as.character(x[[2]])
  } 
  if (!is.null(n)) new('.Info',names=n)
}


#---------
.exFunc <- function(x) {
  xx <- as.character(x[[1]])
  mn <- .sdmFormulaFuncs$funcNames
  names(mn) <- mn
  if (xx %in% mn) xx <- names(mn)[mn == xx]
  else {
    mnlist <- lapply(mn,function(x) .sdmFormulaFuncs$funcs[[x]]@name)
    u <- unlist(lapply(mnlist,function(x) xx %in% x))
    if (any(u)) xx <- names(u)[which(u)]
    else xx <- NULL
  }
  if (!is.null(xx)) {
    ss <- .sdmFormulaFuncs$getFuncs(xx)
    a <- ss[[xx]]@args
    cls <- new(ss[[xx]]@cls[[2]])
    n <- names(x)
    n <- n[2:length(n)]
    
    if (!is.null(n)) {
      n <- .pmatch(n,a)
      if (length(n[n != ''] > 0) && !all(n[n != ''] %in% a)) stop(paste0('some arguments in function ',xx, ' is unknown!'))
      for (i in 1:length(n)) {
        if (n[i] != '') {
          if (class(x[[i+1]]) == 'name') slot(cls,n[i]) <- as.character(x[[i+1]])
          else slot(cls,n[i]) <- x[[i+1]]
        } else {
          if (class(x[[i+1]]) == 'name') slot(cls,a[i]) <- as.character(x[[i+1]])
          else slot(cls,a[i]) <- x[[i+1]]
        } 
      }
    } else {
      for (i in 2:length(x)) {
        if (class(x[[i]]) == 'name') slot(cls,a[i-1]) <- as.character(x[[i]])
        else slot(cls,a[i-1]) <- x[[i]]
      }
    }
  } else {
    if (exists(as.character(x[[1]]),mode='function')) {
      cls <- new('.simple.func')
    } else stop(paste(as.character(x[[1]]),'is not a known function!'))
  }
  
  cls@term <- x
  cls
}
#--------
.exProduct <- function(x) {
  cls <- new('.product')
  if (as.character(x[[1]]) %in% c('*',':')) {
    cls@x <- as.character(.split.formula(x,as.character(x[[1]])))
  } else {
    if (length(x) == 2) {
      xx <- x[[2]]
      cls@x <- as.character(.split.formula(xx,as.character(xx[[1]])))
    } else {
      xx <- c()
      for (i in 2:length(x)) xx <- c(xx,as.character(x[[i]]))
      cls@x <- xx
    }
  }
  
  cls@term <- x
  cls
}
#--------
.exCoords <- function(x) {
  if (length(x[[2]]) > 1 && as.character(x[[2]][[1]]) %in% c('|','+')) {
    xy <- as.character(x[[2]])[2:3]
  } else if (length(x) == 3) {
    xy <- c(as.character(x[[2]]),as.character(x[[3]]))
  } else stop('in formula, coordinates are not properly defined; Example: ~...+coords(x+y)+...')
  new('.coord.vars',xy=xy)
}


#-------------
.select.terms <- function(x) {
  a <- c('formula','n','stat','th','keep')
  s <- list()
  n <- names(x)
  if (length(n) > 0) {
    for (i in seq_along(n)) {
      if (n[i] != '') {
        if (any(!is.na(pmatch(c("n"),tolower(n[i]))))) n[i] <- 'n'
        else if (any(!is.na(pmatch(c("st"),tolower(n[i]))))) n[i] <- 'stat'
        else if (any(!is.na(pmatch(c("th"),tolower(n[i]))))) n[i] <- 'th'
        else if (any(!is.na(pmatch(c("k"),tolower(n[i]))))) n[i] <- 'keep'
      }
    }
  } else n <- rep('',length(x))
  
  if (length(n[n != ''] > 0) && !all(n[n != ''] %in% a)) stop('some arguments in select function is unknown!')
  if (length(x) > 5) stop('the arguments in select function are not match!')
  for (i in 2:length(n)) {
    if (n[i] != '') s[[n[i]]] <- x[[i]]
    else s[[a[i-1]]] <- x[[i]]
  }
  
  if (length(s[['formula']]) > 1) {
    if (s[['formula']][[1]] != '|' && s[['formula']][[1]] != 'select') stop('something in `select` is wrong, check the help to see how the select function in formula should be used...')
    else {
      l <- .split.formula(s[['formula']],'|')
      if (any(unlist(lapply(l,function(x) {length(x) > 1 && x[[1]] == '+'})))) stop('something wrong with select; example: select(var1|var2|m(var1+var2),n=2,stat="auc"')
    }
  } else l <- list(s[['formula']])
  
  n <- new('.selectFrame')
  
  if (!is.null(s[['n']])) n@n <- s[['n']]
  if (!is.null(s[['stat']])) n@stat <- s[['stat']]
  if (!is.null(s[['keep']])) {
    k <- as.character(s[['keep']])
    if (k[1] == '~') {
      k <- strsplit(k[2],'\\+')[[1]]
      for (i in seq_along(k)) k[i] <- .trim(k[i])
    } else if (k[1] == '+') {
      k <- .split.formula(s[['keep']],'+')
    } else if (k[1] == 'c') k <- k[2:length(k)]
    n@keep <- k
  }
  n@sets <- lapply(l,.term)
  n
}
#--------
.trim <- function(x) {
  x <- strsplit(x,'')[[1]]
  paste(x[x != ' '],collapse='')
}
#-------

########

# .exFormula extract terms in formula and detect what each term is. it may be a model.term (including a
# variable, a function, a nested model, etc.) or a data.term (including coordinates, select function, group, etc.)
.exFormula <- function(f,data,detect=TRUE) {
  f <- .fixFormula(f)
  v <- colnames(data)
  
  nall <- n <- all.vars(f)
  
  if ('.' %in% n) {
    n <- n[-which(n == '.')]
    nall <- unique(c(v,n))
  }
  
  nFact <- v[.where(is.factor,data)]
  if (length(nFact) == 0) nFact <- NULL
  else {
    if (any(nFact %in% nall)) nFact <- nFact[nFact %in% nall]
    else nFact <- NULL
  }
  
  sf <- .split.formula(f)
  lhs <- rhs <- NULL
  if (length(sf) > 2) stop('in the right hand side of the formula, the `~` can only be in m(...)')
  f <- new('sdmFormula',formula=f)
  if (length(sf) == 1) rhs <- sf[[1]]
  else {
    lhs <- sf[[1]]
    rhs <- sf[[2]]
  }
  
  if (!is.null(lhs)) {
    lhs <- .split.formula(lhs,'+')
    nsp <- as.character(lhs)
    nall <- .excludeVector(nall,nsp)
    n <- .excludeVector(n,nsp)
    nFact <- .excludeVector(nFact,nsp)
  } else {
    if (detect) {
      w <- which(unlist(lapply(data,.isBinomial)))
      if (length(w) > 0) {
        nsp <- v[w]
        nsp <- nsp[!nsp %in% n]
        nall <- .excludeVector(nall,nsp)
        nFact <- .excludeVector(nFact,nsp)
        lhs <- as.list(nsp)
      } else nsp <- NULL
    } else nsp <- NULL
    
  }
  f@vars <- nall
  f@species <- as.character(lhs)
  
  if (length(rhs) == 2) rhsi <- list(rhs)
  else rhsi <- .split.formula(rhs,'+')
  
  nxy <- NULL
  temp <- unlist(lapply(rhsi,function(x) as.character(x)[[1]] == 'coords'))
  if (any(temp)) nxy <- as.character(.split.formula(rhsi[[which(temp)]][[2]],'+'))
  
  vars <- .excludeVector(nall,c(n,nxy))
  
  w <- unlist(lapply(rhsi,function(x) x == '.'))
  if (any(w)) {
    if (!is.null(vars)) rhsi <- c(rhsi[!w],.split.formula(as.formula(paste('~',paste(vars,collapse='+')))[[2]],'+')) 
    else rhsi <- rhsi[!w]
  }
  
  if (!is.null(nFact)) {
    for (i in seq_along(rhsi)) {
      if (length(as.character(rhsi[[i]])) == 1) {
        if (any(nFact == as.character(rhsi[[i]]))) {
          rhsi[[i]] <- as.formula(paste('~f(',nFact[nFact == as.character(rhsi[[i]])],')'))[[2]]
        }
      }
    }
  }
  
  
  func.cls <- unlist(lapply(.sdmFormulaFuncs$funcNames,function(x) .sdmFormulaFuncs$funcs[[x]]@cls[[2]]))
  temp <- lapply(rhsi,.term)
  w <- unlist(lapply(temp,class))
  wt <- which(w %in% c('.var','.nestedModel',func.cls))
  if (length(wt) > 0) f@model.terms <- temp[wt]
  wt <- which(w %in% c('.coord.vars','.grouping','.Info','.time'))
  if (length(wt) > 0) f@data.terms <- c(f@data.terms,temp[wt])
  wt <- which(w %in% c('.selectFrame'))
  if (length(wt) > 0) {
    for (i in wt) {
      if (temp[[i]]@stat %in% c('vif','pca')) f@data.terms <- c(f@data.terms,temp[wt])
      else f@model.terms <- c(f@model.terms,temp[wt])
    }
  }
  f
}
#-----------
###################

