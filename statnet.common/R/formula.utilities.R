#  File R/formula.utils.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
###################################################################
## This file has utilities whose primary purpose is examining or ##
## manipulating ERGM formulas.                                   ##
###################################################################

## This function appends a list of terms to the RHS of a formula. If
## the formula is one-sided, the RHS becomes the LHS, if
## keep.onesided==FALSE (the default).
## For example,
## append.rhs.formula(y~x,list(as.name("z1"),as.name("z2"))) -> y~x+z1+z2
## append.rhs.formula(~y,list(as.name("z"))) -> y~z
## append.rhs.formula(~y+x,list(as.name("z"))) -> y+x~z
## append.rhs.formula(~y,list(as.name("z")),TRUE) -> ~y+z
append.rhs.formula<-function(object,newterms,keep.onesided=FALSE){
  for(newterm in newterms){
    if(length(object)==3) object[[3]]<-call("+",object[[3]],newterm)
    else if(keep.onesided) object[[2]]<-call("+",object[[2]],newterm)
    else object[[3]]<-newterm
  }
  object
}

# A reimplementation of update.formula() that does not simplify.  Note
# that the resulting formula's environment is set as follows: If
# from.new==FALSE, it is set to that of object. Otherwise, a new
# sub-environment of object, containing, in addition, variables in new
# listed in from.new (if a character vector) or all of new (if TRUE).
nonsimp.update.formula<-function (object, new, ..., from.new=FALSE){
  old.lhs <- if(length(object)==2) NULL else object[[2]]
  old.rhs <- if(length(object)==2) object[[2]] else object[[3]]
  
  new.lhs <- if(length(new)==2) NULL else new[[2]]
  new.rhs <- if(length(new)==2) new[[2]] else new[[3]]
  
  sub.dot <- function(c, dot){
    if(is.null(dot)) c # If nothing to substitute with, just return it.
    else if(is.call(c)) as.call(c(list(c[[1]]), lapply(c[-1], sub.dot, dot))) # If it's a call, construct a call consisting of the call and each of the arguments with the substitution performed, recursively.
    else if(is.name(c) && c==".")  dot # If it's a dot, return substitute.
    else c # If it's anything else, just return it.
  }
  
  deparen<- function(c, ops = c("+","*")){
    if(is.call(c)){
      if(as.character(c[[1]]) %in% ops){
        op <- as.character(c[[1]])
        if(length(c)==2 && is.call(c[[2]]) && c[[2]][[1]]==op)
          return(deparen(c[[2]], ops))
        else if(length(c)==3 && is.call(c[[3]]) && c[[3]][[1]]==op)
          return(call(op, call(op, deparen(c[[2]],ops), deparen(c[[3]][[2]],ops)), deparen(c[[3]][[3]],ops)))
      }
      return(as.call(c(list(c[[1]]), lapply(c[-1], deparen, ops)))) # If it's a non-reducible call, construct a call consisting of the call and each of the arguments with the substitution performed, recursively.
    }else return(c)
  }
  
  out <- if(length(new)==2) call("~", deparen(sub.dot(new.rhs, old.rhs))) else call("~", deparen(sub.dot(new.lhs, old.lhs)), deparen(sub.dot(new.rhs, old.rhs)))

  #  a new sub-environment for the formula, containing both
  # the variables from the old formula and the new.
  
  if(identical(from.new,FALSE)){ # The new formula will use the environment of the original formula (the default).
    e <- environment(object)
  }else{
    # Create a sub-environment also containing variables from environment of new.
    e <- new.env(parent=environment(object))
    
    if(identical(from.new,TRUE)) from.new <- ls(pos=environment(new)) # If TRUE, copy all of them (dangerous!).
    
    for(name in from.new)
      assign(name, get(name, pos=environment(new)), pos=e)
  }

  as.formula(out, env = e)
}

term.list.formula<-function(rhs, sign=+1){
  if(length(rhs)==1) {attr(rhs,"sign")<-sign; list(rhs)}
  else if(length(rhs)==2 && rhs[[1]]=="+") term.list.formula(rhs[[2]],sign)
  else if(length(rhs)==2 && rhs[[1]]=="-") term.list.formula(rhs[[2]],-sign)
  else if(length(rhs)==3 && rhs[[1]]=="+") c(term.list.formula(rhs[[2]],sign),term.list.formula(rhs[[3]],sign))
  else if(length(rhs)==3 && rhs[[1]]=="-") c(term.list.formula(rhs[[2]],sign),term.list.formula(rhs[[3]],-sign))
  else if(rhs[[1]]=="(") term.list.formula(rhs[[2]], sign)
  else {attr(rhs,"sign")<-sign; list(rhs)}
}

