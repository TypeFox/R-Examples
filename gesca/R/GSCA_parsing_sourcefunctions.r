###########
# GeSCA R packging project
# Formula parsing script
###########
# Author: Sungyoung Lee
# Date  : Nov. 11 2015
# E-mail: me@lsy.io
###########

#############
#
# EXTERNAL FUNCTION PART
#
#############

###########
# gesca.run
###########
# Objective : Converts string formula into GeSCA matrices
#            with given data matrix
# Accepts:
#  - s.fml      : String, multiple lines of formulas with specified format
#  - data     : n * p matrix, where n=#sample, p=#manifest
#  - n.param    : Real value, must specify the value to denote 'parameter'
#  - group.equal: String vector, must among 'loadings' and 'paths'
#  - nbt        : # of bootstraps to perform
#  - ceps       : Convergence threshold
#  - itmax      : Maximum number of iterations to try
#
# Requirements:
#  - Formulas in s.fml DO NOT need to be 'ordered'
#  - Names of latent variables must be compatible with C-style naming
#  - Lines with whitespaces only or with comments will be skipped
#  - Whitespaces within/between formula are allowed
#  - Comments must starts with # character
#  - Not all of variables in data are required to consist formula
#
# Returns: Return values of gsca.mg / gsca.mg.ho2
#
###########
# Note: This function halts execution when there is an error during
# formula parsing - so an error can be caught using try() for example.
###########

gesca.run <- function(myModel, data, group.name=NULL, group.equal=NULL, nbt=100, itmax=100, ceps=1e-5, moption=0, missingvalue=-9999) {
	if (!is.character(myModel) || length(myModel) > 1 || nchar(myModel) == 0)
		stop("formula should be a string")
		
	if (!is.data.frame(data) && !is.matrix(data) || !nrow(data) || !ncol(data))
		stop("data should be a non-empty matrix or dataframe")

	if (moption != 0) {
		if (moption != 1 && moption != 2 && moption != 3)
			stop("moption should be 1, 2 or 3")
		cat("Missing value is defined to [", missingvalue, "]\n")
	} else
		moption = 0
		
	if (!is.null(group.name)) {
		orig.gname <- group.name
		if (length(group.name) > 1)
			stop("group.name should be a column name or an index")
			
		if (is.na(match(group.name, colnames(data)))) {
			if ((group.name - as.integer(group.name)) != 0 || ncol(data) < group.name || group.name < 0)
				stop("group.name should be a column name or an index")
		} else
			group.name <- match(group.name, colnames(data))

		# Get group var    
		group.var <- as.character(data[,group.name])
		# Sort
		data = data[order(group.var),]
	} else {
		group.var <- rep(1, nrow(data))
		orig.gname <- NULL
	}
	
	ret <- parse.formula(myModel, m.data=data, group.equal=group.equal)
  
	# Have high-order
	if (!is.null(ret$W002))
		msg <- capture.output(res <- gsca.mg.ho2(ret$Z0, group.var, ret$W00, ret$W002, ret$C00, ret$B00,
		ret$loadtype, ret$loadtype2, ceps=ceps, nbt=nbt, itmax=itmax, moption=moption, missingvalue=missingvalue))
	else
		msg <- capture.output(res <- gsca.mg(ret$Z0, group.var, ret$W00, ret$C00, ret$B00, ret$loadtype, ceps=ceps, nbt=nbt, itmax=itmax, moption=moption, missingvalue=missingvalue))

	res$niter        <- as.numeric(strsplit(msg, " ")[[1]][6])
	res$eps          <- ceps
	res$wname        <- colnames(ret$Z0)
	
	if (!is.null(ret$W002))
		res$lname    <- c(names(ret$loadtype), names(ret$loadtype2))
	else
		res$lname    <- names(ret$loadtype)
		
	res$grp          <- names(table(group.var))
	res$gname        <- orig.gname
	res$B00          <- ret$B00
	res$C00          <- ret$C00
	res$W00          <- ret$W00
	res$nbt          <- nbt
	
	if (!is.null(ret$W002)) {
		res$loadtype <- c(ret$loadtype, ret$loadtype2)
		res$W002     <- ret$W002
	} else
		res$loadtype <- ret$loadtype
		
	class(res)       <- c("gesca", class(res))
	invisible(res)
}

###########
# parse.func
###########
# Objective : Converts string element in a formula into
#              parsed return value by considering it as a function call
# Accepts:
#  - fun.cache: 
#  - str      : 
#  - idx      : 
#
# Requirements:
#  - Formulas in myModel DO NOT need to be 'ordered'
#  - Names of latent variables must be compatible with C-style naming
#  - Lines with whitespaces only or with comments will be skipped
#  - Whitespaces within/between formula are allowed
#  - Comments must starts with # character
#
# Acronyms
#  - #l : Number of latent variables
#  - #L : Number of high-order latent variables
#  - #m : Number of manifest variables
#
# Returns: a list consist of multiple variables
#  - idxmanifest: Length (#l) list, indices of manifest variables
#                 for each latent variable
#  - idxlatent  : Length (#L) list, indices of latent variables for
#                 each high-order latent varaible
#  - loadtype   : Length (#l) vector, type of latent
#                 (1=reflective, 0=formative)
#  - loadtype2  : Length (#L) vector, type of high-order latent
#                 (1=reflective, 0=formative)
#  - W00        : (#m) * (#l) matrix for W001
#  - W002       : (#l) * (#L) matrix for W002
#  - B00        : (#l+#L) * (#l+#L) matrix for B00
#  - C00        : t(W00)
###########
# Note: This function halts execution when there is an error during
# formula parsing - so an error can be caught using try() for example.
###########
parse.func <- function(fun.cache, str, idx) {
	ret <- NULL

	# Extract function part
	# tfun[1] == function name
	# tfun[2] == function argument
	tfun <- strsplit(gsub("^(\\w+)\\(([^\\)]+)\\)$", "\\1||\\2", str), "\\|\\|")[[1]]


	if (tfun[1] == "c") {
		# tfun[2] must x,x
		if (nchar(tfun[2]) != 3 || substr(tfun[2], 2, 2) != "," ||
		substr(tfun[2], 1, 1) != substr(tfun[2], 3, 3) ||
		!is.na(match(substr(tfun[2], 1, 1), fun.cache[[tfun[1]]])))
			stop("Invalid function [", str, "]")
			
		fun.cache[[tfun[1]]] <- c(fun.cache[[tfun[1]]], substr(tfun[2], 1, 1))
		ret <- idx
		idx <- idx + 1
	} else
		stop("Invalid function [", str, "]")

	return(list(
		val = ret,
		idx = idx,
		fun.cache = fun.cache
	))
}

###########
# parse.formula
###########
# Objective : Converts string formula into GeSCA matrices
#            with given data matrix
# Accepts:
#  - s.fml      : String, multiple lines of formulas with specified format
#  - m.data     : n * p matrix, where n=#sample, p=#manifest
#  - b.debug    : TRUE/FALSE, generating debug output if TRUE
#  - n.param    : Real value, must specify the value to denote 'parameter'
#  - group.equal: String vector, must among 'loadings' and 'paths'
#
# Requirements:
#  - Formulas in s.fml DO NOT need to be 'ordered'
#  - Names of latent variables must be compatible with C-style naming
#  - Lines with whitespaces only or with comments will be skipped
#  - Whitespaces within/between formula are allowed
#  - Comments must starts with # character
#
# Acronyms
#  - #l : Number of latent variables
#  - #L : Number of high-order latent variables
#  - #m : Number of manifest variables
#
# Returns: a list consist of multiple variables
#  - idxmanifest: Length (#l) list, indices of manifest variables
#                 for each latent variable
#  - idxlatent  : Length (#L) list, indices of latent variables for
#                 each high-order latent varaible
#  - loadtype   : Length (#l) vector, type of latent
#                 (1=reflective, 0=formative)
#  - loadtype2  : Length (#L) vector, type of high-order latent
#                 (1=reflective, 0=formative)
#  - W00        : (#m) * (#l) matrix for W001
#  - W002       : (#l) * (#L) matrix for W002
#  - B00        : (#l+#L) * (#l+#L) matrix for B00
#  - C00        : t(W00)
###########
# Note: This function halts execution when there is an error during
# formula parsing - so an error can be caught using try() for example.
###########

parse.formula <- function(
  s.fml,
  m.data      = NULL,
  b.debug     = FALSE,
  n.param     = 99,
  group.equal = NULL
  ) {
  opt.debug <- options("b.debug")[[1]]
  if (length(opt.debug) && opt.debug == TRUE) b.debug <- TRUE
  if (is.null(m.data) && !is.null(options("m.data"))) m.data <- options("m.data")[[1]]
  a.col <- colnames(m.data)
  if (length(a.col) == 0)
    a.col <- paste("V", 1:ncol(m.data), sep="")
    
  # group.equal check
  b.equalL <- FALSE
  b.equalP <- FALSE
  if (!is.null(group.equal)) {
    i.valid <- match(group.equal, c("loadings", "paths"))
	if (any(is.na(i.valid)))
	  stop(paste("Invalid group.equal value [", paste(group.equal, collapse=" "), "]"))
	if (any(i.valid == 1)) b.equalL <- TRUE
	if (any(i.valid == 2)) b.equalP <- TRUE
  }
  
  # Divide and remove whitespaces
  a.fml <- gsub("\\s", "", strsplit(s.fml, "\n")[[1]])
  a.latent <- c()
  
  # Remove comments
  i.cmt <- grep("^#", a.fml)
  if (length(i.cmt)) {
    a.fml[i.cmt] <- ""
    if (b.debug) cat("Comments removed:", i.cmt, "\n")
  }

  # Survived data
  if (b.debug) print(a.fml)

  # Parse constant (==)
  i.const <- grep("==", a.fml)
  a.const <- c()
  if (length(i.const)) {
    # Divide
	a.const <- getFml(a.fml[i.const], "==")

	####
	# Left-side variables check
	####
	a.left <- unlist(lapply(a.const, function(a)a[1]))
    a.nice <- grep("^[a-zA-Z_]\\w*$", a.left)
    
    # stop: invalid naming
    if (length(a.left) != length(a.nice)) {
	  a.err <- a.left[-a.nice]
	  stop(paste("The following constants are invalid:", paste(a.err, collapse=" ")))
	}

	####
	# Right-side variables check
	####
	a.right  <- as.numeric(unlist(lapply(a.const, function(a)a[2])))
    i.notnum <- which(is.na(a.right))
	if (length(i.notnum)) {
	  a.err <- a.right[i.notnum]
	  stop(paste("The following constants are not numeric:", paste(a.err, collapse=" ")))
	}

	####
	# Reallocate
	####
	a.const <- a.right
	names(a.const) <- a.left

	if (b.debug == TRUE) {
	  cat("The following constant(s) defined:\n")
	  a.nameconst <- names(a.const)
	  for (i in 1:length(a.const)) {
	    cat("  ", a.nameconst[i], "=", a.const[i], "\n")
	  }
	}
    
    # Remove parsed data
    a.fml[i.const] <- ""
  }

  # Cache for function-wise static variables
  fun.cache <- list()

  # Parse =~ first
  i.outer <- grep("=~", a.fml)
  a.outer <- c()
  a.ridx  <- list()
  idx.C00 <- 1
  if (length(i.outer)) {
    # Divide
    a.outer <- getFml(a.fml[i.outer], "=~")

    ####
    # Left-side variables check
    ####
    a.left <- unlist(lapply(a.outer, function(a)a[1]))
    a.nice <- unique(c(
      grep("^[a-zA-Z_]\\w*$", a.left),
      grep("^[a-zA-Z_]\\w*\\([01]\\)$", a.left)
    ))

    # stop: Some variable is not right
    if (length(a.nice) != length(a.left)) {
      a.err <- a.left[-a.nice]
      stop(paste("The following latents are invalid:", paste(a.err, collapse=" ")))
    }
  
    # Make left guys to latent
    a.latent <- rep(1, length(a.left))
    names(a.latent) <- gsub("^([a-zA-Z_]\\w*).*$", "\\1", a.left)
    
    # Determine 0/1
    i.0 <- which(gsub("^[a-zA-Z_]\\w*", "", a.left) == "(0)")
    a.latent[i.0] <- 0
    
    # Check left duplication
    a.tbleft <- table(names(a.latent))
    
    # stop: duplicated naming of latent
    if (any(a.tbleft > 1)) {
      stop(paste("Duplicated latent(s) definition found:",
	    names(a.tbleft)[a.tbleft>1]))
    }
    
    # Check left-mv duplication
    i.dup <- which(!is.na(match(names(a.latent), a.col)))
    
    # stop: latent variable name(s) dup. w/ manifest variables
    if (length(i.dup))
      stop(paste("Duplicated latent naming(s) with manifest variables found:", paste(names(a.latent)[i.dup], collapse=" ")))
  
    ####
    # Right-side variables check
    ####
    a.right <- unlist(lapply(a.outer, function(a)a[2]))
    a.rcol <- lapply(strsplit(a.right, "\\+"), function(a)unlist(lapply(strsplit(a, "\\*"), function(b)b[length(b)])))
    a.rall <- unique(unlist(a.rcol))
    i.namechk <- which(is.na(match(a.rall, a.col)))
    
    # stop: invalid column
    if (length(i.namechk))
      stop(paste("The following column(s) in the formula do not exist:", paste(a.rall[i.namechk], collapse=" "), "\nColumns:", paste(a.col, collapse=" ")))

    ### RETURN : W00
    r.C00 <- matrix(0, nrow=length(a.latent), ncol=length(a.col))
    I <- 1
    for (i in a.outer) {
      i.x <- I
      I <- I + 1
      
      # For each element in the right part of the equation
	  for (j in strsplit(i[2], "\\+")[[1]]) {
	    #cat("J value\n")
		#print(j)
        
        # If the element contains *, it should be a constraint
	    if (length(grep("\\*", j))) {
          # Assumes the former part is contraint/function, while the latter part is the variable
		  tmp <- unlist(strsplit(j, "\\*"))

		  # Try to convert it num
		  tvar <- suppressWarnings(as.numeric(tmp[1]))
          
		  # Try to convert it const if fails
		  if (is.na(tvar))
		    tvar <- a.const[tmp[1]]
            
		  # Try to convert it function
		  if (length(grep("^\\w+\\([^\\)]+\\)$", tmp[1]))) {
            tmp.ret   <- parse.func(fun.cache, tmp[1], idx.C00)
            #print(tmp.ret)
            tvar      <- tmp.ret$val
            idx.C00   <- tmp.ret$idx
            fun.cache <- tmp.ret$fun.cache
		  }
		  if (is.null(tvar) || is.na(tvar))
            stop(paste("The modifier [", tmp[1], "] does not recognized!"))

		  i.y   <- match(tmp[2], a.col)
		  n.par <- tvar
		} else {
		  tvar  <- ifelse(b.equalL, idx.C00, n.param)
		  if (tvar != 99 && tvar != 0 && a.latent[i.x] != 0)
			idx.C00 <- idx.C00 + 1
		  i.y   <- match(j, a.col)
		  n.par <- tvar
		}
        if (a.latent[i.x] == 0) # Set 0 if formative
          n.par <- 0
#		else if (n.par == 99)
#		  idx.C00 <- idx.C00 - 1 # Actually make it increase
        r.C00[i.x, i.y] <- n.par
      }
    }
    #browser()
    #print(a.col)
    
    a.ridx <- lapply(a.rcol, function(a)match(a, a.col))
    names(a.ridx) <- names(a.latent)
    for (i in 1:length(a.ridx))
      attr(a.ridx[[i]], "reflective") <- a.latent[i]

    ### RETURN : loadtype0
    r.loadtype <- a.latent
    ### RETURN : W00
    r.W00 <- matrix(0, nrow=length(a.latent), ncol=length(a.col))
    for (i in 1:length(a.ridx))
      r.W00[i,a.ridx[[i]]] <- n.param
    r.W00 <- t(r.W00)
    ### RETURN : C00
	#if (b.equalL == TRUE) {
	#  vidx <- which(r.W00 == n.param)
	#  r.C00 <- r.W00
	#  r.C00[vidx] <- 1:length(vidx)
	#  r.C00 <- t(r.C00)
	#} else
	#  r.C00 <- t(r.W00)
    #print(r.W00)
    
    # Remove parsed data
    a.fml[i.outer] <- ""
  } else {
    r.loadtype <- NULL
    r.W00      <- NULL
    r.C00      <- NULL
  }
  n.W00 <- sum(r.W00 == n.param)
  
  # Parse =: next
  i.higher <- grep("=:", a.fml)
  a.higher <- c()
  a.hidx  <- list()
  if (length(i.higher)) {
    # Divide
    a.higher <- getFml(a.fml[i.higher], "=:")

    ####
    # Left-side variables check
    ####
    a.left <- unlist(lapply(a.higher, function(a)a[1]))
    a.nice <- unique(c(
      grep("^[a-zA-Z_]\\w*$", a.left),
      grep("^[a-zA-Z_]\\w*\\([01]\\)$", a.left)
    ))

    # stop: Some variable is not right
    if (length(a.nice) != length(a.left)) {
      a.err <- a.left[-a.nice]
      stop(paste("The following high-order latents are invalid:", paste(a.err, collapse=" ")))
    }
  
    # Make left guys to latent
    a.holatent <- rep(1, length(a.left))
    names(a.holatent) <- gsub("^([a-zA-Z_]\\w*).*$", "\\1", a.left)
    
    # Determine 0/1
    i.0 <- which(gsub("^[a-zA-Z_]\\w*", "", a.left) == "(0)")
    a.holatent[i.0] <- 0
  
    # Check left duplication
    a.tbleft <- table(names(a.holatent))
    
    # stop: duplicated naming of latent
    if (any(a.tbleft > 1)) {
      stop(paste("Duplicated high-order latent(s) definition found:",
	    names(a.tbleft)[a.tbleft>1]))
    }
    #print(a.holatent)
    
    # Check left-mv duplication / left-existing lv duplication
    i.dup <- union(which(!is.na(match(names(a.holatent), a.col))),
      which(!is.na(match(names(a.holatent), names(a.ridx)))))
    
    # stop: latent variable name(s) dup. w/ manifest variables
    if (length(i.dup))
      stop(paste("Duplicated high-order latent naming(s) with manifest/latent variables found:", paste(names(a.holatent)[i.dup], collapse=" ")))
  
    ####
    # Right-side variables check
    ####
    a.right <- unlist(lapply(a.higher, function(a)a[2]))
    a.rval <- strsplit(a.right, "\\+")
    a.rall <- unique(unlist(a.rval))
        
    i.namechk <- which(is.na(match(a.rall, names(a.ridx))))
    
    # stop: invalid column
    if (length(i.namechk))
      stop(paste("The following latents in the formula do not exist:", paste(a.rall[i.namechk], collapse=" ")))
      
    a.hridx <- lapply(strsplit(a.right, "\\+"), function(a)match(a, names(a.ridx)))
    names(a.hridx) <- names(a.holatent)
    for (i in 1:length(a.hridx))
      attr(a.hridx[[i]], "reflective") <- a.holatent[i]

    ### RETURN : loadtype0
    r.loadtype2 <- a.holatent
    ### RETURN : W00
    r.W002 <- matrix(0, nrow=length(a.hridx), ncol=length(a.latent))
    for (i in 1:length(a.hridx))
      r.W002[i,a.hridx[[i]]] <- n.param
    r.W002 <- t(r.W002)
    
    # Remove parsed data
    a.fml[i.higher] <- ""
  } else
    a.hridx <- c()

  # Parse ~ second
  idx.B00 <- idx.C00
  idx.const.B00 <- matrix(0, 0, 2)
  #browser()
  a.ltname <- c(names(a.ridx), names(a.hridx))
  i.inner <- setdiff(grep("~", a.fml), i.outer)
  r.B00 <- matrix(0, nrow=length(a.ltname), ncol=length(a.ltname))
  if (length(i.inner)) {
    a.inner <- getFml(a.fml[i.inner], "~")

    ####
    # Right-side variables check
    ####
    a.right <- unlist(lapply(a.inner, function(a)a[2]))
	a.rval <- strsplit(a.right, "\\+")
	a.rvalwoV <- lapply(a.rval, function(v)gsub("^[^\\*]+\\*", "", v))
    a.rall <- unique(unlist(a.rvalwoV))

    ####
    # Variables check
    ####
    a.lr <- unique(c(unlist(lapply(a.inner, function(a)a[1])), a.rall))
    a.nice <- grep("^[a-zA-Z_]\\w*$", a.lr)
    
    # stop: invalid naming
    if (length(a.lr) != length(a.nice))
      stop(paste("The following latent(s) have invalid naming:", paste(a.lr[-a.nice], collapse=" ")))
    #print(a.lr)
    
    # Duplication check w/ manifest variables
    a.lr.uniq <- unique(a.lr)
    i.lr.dup <- which(!is.na(match(a.lr.uniq, a.col)))
    
    # stop: latent variable name(s) dup. w/ manifest variables
    if (length(i.lr.dup))
      stop(paste("Duplicated latent naming(s) with manifest variables found:", paste(a.lr.uniq[i.lr.dup], collapse=" ")))

    # All must exist on l-m
    b.sane <- !is.na(match(a.lr, a.ltname))
    i.invalid <- which(b.sane == FALSE)
    
    # stop: Invalid left-side variables found
    if (length(i.invalid)) {
      stop(paste("Invalid latent(s) found:", paste(unique(a.lr[i.invalid]), collapse=" ")))
    }
    
    # Remove parsed data
    a.fml[i.inner] <- ""
    
    ### RETURN : B00
	#print(a.inner)
    for (i in a.inner) {
      i.x <- match(i[1], a.ltname)
      
      # For each element in the right part of the equation
	  for (j in strsplit(i[2], "\\+")[[1]]) {
	    #cat("J value\n")
		#print(j)
        
        # If the element contains *, it should be a constraint
	    if (length(grep("\\*", j))) {
          # Assumes the former part is contraint/function, while the latter part is the variable
		  tmp <- unlist(strsplit(j, "\\*"))

		  # Try to convert it num
		  tvar <- suppressWarnings(as.numeric(tmp[1]))
          
		  # Try to convert it const if fails
		  if (is.na(tvar))
		    tvar <- a.const[tmp[1]]
            
		  # Try to convert it function
		  if (length(grep("^\\w+\\([^\\)]+\\)$", tmp[1]))) {
            tmp.ret   <- parse.func(fun.cache, tmp[1], idx.B00)
            tvar      <- tmp.ret$val
            idx.B00   <- tmp.ret$idx
            fun.cache <- tmp.ret$fun.cache
			idx.const.B00 <- rbind(idx.const.B00, c(i.x, match(tmp[2], a.ltname)))
          }
		  if (is.null(tvar) || is.na(tvar))
            stop(paste("The modifier [", tvar, "] does not recognized!"))

		  i.y <- match(tmp[2], a.ltname)
		  n.par <- tvar
		} else {
		  if (b.equalP) {
		    tvar <- idx.B00
			idx.B00 <- idx.B00 + 1		    
			idx.const.B00 <- rbind(idx.const.B00, c(i.x, match(j, a.ltname)))
		  } else
		    tvar <- n.param
		  i.y <- match(j, a.ltname)
		  n.par <- tvar
		}
		if (b.debug == TRUE)
			browser()
        r.B00[i.x, i.y] <- n.par
      }
    }
	# B00 check
	if (length(nzidx <- which(diag(r.B00) > 0)))
	  stop(paste("The following latent(s) have recursion!", paste(a.ltname[nzidx], collapse=" ")))
    r.B00 <- t(r.B00)
  }
#  browser()
  # Update B00 if high-order exist and 'reflective'
  if (length(a.hridx) > 0) {
    a.idx.B00 <- which(r.W002!=0, arr.ind=TRUE)
    n.nonho <- length(a.ltname)-length(a.hridx)
    for (j in 1:nrow(a.idx.B00)) {
      tmp <- a.idx.B00[j,]	
      if (attr(a.hridx[[tmp[2]]], "reflective") != 0 &&
	    r.B00[tmp[2]+n.nonho,tmp[1]] == 0) {
		if (a.holatent[tmp[2]] == 0) next()
		if (b.equalL == FALSE) val <- 99
		else {
		  val <- idx.B00
		  idx.B00 <- idx.B00 + 1
		}
        r.B00[tmp[2]+n.nonho,tmp[1]] <- val
      }
    }
  }
#  browser()
  if (b.debug == TRUE) {
    cat("Browser open by debug mode\n")
    browser()
  }
  
  a.remained <- which(nchar(a.fml) > 0)
  if (length(a.remained)) {
    stop(paste("The following lines were not parsed:", paste(a.remained, collapse=" ")))
  } else if (b.debug == TRUE)
    cat("## SUCCEEDED [", gsub("\\n", " ", s.fml), "]\n")
	
  # Re-arrange B00 if required
  if (b.equalP == TRUE) {
  	a.idx.B00 <- sort((idx.const.B00[,1]-1)*ncol(r.B00) + idx.const.B00[,2])
#	a.idx.B00 <- which(r.B00 != 0 & r.B00 != 99)
	if (length(a.idx.B00) > 0) {
	  i.idx.B00 <- min(r.B00[a.idx.B00])
	  for (i in a.idx.B00) {
	    r.B00[i] <- i.idx.B00
	    i.idx.B00 <- i.idx.B00 + 1
	  }
	}
  }
    
  # Which row is empty?
  i.val <- rep(1, nrow(r.W00))
  i.val[which(rowSums(r.W00) == 0)] <- 0
  i.cum <- cumsum(i.val)
  names(i.cum) <- 1:length(i.cum)
  #print(i.cum)
  m.data.ret <- m.data[,i.val != 0]
  #print(m.data.ret)
  #print(head(m.data.ret))
  for (i in names(a.ridx)) {
    a.ridx[[i]] <- i.cum[a.ridx[[i]]]
    names(a.ridx[[i]]) <- NULL
  }
  r.W00 <- r.W00[i.val != 0,,drop=F]
  r.C00 <- r.C00[,i.val != 0,drop=F]
  
  # Make sure Z0 as numeric
  m.data.ret <- matrix(as.numeric(as.matrix(m.data.ret)), ncol=sum(i.val!=0))
  colnames(m.data.ret) <- colnames(m.data)[i.val != 0]
#  for (i in 1:ncol(m.data.ret))
#    m.data.ret[,i] <- as.numeric(m.data.ret[,i])

  if (length(i.higher)) return(list(
    idxmanifest = a.ridx,
    idxlatent = a.hridx,
    loadtype = r.loadtype,
    loadtype2 = r.loadtype2,
    Z0 = m.data.ret,
    W00 = r.W00,
    W002 = r.W002,
    C00 = r.C00,
    B00 = r.B00
  )) else return(list(
    idxmanifest = a.ridx,
    loadtype = r.loadtype,
    Z0 = m.data.ret,
    W00 = r.W00,
    C00 = r.C00,
    B00 = r.B00
  ))
}

#############
#
# INTERNAL FUNCTION PART
#
#############

idxListLenNotEq <- function(v, l) which(unlist(lapply(v, length)) != l)
getFml <- function(v, d) {
  # Divide
  a.ret <- strsplit(v, d)
  # stop: invalid syntax (divided length by =~ is not 2)
  if (length(i.lenchk <- idxListLenNotEq(a.ret, 2)))
    stop(paste("The following equations have problem:\n",
      paste(paste(" ", i.lenchk), collapse="\n")
    ))
  a.ret
}

do.test <- function(s.name, b.expectErr=F, f.fun=NULL, ...) {
  if (!is.null(options("f.fun"))) f.fun <- options("f.fun")[[1]]
  if (!is.function(f.fun)) {
    print(f.fun)
    stop("ERROR")
  }
  ret <- try(f.fun(...), T)
  ret <- try(f.fun(...), T)
  b.fail <- class(ret) == 'try-error'
  if (b.expectErr == T) b.fail <- !b.fail
  if (b.fail) {
    cat("Test [", s.name, "] failed\n")
    print(ret)
  } else {
    if (b.expectErr == T) cat("Test [", s.name, "] ok (error raised)\n")
    else cat("Test [", s.name, "] ok\n")
  }
}

prt.latent <- function(lname, samplesizes, prt3, prt4, prt5, v) {
	nw <- ncol(samplesizes)
	if (nrow(v) > 1) {
	  cat(sprintf(prt3, ""))
	  for (i in 1:nrow(v)) cat(sprintf(prt5, paste("Group", i)))
	  cat("\n")
	}
	for (i in 1:ncol(v)) {
	  cat(sprintf(prt3, lname[i]))
	  for (j in 1:nrow(v)) cat(sprintf(prt4, round(v[j,i], 4)))
	  cat("\n")
	}
	cat("\n")
}

prt.msd <- function(gname, grp, n1, n2, samplesizes, prt3, prt4, prt5, prt6, v1, v2=NULL) {
	n.group <- ncol(samplesizes)
	
	if (all(v1==0)) {
		cat(sprintf(prt6, "NULL"), "\n\n")
		return()
	}

	# Get SD & remap
	nzidx <- which(v1 != 0, arr.ind=T)
	if (is.null(v2)) {
		tesd <- test.mat <- NULL
	} else {
		tesd <- apply(t(v2), 1, sd)
		tesd.mat <- matrix(0, nrow(v1), ncol(v1))
		for (i in 1:nrow(nzidx))
		  tesd.mat[nzidx[i,1],nzidx[i,2]] <- tesd[i]
	}

	pl <- 0
	for (j in 1:n.group) {
	  if (n.group > 1)
	    cat("  Group", j, "[", gname, "=", grp[j], "]:\n")
	  # Header
	  cat(sprintf(prt6, ""))
	  for (i in 1:length(n1)) cat(sprintf(prt5, n1[i]))
	  cat("\n")

	  for (i in 1:length(n2)) {
		cat(sprintf(prt6, ifelse(is.null(tesd),n2[i],"")))

	    # Print mean first
		for (j in 1:length(n1)) {
		  v <- v1[j+pl,i]
		  if (v == 0)
		    cat(sprintf(prt5, ""))
		  else
		    cat(sprintf(prt4, round(v, 4)))
		}
		cat("\n")

		# Print sd then
		if (!is.null(tesd)) {
			cat(sprintf(prt6, n2[i]))
			for (j in 1:length(n1)) {
			  v <- v1[j+pl,i]

			  if (v == 0)
				cat(sprintf(prt5, ""))
			  else
				cat(sprintf(prt5, sprintf("(%g)", round(tesd.mat[j+pl,i], 4))))
			}
			cat("\n")
		}
	  }
	  pl <- pl + length(n1)
	  cat("\n")
	}

}

qualmeasures <- function(object) {
	cat("  Type of indicators per latent variable:\n")
	cat("   (0=formative, 1=reflective)\n\n")

	n.prtcollen <- 9
	nv <- max(nchar(object$lname)) + 2
	prt3 <- paste("%", nv, "s", sep="")
	prt4 <- paste("  %", nv, "g", sep="")
	prt5 <- paste("  %", nv, "s", sep="")
    
    prt.lname <- object$lname
	prt.loadtype <- object$loadtype
    while (length(prt.lname)) {
      lenrow <- 80
      for (i in 1:length(prt.lname)) {
        cur <- sprintf(prt5, prt.lname[i])
        if (nchar(cur)>lenrow) {
          i <- i - 1
          break()
        }
	    cat(cur)
        lenrow <- lenrow - nchar(cur)
      }
	  cat("\n")
	  for (j in 1:i)
	    cat(sprintf(prt4, prt.loadtype[j]))
	  cat("\n")
      prt.lname <- prt.lname[-(1:i)]
      prt.loadtype <- prt.loadtype[-(1:i)]
    }
    cat("\n")

	cat("  Cronbach's alpha:\n")
	prt.latent(object$lname, object$samplesizes, prt3, prt4, prt5, object$Alpha)

	cat("  Dilon-Goldstein's rho:\n")
	prt.latent(object$lname, object$samplesizes, prt3, prt4, prt5, object$rho)

	cat("  Average Variance Extracted (AVE):\n")
	prt.latent(object$lname, object$samplesizes, prt3, prt4, prt5, object$AVE)

	cat("  Number of eigenvalues greater than one per block of indicators:\n")
	prt.latent(object$lname, object$samplesizes, prt3, prt4, prt5, object$Dimension)
}

latentmeasures <- function(object) {
	n.group <- ncol(object$samplesizes)
	n.prtcollen <- 9
	nv <- max(nchar(object$lname)) + 2
	prt3 <- paste("%", nv, "s", sep="")
	prt4 <- paste("  %", n.prtcollen, "g", sep="")
	prt5 <- paste("  %", n.prtcollen, "s", sep="")

	cat("  Means of latent variables:\n")
	prt.latent(object$lname, object$samplesizes, prt3, prt4, prt5, object$LV_MEAN)

	cat("  Variances of latent variables:\n")
	prt.latent(object$lname, object$samplesizes, prt3, prt4, prt5, object$LV_VAR)

	cat("  Correlations of latent variables:\n\n")
	pl <- 0
	for (j in 1:n.group) {
	  if (n.group > 1)
	    cat("  Group", j, "[", object$gname, "=", object$grp[j], "]:\n")

	  for (i in 1:length(object$lname)) {
	    cat(sprintf(prt3, object$lname[i]))
	    for (j in 1:i)
		  cat(sprintf(prt4, round(object$latentcorr[i+pl,j+pl], 4)))
		cat("\n")
	  }
	  pl <- pl + length(object$lname)
	  cat("\n")
	}
}

effectmeasures <- function(object) {
	n.prtcollen <- 9
	nv <- max(max(nchar(object$lname)) + 2, n.prtcollen)
	nv2 <- max(max(nchar(object$wname)) + 2, n.prtcollen)
	prt3 <- paste("%", nv, "s", sep="")
	prt4 <- paste("  %", nv, "g", sep="")
	prt5 <- paste("  %", nv, "s", sep="")
	prt6 <- paste("%", nv2, "s", sep="")

	if (is.null(object$MatTE_S)) {
		cat("Total effects of latent variables:\n\n")
		prt.msd(object$gname, object$grp, object$lname, object$lname, object$samplesizes,
			prt3, prt4, prt5, prt3, object$TE_S)

		cat("Indirect effects of latent variables:\n\n")
		prt.msd(object$gname, object$grp, object$lname, object$lname, object$samplesizes,
			prt3, prt4, prt5, prt3, object$ID_S)

		cat("Total effects of latent variables on indicators:\n\n")
		prt.msd(object$gname, object$grp, object$lname, object$wname, object$samplesizes,
			prt3, prt4, prt5, prt6, object$TE_M)

		cat("Indirect effects of latent variables on indicators:\n\n")
		prt.msd(object$gname, object$grp, object$lname, object$wname, object$samplesizes,
			prt3, prt4, prt5, prt6, object$ID_M)
	} else {
		cat("Total effects of latent variables(Std.Error):\n\n")
		prt.msd(object$gname, object$grp, object$lname, object$lname, object$samplesizes,
			prt3, prt4, prt5, prt3, object$TE_S, object$MatTE_S)

		cat("Indirect effects of latent variables(Std.Error):\n\n")
		prt.msd(object$gname, object$grp, object$lname, object$lname, object$samplesizes,
			prt3, prt4, prt5, prt3, object$ID_S, object$MatID_S)

		cat("Total effects of latent variables on indicators(Std.Error):\n\n")
		prt.msd(object$gname, object$grp, object$lname, object$wname, object$samplesizes,
			prt3, prt4, prt5, prt6, object$TE_M, object$MatTE_M)

		cat("Indirect effects of latent variables on indicators(Std.Error):\n\n")
		prt.msd(object$gname, object$grp, object$lname, object$wname, object$samplesizes,
			prt3, prt4, prt5, prt6, object$ID_M, object$MatID_M)
	}
}

fitmeasures <- function(object) {
	cat("  Number of parameters       ", object$NPAR, "\n")
	cat("  Number of bootstrap samples", object$nbt, "\n\n")
	
	lb <- object$lb
	ub <- object$ub

	vars <- c("FIT", "Adjusted FIT (AFIT)", "GFI",
		"Standardized Root Mean Square (SRMR)",
		"FIT_M", "FIT_S")
	n.prtcollen <- 9
	nv <- max(nchar(vars)) + 2
	if (!is.null(object$vec_FIT))
	  res1 <- rbind(
		round(
			cbind(object$FIT, as.matrix(apply(t(object$vec_FIT),1,sd)),
			as.matrix(t(object$sortFIT)[,lb]), as.matrix(t(object$sortFIT)[,ub])),4),
		round(
			cbind(object$AFIT, as.matrix(apply(t(object$vec_AFIT),1,sd)),
			as.matrix(t(object$sortAFIT)[,lb]), as.matrix(t(object$sortAFIT)[,ub])),4),
		round(
			cbind(object$GFI, as.matrix(apply(t(object$vec_GFI),1,sd)),
			as.matrix(t(object$sortGFI)[,lb]), as.matrix(t(object$sortGFI)[,ub])),4),
		round(
			cbind(object$SRMR, as.matrix(apply(t(object$vec_SRMR),1,sd)),
			as.matrix(t(object$sortSRMR)[,lb]), as.matrix(t(object$sortSRMR)[,ub])),4),
		round(
			cbind(object$FIT_M, as.matrix(apply(t(object$vec_FIT_m),1,sd)),
			as.matrix(t(object$sortFIT_m)[,lb]), as.matrix(t(object$sortFIT_m)[,ub])),4),
		round(
			cbind(object$FIT_S, as.matrix(apply(t(object$vec_FIT_s),1,sd)),
			as.matrix(t(object$sortFIT_s)[,lb]), as.matrix(t(object$sortFIT_s)[,ub])),4)
	) else
	  res1 <- rbind(
		round(cbind(object$FIT),4),
		round(cbind(object$AFIT),4),
		round(cbind(object$GFI),4),
		round(cbind(object$SRMR),4),
		round(cbind(object$FIT_M),4),
		round(cbind(object$FIT_S),4)
	)
	
	prt1 <- paste("%", nv, "s   %", n.prtcollen, "g   %", n.prtcollen, "g   %", n.prtcollen, "g   %", n.prtcollen, "g\n", sep="")
	prt2 <- paste("%", nv, "s   %", n.prtcollen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s\n", sep="")

	prt1s <- paste("%", nv, "s   %", n.prtcollen, "g\n", sep="")
	prt2s <- paste("%", nv, "s   %", n.prtcollen, "s\n", sep="")

	# Print header
	if (ncol(res1) == 4) {
		cat(sprintf(prt2, "", "Measure", "Std.Error", "95%CI_LB", "95%CI_UB"))
		# Print contents
		for (i in 1:length(vars))
		  cat(sprintf(prt1, vars[i], res1[i,1], res1[i,2], res1[i,3], res1[i,4]))
	} else {
		cat(sprintf(prt2s, "", "Measure"))
		# Print contents
		for (i in 1:length(vars))
		  cat(sprintf(prt1s, vars[i], res1[i,1]))
	}
}

summary.gesca <- function(object, ...) {
	n.group     <- ncol(object$samplesizes)     # Number of groups in the result
	n.rndigit   <- 4                            # Number of digits to preserve in real number
    n.prtcollen <- 11                           # Print length of column width
    n.prtvarlen <- max(nchar(object$wname))+2   # Print length of variable width
    n.prtvarlen2 <- max(nchar(object$lname))+2  # Print length of variable width
    n.prtrellen <- max(nchar(object$lname))*2+2 # Print length of variable width
    n.prtgrplen <- max(nchar(object$grp))       # Print length of group width
	
	wname <- object$wname
	lb    <- object$lb
	ub    <- object$ub
	gname <- object$gname
	grp   <- object$grp

    ##########################
    # Define printing formats
    ##########################

	prt1 <- paste("%", n.prtvarlen, "s   %", n.prtcollen, "g   %", n.prtcollen, "g   %", n.prtcollen, "g   %", n.prtcollen, "g\n", sep="")
	prt1s <- paste("%", n.prtvarlen, "s   %", n.prtcollen, "g\n", sep="")
	prt2 <- paste("%", n.prtvarlen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s\n", sep="")
	prt1x <- paste("%", n.prtrellen, "s   %", n.prtcollen, "g   %", n.prtcollen, "g   %", n.prtcollen, "g   %", n.prtcollen, "g\n", sep="")
	prt1xs <- paste("%", n.prtrellen, "s   %", n.prtcollen, "g\n", sep="")
	prt2x <- paste("%", n.prtrellen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s\n", sep="")
	prt2y <- paste("%", n.prtvarlen2, "s   %", n.prtcollen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s   %", n.prtcollen, "s\n", sep="")
	prt2ys <- paste("%", n.prtvarlen2, "s   %", n.prtcollen, "s\n", sep="")
	prt3 <- paste("%", n.prtvarlen2, "s", sep="")
	prt5 <- paste("  %", n.prtcollen, "s", sep="")
	prt4 <- paste("  %", n.prtcollen, "g", sep="")

    ##########################
    # (1) Print default information
    ##########################

	cat("The ALS algorithm converged in", object$niter,
		"iterations ( convergence criterion =", object$eps, ")\n\n")

	if (n.group > 1) {
	    cat("Number of observations per group [", gname, "]\n\n")
		for (i in 1:n.group) {
          n.space <- n.prtgrplen - nchar(grp[i])
		  cat("  Group ", i, " [ ", gname, " = ", grp[i], " ]",
			rep(" ", n.space), "  ", object$samplesizes[1,i], "\n", sep="")
        }
		cat("\n")
	} else
        cat("  Number of observations      ", object$samplesizes[1,1] , "\n")
	cat("  Number of parameters        ", object$NPAR, "\n")
	cat("  Number of bootstrap samples ", object$nbt, "\n\n")

    ##########################
    # (2) Print model fit measures
    ##########################

	cat("Model fit measures:\n\n")
	cat("  FIT                                 ", round(object$FIT, n.rndigit), "\n")
	cat("  Adjusted FIT (AFIT)                 ", round(object$AFIT, n.rndigit), "\n")
	cat("  GFI                                 ", round(object$GFI, n.rndigit), "\n")
	cat("  Standardized Root Mean Square (SRMR)", round(object$SRMR, n.rndigit), "\n")
	cat("  FIT_M                               ", round(object$FIT_M, n.rndigit), "\n")
	cat("  FIT_S                               ", round(object$FIT_S, n.rndigit), "\n\n")

    ##########################
    # (3) Print estimates of weights
    ##########################
	sub.wname <- wname[which(object$W00 != 0, arr.ind=T)[,1]]
    if (is.null(object$mW)) {
		if (is.null(object$W002))
			cat("Estimates of Weights:\n\n")
		else
			cat("Estimates of First-order Weights:\n\n")
        if (is.null(object$MatW1))
            # (10)
            res10 <- round(cbind(object$mW1),4)
        else
            # (10)
            res10 <- round(cbind(object$mW1,
				as.matrix(apply(t(object$MatW1),1,sd)),
				as.matrix(t(object$sortw1)[,lb]), as.matrix(t(object$sortw1)[,ub])),4)
			
        pl <- 0
        for (j in 1:n.group) {
          if (n.group > 1)
            cat("  Group", j, "[", gname, "=", grp[j], "]:\n")
          if (ncol(res10) == 4)
            cat(sprintf(prt2, "", "Estimate", "Std.Error", "95%CI_LB", "95%CI_UB"))
          else
            cat(sprintf(prt2, "", "Estimate", "", "", ""))
          
          for (i in 1:length(sub.wname))
            if (ncol(res10) == 4) cat(sprintf(prt1, sub.wname[i],
              res10[i+pl,1], res10[i+pl,2],
              res10[i+pl,3], res10[i+pl,4]))
            else cat(sprintf(prt1s, sub.wname[i],
              res10[i+pl,1]))
          pl <- pl + length(sub.wname)
          cat("\n")
        }
		# Print second-order weights
		if (!is.null(object$mW2)) {
			cat("Estimates of Second-order Weights:\n\n")
            # (10)
			if (is.null(object$MatW1))
				res10.2 <- round(cbind(object$mW2),4)
			else
				res10.2 <- round(cbind(object$mW2,
					as.matrix(apply(t(object$MatW2),1,sd)),
					as.matrix(t(object$sortw2)[,lb]), as.matrix(t(object$sortw2)[,ub])),4)

			pl <- 0
			for (j in 1:n.group) {
			  if (n.group > 1)
				cat("  Group", j, "[", gname, "=", grp[j], "]:\n")
			  if (ncol(res10.2) == 4)
				cat(sprintf(prt2y, "", "Estimate", "Std.Error", "95%CI_LB", "95%CI_UB"))
			  else
				cat(sprintf(prt2y, "", "Estimate", "", "", ""))
				
			  lname.a <- na.omit(object$lname[rowSums(object$W002)!=0][1:nrow(object$W002)])
			  
			  for (i in 1:length(lname.a))
				if (ncol(res10.2) == 4) cat(sprintf(prt2y, lname.a[i],
				  res10.2[i+pl,1], res10.2[i+pl,2],
				  res10.2[i+pl,3], res10.2[i+pl,4]))
				else cat(sprintf(prt2ys, lname.a[i],
				  res10.2[i+pl,1]))
			  pl <- pl + length(lname.a)
			  cat("\n")
			}     
		}
    } else {
        cat("Estimates of Weights:\n\n")
		if (is.null(object$MatW))
			# (10)
			res10 <- round(cbind(object$mW),4)
		else
			# (10)
			res10 <- round(cbind(object$mW,
				as.matrix(apply(t(object$MatW),1,sd)),
				as.matrix(t(object$sortw)[,lb]), as.matrix(t(object$sortw)[,ub])),4)
        pl <- 0
        for (j in 1:n.group) {
          if (n.group > 1)
            cat("  Group", j, "[", gname, "=", grp[j], "]:\n")
		  if (ncol(res10) == 4)
            cat(sprintf(prt2, "", "Estimate", "Std.Error", "95%CI_LB", "95%CI_UB"))
		  else
		    cat(sprintf(prt2, "", "Estimate", "", "", ""))
          for (i in 1:length(sub.wname))
            if (ncol(res10) == 4) cat(sprintf(prt1, sub.wname[i],
              res10[i+pl,1], res10[i+pl,2],
              res10[i+pl,3], res10[i+pl,4])) 
			else cat(sprintf(prt1s, sub.wname[i],
              res10[i+pl,1])) 
          pl <- pl + length(sub.wname)
          cat("\n")
        }
    }

    ##########################
    # (4) Print estimates of loadings
    ##########################

    cat("Estimates of Loadings:\n\n")
    sub.wname <- wname[which(t(object$C00) != 0,arr.ind=TRUE)[,1]]
    if (length(sub.wname) > 0) {
        nv <- max(nchar(sub.wname))
        if (is.null(object$Matload))
          # (11)
          res11 <- round(cbind(object$mC), 4)
        else
          # (11)
          res11 <- round(cbind(object$mC,
			as.matrix(apply(t(object$Matload),1,sd)),
			as.matrix(t(object$sortload)[,lb]), as.matrix(t(object$sortload)[,ub])),4)
        pl <- 0
        for (j in 1:n.group) {
          if (n.group > 1)
            cat("  Group", j, "[", gname, "=", grp[j], "]:\n")
          if (ncol(res11) == 4)
            cat(sprintf(prt2, "", "Estimate", "Std.Error", "95%CI_LB", "95%CI_UB"))
          else
            cat(sprintf(prt2, "", "Estimate", "", "", ""))
          for (i in 1:length(sub.wname)) {
            if (ncol(res11) == 4)
              cat(sprintf(prt1, sub.wname[i],
                res11[i+pl,1], res11[i+pl,2],
                res11[i+pl,3], res11[i+pl,4]))
            else
              cat(sprintf(prt1s, sub.wname[i],
                res11[i+pl,1]))
          }
          pl <- pl + length(sub.wname)
          cat("\n")
        }
    } else
      cat("  NULL\n\n")

	cat("Estimates of Path Coefficients:\n\n")
	# (12)
    idx  <- which(t(object$B00) != 0, arr.ind=TRUE)
    if (nrow(idx)) {
	  idx  <- cbind(idx[,2], idx[,1])
	  if (is.null(object$Matbeta))
	    res12 <- round(cbind(object$mB),4)
	  else
        res12 <- round(cbind(object$mB,
			as.matrix(apply(t(object$Matbeta),1,sd)),
			as.matrix(t(object$sortbeta)[,lb]), as.matrix(t(object$sortbeta)[,ub])),4)
	  pl <- 0
	  for (j in 1:n.group) {
  	  if (n.group > 1)
      cat("  Group", j, "[", gname, "=", grp[j], "]:\n")
      if (n.prtrellen > 20) {
	    if (ncol(res12) == 4)
          cat(sprintf(prt2y, "", "Estimate", "Std.Error", "95%CI_LB", "95%CI_UB"))
		else
		  cat(sprintf(prt2y, "", "Estimate", "", "", ""))
        for (i in 1:nrow(idx)) {
            cat(sprintf(prt3, paste(object$lname[idx[i,2]], "~", sep="")), "\n")
			if (ncol(res12) == 4) cat(sprintf(prt2y, object$lname[idx[i,1]],
                res12[i+pl,1], res12[i+pl,2],
                res12[i+pl,3], res12[i+pl,4]))
			else cat(sprintf(prt2ys, object$lname[idx[i,1]],
                res12[i+pl,1]))
        }
      } else {
	    if (ncol(res12) == 4) cat(sprintf(prt2x, "", "Estimate", "Std.Error", "95%CI_LB", "95%CI_UB"))
		else cat(sprintf(prt2x, "", "Estimate", "", "", ""))
        for (i in 1:nrow(idx)) {
		  if (ncol(res12) == 4) cat(sprintf(prt1x,
            paste(object$lname[idx[i,2]], "~", object$lname[idx[i,1]], sep=""),
              res12[i+pl,1], res12[i+pl,2],
              res12[i+pl,3], res12[i+pl,4]))
		  else cat(sprintf(prt1xs,
            paste(object$lname[idx[i,2]], "~", object$lname[idx[i,1]], sep=""),
              res12[i+pl,1]))
        }
      }
	  pl <- pl + nrow(idx)
	  cat("\n")
	}
	}

	cat("R-squared Values of Endogenous Latent Variables:\n\n")
	prt.latent(object$lname, object$samplesizes, prt3, prt4, prt5, object$R2)
}
