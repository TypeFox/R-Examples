# Trotter

# Helpers
# Returns the kth permutation of elements
perm.worker <- function(k, elements) {
  n <- length(elements)
  if (n == 1) elements
  else {
    group <- as.integer((k - 1) / n)
    item <- (k - 1) %% n
    position <- if (group %% 2 == 0) n - item - 1 else item
    append(
      perm.worker(
        group + 1,
        elements[1:(n - 1)]
      ),
      elements[n],
      after = position
    )
  }
}

# Returns the kth r-combination of elements 
combination <- function(k, r, elements) {
  n <- length(elements)
  position <- 0
  d <- choose(n - position - 1, r - 1)
  
  while ((k - 1) >= d) {
    k <- k - d
    position <- position + 1
    d <- choose(n - position - 1, r - 1)
  }
  
  if (r <= 1) elements[position + 1]
  else {
    right.tail <- elements[(position + 2):length(elements)]
    c(elements[position + 1], combination(k, r - 1, right.tail))
  }
}

# Returns the kth r-selection of elements
selection <- function(k, r, elements) {
  n <- length(elements)
  position <- 0
  d <- choose(n + r - position - 2, r - 1)
  
  while ((k - 1) >= d) {
    k <- k - d
    position <- position + 1
    d <- choose(n + r - position - 2, r - 1)
  }
  
  if (r <= 1) elements[position + 1] ###
  else {
    tail <- elements[(position + 1):length(elements)]
    c(elements[position + 1], selection(k, r - 1, tail))
  }
}

# Returns the kth r-permutation of elements
permutation <- function(k, r, elements) {
  n <- length(elements)
  f <- factorial(r)
  group <- as.integer((k - 1) / f)
  item <- (k - 1) %% f
  
  comb <- combination(group + 1, r, elements)
  perm.worker(item + 1, comb)
}

# Returns the kth r-amalgam of elements
amalgam <- function(k, r, elements) {
  k <- k - 1
  sapply(
    1:r,
    function (i) {
      p <- length(elements) ^ (r - i)
      index <- as.integer(k / p)
      k <<- k %% p
      elements[index + 1]
    }
  )
}

# Returns the kth subset of elements
k.subset <- function(k, elements) {
  r <- c()
  for (i in 0:(length(elements) - 1)) 
    if (bitwAnd(k - 1, 2 ^ i) != 0) r <- c(r, elements[i + 1])
  r
}

# Index checks and adjustments
index.check <- function(i, n) {
  if (missing(i) || length(i) == 0) {
    cat("Warning: Missing an index. First combination returned.\n")
    i = 1
  }
  
  if (!is.numeric(i)) {
    cat("Warning: Numerical index expected. Index 1 used.\n")
    i = 1
  }
  
  if (length(i) > 1) {
    cat("Warning: A single index expected. Only the first element used.\n")
    i <- i[1]
  }
  
  if (i < 1 || i > n) {
    cat("Warning: Index out of bounds. Wrap-around used.\n")
    while (i < 1) i <- i + n
    i <- i %% (n + 1)
  }
  
  i
}

# Permutations Pseudo-Vector
# Permutations Pseudo-vector

setClass(
  Class = "PPV",
  representation(k = "numeric", items = "vector"),
)

setMethod(
  f = "show",
  signature = "PPV",
  definition = function(object) cat(
    sprintf(
      paste(
        "Instance of class PPV\n",
        "Pseudo-vector containing %s %s-permutations\n",
        "of items taken from the list:\n[%s]",
        collapse = ""
      ),
      length(object),
      object@k,
      paste(object@items, collapse = ", ")
    )
  )
)

#' Permutations Pseudo-Vector Length
#' @description
#' Get the length of a \code{PPV} instance.
#' @param x an instance of \code{PPV}
#' @return the number of permutations in pseudo-vector \code{x}
#' @details
#' Since \code{x} contains all the \code{k}-permutations of objects in vector
#' \code{items}, \code{length(x)} will return 
#' \code{choose(length(items), k) * factorial(k)}.
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 

setMethod(
  f = "length",
  signature = "PPV",
  definition = function(x) choose(length(x@items), x@k) * factorial(x@k)
)

#' Retrieve a Permutation by Index
#' @description
#' Access a permutation stored in a \code{PPV} instance by index.
#' @param x an instance of \code{PPV}.
#' @param i an index specifying the position of the sought permutation.
#' @param j not used.
#' @param drop not used.
#' @return the permutation located at position \code{i} in pseudo-vector \code{x}
#' @details
#' The permutation at index \code{i} of pseudo-vector \code{x} is not actually 
#' stored in memory but calculated as needed. The extract method is used solely
#' for interpretation.
#' 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 

setMethod(
  f = "[",
  signature = "PPV",
  definition = function(x, i, j, drop) {
    i <- index.check(i, length(x))
    if (!missing(j)) cat("Warning: Only first coordinate used.\n")
    permutation(i, x@k, x@items)
  }
)

# Export #######################################################################

#' Permutations Pseudo-Vector Constructor 
#' @description
#' The \code{PPV} class defines a pseudo-vector containing all 
#' the \code{k}-permutations of the objects stored
#' in \code{items}. The function \code{ppv} is a constructor for this class.
#' @aliases
#' permutation
#' permutations
#' @param k the number of objects taken at a time.
#' @param items a vector of objects to be permuted.
#' @return an instance of \code{PPV}.
#' @author Richard Ambler
#' @details
#' The arrangement of permutations is similar, but in many cases not identical, 
#' to that obtained from the
#' Steinhaus-Johnson-Trotter algorithm (see references).
#' @examples
#' # create a pseudo-vector of 5-permutations from the first 10 letters
#' p <- ppv(5, letters[1:10])
#' # generate a description
#' print(p)
#' # compatable with length
#' length(p)
#' # inspect a few of the permutations "stored" in p
#' p[1]
#' p[1000]
#' p[30240]
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 
#' @references
#' Steinhaus-Johnson-Trotter algorithm. (2014, April 29).
#' In \emph{Wikipedia, The Free Encyclopedia}.
#' Retrieved 13:24, September 5, 2014
#' @export
#' @import methods

ppv <- function(k, items) new(
  Class = "PPV", 
  k = k, 
  items = items
)


# Combinations Pseudo-Vector
# Combinations Pseudo-vector ###################################################

setClass(
  Class = "CPV",
  representation(k = "numeric", items = "vector"),
)

setMethod(
  f = "show",
  signature = "CPV",
  definition = function(object) cat(
    sprintf(
      paste(
        "Instance of class CPV\n",
        "Pseudo-vector containing %s %s-combinations\n",
        "of items taken from the list:\n[%s]",
        collapse = ""
      ),
      length(object),
      object@k,
      paste(object@items, collapse = ", ")
    )
  )
)

#' Combinations Pseudo-Vector Length
#' @description
#' Get the length of a \code{CPV} instance.
#' @param x an instance of \code{CPV}
#' @return the number of combinations in pseudo-vector \code{x}
#' @details
#' Since \code{x} contains all the \code{k}-combinations of objects in vector
#' \code{items}, \code{length(x)} will return \code{choose(length(items), k)}.
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 

setMethod(
  f = "length",
  signature = "CPV",
  definition = function(x) choose(length(x@items), x@k)
)

#' Retrieve a Combination by Index
#' @description
#' Access a combination stored in a \code{CPV} instance by index.
#' @param x an instance of \code{CPV}.
#' @param i an index specifying the position of the sought combination.
#' @param j not used.
#' @param drop not used.
#' @return the combination located at position \code{i} in pseudo-vector \code{x}
#' @details
#' The combination at index \code{i} of pseudo-vector \code{x} is not actually 
#' stored in memory but calculated as needed. The extract method is used solely
#' for interpretation.
#' 
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 

setMethod(
  f = "[",
  signature = "CPV",
  definition = function(x, i, j, drop) {
    i <- index.check(i, length(x))
    if (!missing(j)) cat("Warning: Only first coordinate used.\n")
    combination(i, x@k, x@items)
  }
)

# Export #######################################################################

#' Combinations Pseudo-Vector Constructor 
#' @description
#' The \code{CPV} class defines a pseudo-vector containing all 
#' the arranged \code{k}-combinations of the objects stored
#' in \code{items}. The function \code{cpv} is a constructor for this class.
#' @aliases
#' combination
#' combinations
#' @param k the number of objects taken at a time.
#' @param items a vector of objects to be combined.
#' @return an instance of \code{CPV}.
#' @author Richard Ambler
#' @details
#' The combinations are arranged according to the order in which the objects
#' appear in \code{items}. Combinations containing the first object in 
#' \code{items} are followed by combinations that contain the second object
#' but not the first, which are followed by combinations that contain the third
#' but neither the first or the second, etc.
#' @examples
#' # create a pseudo-vector of 10-combinations from the first 15 letters
#' c <- cpv(10, letters[1:15])
#' # generate a description
#' print(c)
#' # compatable with length
#' length(c)
#' # inspect a few of the combinations "stored" in c
#' c[1]
#' c[1000]
#' c[3003]
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 
#' @references
#' Steinhaus-Johnson-Trotter algorithm. (2014, April 29).
#' In \emph{Wikipedia, The Free Encyclopedia}.
#' Retrieved 13:24, September 5, 2014
#' @export
#' @import methods

cpv <- function(k, items) new(
  Class = "CPV", 
  k = k, 
  items = items
)

# Selections Pseudo-Vector
# Selections Pseudo-vector ###################################################

setClass(
  Class = "SPV",
  representation(k = "numeric", items = "vector"),
)

setMethod(
  f = "show",
  signature = "SPV",
  definition = function(object) cat(
    sprintf(
      paste(
        "Instance of class SPV\n",
        "Pseudo-vector containing %s %s-selections (combinations with replacement)\n",
        "of items taken from the list:\n[%s]",
        collapse = ""
      ),
      length(object),
      object@k,
      paste(object@items, collapse = ", ")
    )
  )
)

#' Selections Pseudo-Vector Length
#' @description
#' Get the length of an \code{SPV} instance.
#' @param x an instance of \code{SPV}
#' @return the number of selections (combinations with replacement) in pseudo-vector \code{x}
#' @details
#' Since \code{x} contains all the \code{k}-selections of objects in vector
#' \code{items}, \code{length(x)} will return \code{choose(length(items) + k - 1, k)}.
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}}  
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 

setMethod(
  f = "length",
  signature = "SPV",
  definition = function(x) choose(length(x@items) + x@k - 1, x@k)
)

#' Retrieve a Selection by Index
#' @description
#' Access a selection (combination with replacement) stored in an \code{SPV} instance by index.
#' @param x an instance of \code{SPV}.
#' @param i an index specifying the position of the sought selection.
#' @param j not used.
#' @param drop not used.
#' @return the selection located at position \code{i} in pseudo-vector \code{x}
#' @details
#' The selection at index \code{i} of pseudo-vector \code{x} is not actually 
#' stored in memory but calculated as needed. The extract method is used solely
#' for interpretation.
#' 
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 

setMethod(
  f = "[",
  signature = "SPV",
  definition = function(x, i, j, drop) {
    i <- index.check(i, length(x))
    if (!missing(j)) cat("Warning: Only first coordinate used.\n")
    selection(i, x@k, x@items)
  }
)

# Export #######################################################################

#' Selections Pseudo-Vector Constructor
#' @description
#' The \code{SPV} class defines a pseudo-vector containing all 
#' the arranged \code{k}-selections (combinations with replacement) of the objects stored
#' in \code{items}. The function \code{spv} is a constructor for this class.
#' @aliases
#' selection
#' selections
#' @param k the number of objects taken at a time.
#' @param items a vector of objects to be selected.
#' @return an instance of \code{SPV}.
#' @author Richard Ambler
#' @details
#' The selections are arranged according to the order in which the objects
#' appear in \code{items}. The arrangement is very similar to the arrangement
#' of combinations (see \link{cpv}) except that objects may be repeatedly selected.
#' 
#' @examples
#' # create a pseudo-vector of 10-selections from the first 15 letters
#' s <- spv(10, letters[1:15])
#' # generate a description
#' print(s)
#' # compatable with length
#' length(s)
#' # inspect a few of the combinations "stored" in s
#' s[1]
#' s[1000]
#' s[1961256]
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 
#' @references
#' Steinhaus-Johnson-Trotter algorithm. (2014, April 29).
#' In \emph{Wikipedia, The Free Encyclopedia}.
#' Retrieved 13:24, September 5, 2014
#' @export
#' @import methods

spv <- function(k, items) new(
  Class = "SPV", 
  k = k, 
  items = items
)

# Amalgams Pseudo-Vector
# Amalgams Pseudo-vector ###################################################

setClass(
  Class = "APV",
  representation(k = "numeric", items = "vector"),
)

setMethod(
  f = "show",
  signature = "APV",
  definition = function(object) cat(
    sprintf(
      paste(
        "Instance of class APV\n",
        "Pseudo-vector containing %s %s-amalgams (permutations with replacement)\n",
        "of items taken from the list:\n[%s]",
        collapse = ""
      ),
      length(object),
      object@k,
      paste(object@items, collapse = ", ")
    )
  )
)

#' Amalgams Pseudo-Vector Length
#' @description
#' Get the length of an \code{APV} instance.
#' @param x an instance of \code{APV}
#' @return the number of amalgams (permutations with replacement) in pseudo-vector \code{x}
#' @details
#' Since \code{x} contains all the \code{k}-amalgams of objects in vector
#' \code{items}, \code{length(x)} will return \code{length(items) ^ k)}.
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 

setMethod(
  f = "length",
  signature = "APV",
  definition = function(x) length(x@items) ^ x@k 
)

#' Retrieve an Amalgam by Index
#' @description
#' Access an amalgam (permutation with replacement) stored in an \code{APV} instance by index.
#' @param x an instance of \code{APV}.
#' @param i an index specifying the position of the sought amalgam
#' @param j not used.
#' @param drop not used.
#' @return the amalgam located at position \code{i} in pseudo-vector \code{x}
#' @details
#' The amalgam at index \code{i} of pseudo-vector \code{x} is not actually 
#' stored in memory but calculated as needed. The extract method is used solely
#' for interpretation.
#' 
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 

setMethod(
  f = "[",
  signature = "APV",
  definition = function(x, i, j, drop) {
    i <- index.check(i, length(x))
    if (!missing(j)) cat("Warning: Only first coordinate used.\n")
    amalgam(i, x@k, x@items)
  }
)

# Export #######################################################################

#' Amalgams Pseudo-Vector Constructor
#' @description
#' The \code{APV} class defines a pseudo-vector containing all 
#' the arranged \code{k}-amalgams (permutations with replacement) of the objects stored
#' in \code{items}. The function \code{apv} is a constructor for this class.
#' @aliases
#' amalgam
#' amalgams
#' @param k the number of objects taken at a time.
#' @param items a vector of objects to be amalgamated.
#' @return an instance of \code{APV}.
#' @author Richard Ambler
#' @details
#' The amalgams are arranged according to the order in which the objects
#' appear in \code{items}. The arrangement is very similar to that used by the \code{PPV} class
#' (see \link{ppv}) except that objects are replaced during permutation creation.
#' 
#' @examples
#' # create a pseudo-vector of 10-amalgams from the first 15 letters
#' a <- apv(10, letters[1:15])
#' # generate a description
#' print(a)
#' # compatable with length
#' length(a)
#' # inspect a few of the combinations "stored" in a
#' a[1]
#' a[1000000]
#' a[576650390625]
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}}
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 
#' @seealso Subsets Pseudo-Vector \code{\link{sspv}} 
#' @references
#' Steinhaus-Johnson-Trotter algorithm. (2014, April 29).
#' In \emph{Wikipedia, The Free Encyclopedia}.
#' Retrieved 13:24, September 5, 2014
#' @export
#' @import methods

apv <- function(k, items) new(
  Class = "APV", 
  k = k, 
  items = items
)

# Subsets Pseudo-Vector
# Subsets Pseudo-vector ###################################################

setClass(
  Class = "SSPV",
  representation(items = "vector"),
)

setMethod(
  f = "show",
  signature = "SSPV",
  definition = function(object) cat(
    sprintf(
      paste(
        "Instance of class SSPV\n",
        "Pseudo-vector containing the %s subsets\n",
        "of items taken from the list:\n[%s]",
        collapse = ""
      ),
      length(object),
      paste(object@items, collapse = ", ")
    )
  )
)

#' Subsets Pseudo-Vector Length
#' @description
#' Get the length of an \code{SSPV} instance.
#' @param x an instance of \code{SSPV}
#' @return the number of subsets in pseudo-vector \code{x}
#' @details
#' Since \code{x} contains all the subsets of objects in vector
#' \code{items}, \code{length(x)} will return \code{2 ^ length(items)}.
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}} 


setMethod(
  f = "length",
  signature = "SSPV",
  definition = function(x) 2 ^ length(x@items) 
)

#' Retrieve a Subset by Index
#' @description
#' Access asubset stored in an \code{SSPV} instance by index.
#' @param x an instance of \code{SSPV}.
#' @param i an index specifying the position of the sought amalgam
#' @param j not used.
#' @param drop not used.
#' @return the subset located at position \code{i} in pseudo-vector \code{x}
#' @details
#' The subset at index \code{i} of pseudo-vector \code{x} is not actually 
#' stored in memory but calculated as needed. The extract method is used solely
#' for interpretation.
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}}

setMethod(
  f = "[",
  signature = "SSPV",
  definition = function(x, i, j, drop) {
    i <- index.check(i, length(x))
    if (!missing(j)) cat("Warning: Only first coordinate used.\n")
    k.subset(i, x@items)
  }
)

# Export #######################################################################

#' Subsets Pseudo-Vector Constructor
#' @description
#' The \code{SSPV} class defines a pseudo-vector containing all 
#' the arranged subsets of the objects stored
#' in \code{items}. The function \code{sspv} is a constructor for this class.
#' @aliases
#' subsets
#' @param items a vector of objects to be subsetted.
#' @return an instance of \code{SSPV}.
#' @author Richard Ambler
#' @details
#' The subsets are arranged according to the order in which the objects
#' appear in \code{items}. The first subset, containing none of the objects, 
#' is \code{NULL}.
#' 
#' @examples
#' # create a pseudo-vector of subsets from the first 15 letters
#' ss <- sspv(letters[1:15])
#' # generate a description
#' print(ss)
#' # compatable with length
#' length(ss)
#' # inspect a few of the combinations "stored" in ss
#' ss[1]
#' ss[1000]
#' ss[32768]
#' @seealso Permutations Pseudo-Vector \code{\link{ppv}} 
#' @seealso Combinations Pseudo-Vector \code{\link{cpv}} 
#' @seealso Amalgams Pseudo-Vector \code{\link{apv}} 
#' @seealso Selections Pseudo-Vector \code{\link{spv}}
#' @export
#' @import methods

sspv <- function(items) new(
  Class = "SSPV",
  items = items
)