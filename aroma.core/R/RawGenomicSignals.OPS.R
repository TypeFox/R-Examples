setMethodS3("subtractBy", "RawGenomicSignals", function(this, ...) {
  applyBinaryOperator(this, ..., FUN=get("-", mode="function"));
})

setMethodS3("addBy", "RawGenomicSignals", function(this, ...) {
  applyBinaryOperator(this, ..., FUN=get("+", mode="function"));
})

setMethodS3("divideBy", "RawGenomicSignals", function(this, ...) {
  applyBinaryOperator(this, ..., FUN=get("/", mode="function"));
})

setMethodS3("multiplyBy", "RawGenomicSignals", function(this, ...) {
  applyBinaryOperator(this, ..., FUN=get("*", mode="function"));
})


setMethodS3("applyBinaryOperator", "RawGenomicSignals", function(this, other, fields=NULL, FUN, ..., sort=FALSE) {
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, class(this)[1]);

  # Argument 'fields':
  if (is.null(fields)) {
    fields <- getColumnNames(this);
  }

  # Argument 'FUN':
  if (!is.function(FUN)) {
    throw("Argument 'FUN' is not a function: ", mode(FUN)[1]);
  }
 

  nbrOfLoci <- nbrOfLoci(this);
  if (nbrOfLoci(other) != nbrOfLoci) {
    throw("The number of loci in argument 'other' does not match the number of loci in this object: ", nbrOfLoci(other), " != ", nbrOfLoci);
  }

  fieldsOther <- getColumnNames(other);
  fields <- intersect(fields, fieldsOther);

  res <- this;

  # Sort by genomic position?
  if (sort) {
    # sort() returns a sorted clone():d object. /HB 2012-03-01
    res <- sort(res);
    other <- sort(other);
  }

  # Has genomic locations?
  if (is.element("x", fields)) {
    # Assert that positions are the same
    if (!all.equal(this$x, other$x)) {
      throw("Cannot subtract argument 'other' from this object, because their genomic locations do not match.");
    }

    # Keep positions
    fields <- setdiff(fields, "x");
  }


  for (field in fields) {
    delta <- FUN(res[[field]], other[[field]]);
    res[[field]] <- delta;
  }

  res;  
}, protected=TRUE)



setMethodS3("+", "RawGenomicSignals", function(e1, e2) {
  # To please R CMD check 
  this <- e1;
  other <- e2;

  addBy(this, other);
}, appendVarArgs=FALSE, validators=NULL);

setMethodS3("-", "RawGenomicSignals", function(e1, e2) {
  # To please R CMD check 
  this <- e1;
  other <- e2;

  subtractBy(this, other);
}, appendVarArgs=FALSE, validators=NULL);

setMethodS3("*", "RawGenomicSignals", function(e1, e2) {
  # To please R CMD check 
  this <- e1;
  value <- e2;

  # Swap 'this' and 'value'?
  if (inherits(value, "RawGenomicSignals")) {
    tmp <- this;
    this <- value;
    value <- tmp;
  }

  value <- Arguments$getDouble(value);

  fields <- getColumnNames(this, translate=FALSE);
  fields <- setdiff(fields, "x");

  res <- clone(this);
  for (field in fields) {
    res[[field]] <- value * res[[field]];
  } 

  res;
}, appendVarArgs=FALSE, validators=NULL);


############################################################################
# HISTORY:
# 2011-12-10
# o ROBUSTNESS: Turned of RCC validation for "+", "-" and "*" methods.
# 2010-09-11
# o Added basic support for operators +, - and * to RawGenomicSignals.
# 2009-05-10
# o Added {add,subtract,multiply,divide}By() to RawGenomicSignals.
# o Added applyBinaryOperator().
# o Created.
############################################################################
