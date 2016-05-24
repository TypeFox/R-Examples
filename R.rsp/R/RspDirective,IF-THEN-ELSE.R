###########################################################################/**
# @RdocClass RspIfDirective
# @alias RspElseDirective
# @alias RspEndifDirective
# @alias RspIfdefDirective
# @alias RspIfndefDirective
# @alias RspIfeqDirective
# @alias RspIfneqDirective
#
# @title "The RspIfDirective class"
#
# \description{
#  @classhierarchy
#
#  An RspIfDirective is an @see "RspDirective" that will include or
#  exclude all @see "RspConstruct":s until the next @see "RspEndifDirective"
#  based on the preprocessing value of the particular if clause.
#  Inclusion/exclusion can be reversed via an @see "RspElseDirective".
# }
#
# @synopsis
#
# \arguments{
#   \item{value}{A @character string.}
#   \item{...}{Arguments passed to the constructor of @see "RspDirective".}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setConstructorS3("RspIfDirective", function(value="if", ...) {
  this <- extend(RspDirective(value, ...), "RspIfDirective")
  if (!missing(value)) {
    requireAttributes(this, c("test"));

    # Test aliases
    test <- getAttribute(this, "test");
    map <- c(
      "==" = "equal-to",
      "!=" = "not-equal-to",
      ">"  = "greater-than",
      ">=" = "greater-than-or-equal-to",
      "<"  = "less-than",
      "<=" = "less-than-or-equal-to"
    );
    testA <- map[test];
    if (!is.na(testA)) {
      this <- setAttribute(this, "test", testA);
    }
  }
  this;
})

setConstructorS3("RspElseDirective", function(value="else", ...) {
  extend(RspDirective(value, ...), "RspElseDirective")
})

setConstructorS3("RspEndifDirective", function(value="endif", ...) {
  extend(RspDirective(value, ...), "RspEndifDirective")
})


# Alias: <%@ifeq ...%> => <%@if test="equals" ...%>
setConstructorS3("RspIfeqDirective", function(value="if", ...) {
  extend(RspIfDirective(value, test="equal-to", ...), "RspIfeqDirective");
})

# Alias: <%@ifneq ...%> => <%@if test="equals" negate="TRUE", ...%>
setConstructorS3("RspIfneqDirective", function(value="if", ...) {
  this <- extend(RspIfeqDirective(value, ...), "RspIfneqDirective");
  # Negate the 'test' result
  negate <- !getAttribute(this, "negate", FALSE);
  this <- setAttribute(this, "negate", negate);
  this;
})


# Alias: <%@ifdef ...%> => <%@if test="exists" ...%>
setConstructorS3("RspIfdefDirective", function(value="if", ...) {
  extend(RspIfDirective(value, test="exists", ...), "RspIfeqDirective");
})

# Alias: <%@ifndef ...%> => <%@if test="exists" negate="TRUE", ...%>
setConstructorS3("RspIfndefDirective", function(value="if", ...) {
  this <- extend(RspIfdefDirective(value, ...), "RspIfneqDirective");
  # Negate the 'test' result
  negate <- !getAttribute(this, "negate", FALSE);
  this <- setAttribute(this, "negate", negate);
  this;
})


##############################################################################
# HISTORY:
# 2013-10-18
# o Added aliases <%@ifdef ...%>/<%@ifndef ...%> for <@if test="exists" ...%>.
# 2013-03-17
# o Now <%@ifeq ...%>/<%@ifneq ...%> is an alias for <@if test="equals" ...%>.
# o Moved all RSP if-then-else directives to one file.
# 2013-02-18
# o Added RspIfeqDirective, RspElseDirective, and RspEndifDirective.
##############################################################################
