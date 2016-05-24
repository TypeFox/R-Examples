###########################################################################/**
# @RdocDefault setConstructorS3
#
# @title "Defines a class in S3/UseMethod style"
#
# \description{
#  Defines a class in R.oo/S3 style.
#  What this function currently does is simply creating a constructor
#  function for the class.
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{The name of the class.}
#   \item{definition}{The constructor defintion. \emph{Note: The constructor
#     must be able to be called with no arguments, i.e. use default values
#     for all arguments or make sure you use \code{missing()} or similar!}}
#   \item{static}{If @TRUE this class is defined to be static,
#      otherwise not. Currently this has no effect expect as an indicator.}
#   \item{abstract}{If @TRUE this class is defined to be abstract,
#      otherwise not. Currently this has no effect expect as an indicator.}
#   \item{private}{If @TRUE this class is defined to be private.}
#   \item{protected}{If @TRUE this class is defined to be protected.}
#   \item{export}{A @logical setting attribute \code{"export"}.}
#   \item{trial}{If @TRUE this class is defined to be a trial class,
#      otherwise not. A trial class is a class that is introduced to be
#      tried out and it might be modified, replaced or even removed in a
#      future release. Some people prefer to call trial versions, beta
#      version. Currently this has no effect expect as an indicator.}
#   \item{deprecated}{If @TRUE this class is defined to be deprecated,
#      otherwise not. Currently this has no effect expect as an indicator.}
#   \item{envir}{The environment for where the class (constructor function)
#      should be stored.}
#   \item{enforceRCC}{If @TRUE, only class names following the R Coding
#      Convention is accepted. If the RCC is violated an RccViolationException
#      is thrown.}
#   \item{...}{Not used.}
#
#   Note: If a constructor is not declared to be private nor protected, it
#   will be declared to be public.
# }
#
# \section{A constructor must be callable without arguments}{
#   The requirement that a constructor function should be callable without
#   arguments (e.g. \code{MyConstructor()}) is because that call is used
#   to create the static instance of a class.  The reason for this is that
#   a static instance of the class is created automatically when the
#   constructor is called \emph{the first time} (only), that is,
#   when the first of object of that class is created.
#   All classes have to have a static instance.
#
#   To make a constructor callable without arguments, one can either make
#   sure all arguments have default values or one can test for missing
#   arguments using \code{missing()}.
#   For instance the following defintion is \emph{not} correct:
#   \code{setConstructorS3("Foo", function(x) extend(Object(), "Foo", x=x))}
#   whereas this one is
#   \code{setConstructorS3("Foo", function(x=NA) extend(Object(), "Foo", x=x))}
# }
#
# \section{Code validation}{
#  If argument \code{enforceRCC} is @TRUE,
#  the class name is validated so it starts with a letter and it
#  also gives a @warning if its first letter is \emph{not} captial. The
#  reason for this is to enforce a naming convention that names classes
#  with upper-case initial letters and methods with lower-case initial
#  letters (this is also the case in for instance Java).
# }
#
# \examples{\dontrun{For a complete example see help(Object).}}
#
# \seealso{
#   To define a method see @see "R.methodsS3::setMethodS3".
#   For information about the R Coding Conventions, see
#   @see "RccViolationException".
#   For a thorough example of how to use this method see @see "Object".
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("setConstructorS3", "default", function(name, definition, private=FALSE, protected=FALSE, export=TRUE, static=FALSE, abstract=FALSE, trial=FALSE, deprecated=FALSE, envir=parent.frame(), enforceRCC=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that RCC naming conventions are followed.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (enforceRCC) {
    # Assert that the class name is a valid class name.
    firstLetter <- substring(name, 1,1);
    if (!is.element(tolower(firstLetter), letters))
      throw(RccViolationException("Class names must begin with a letter: ", name));

    # Check first letter
    if (firstLetter == tolower(firstLetter))
      throw(RccViolationException("Class names should be nouns starting with a capital letter: ", name));

    # Check if name contains . (period)
    if (regexpr("\\.", name) != -1)
      throw(RccViolationException("Class names must not contain . (period): ", name));
  }

  # Check for forbidden names.
  ns <- getNamespace("R.methodsS3");
  R.methodsS3_R.KEYWORDS <- get("R.KEYWORDS", envir=ns);
  if (is.element(name, R.methodsS3_R.KEYWORDS))
    throw(RccViolationException("Class names must not be same as a reserved keyword in R: ", name));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set up the modifiers
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (private)
    protection="private"
  else if (protected)
    protection="protected"
  else
    protection="public";

  # Create the modifiers
  modifiers <- protection;
  if (static == TRUE) modifiers <- c(modifiers, "static");
  if (abstract == TRUE) modifiers <- c(modifiers, "abstract");
  if (deprecated == TRUE) modifiers <- c(modifiers, "deprecated");
  if (trial == TRUE) modifiers <- c(modifiers, "trial");
  modifiers <- c(modifiers, "class");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the constructor function (by default in the parent frame)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create
  expr <- substitute({
      ns <- getNamespace("R.oo");
      if (exists("Class", mode="function", envir=ns, inherits=FALSE)) {
        R.oo_Class <- get("Class", mode="function", envir=ns, inherits=FALSE);
        fcn <- R.oo_Class(name, definition);
        rm(list="R.oo_Class");
      } else {
        # Only used for/by R.oo itself.
        fcn <- Class(name, definition);
      }
      rm(list="ns");
      attr(fcn, "export") <- export;
      attr(fcn, "modifiers") <- modifiers;
    }, list(fcn=as.name(name), name=name, definition=definition,
            export=export, modifiers=modifiers)
  );

  # Assign
  retValue <- eval(expr, envir=envir);

  invisible();
}) # setConstructorS3()




############################################################################
# HISTORY:
# 2012-08-29
# o Now setConstructorS3() no longer requires that R.oo is attached
#   ("loaded") - it's enough that it's namespace is loaded.
# 2012-04-17
# o Added argument 'export' to setConstructorS3().
# o CLEANUP: setConstructorS3() no longer sets attribute "formals".  It
#   has been deprecated since April 2003.
# 2007-06-09
# o Removed (incorrect) argument name 'list' from all substitute() calls.
# 2006-05-30
# o Renamed class from "ANY" to "default".
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2004-03-03
# o Moved deprecated setCl assS3() to its own file.
# 2003-07-17
# o Added Rdoc comments saying that the constructor function must be able
#   to be called without any arguments!
# 2003-04-13
# o Now setConstructorS3() sets the attribute "modifiers". "formals" will
#   still be around, but will be phased out.
# 2002-11-23
# o Renamed setClassS3() to setConstructorS3(), since this is what it is
#   actually doing. Keeping setClassS3() for backward compatibility but made
#   it deprecated.
# 2002-10-17
# o Removed obsolete "modifiers<-"().
# o Added also "Object" to the class attribute to make static methods to
#   work.
# 2002-10-16
# o There are times when
#     generic <- function(...) UseMethod()
#   is not working, for example
#     fcn <- get("generic"); fcn(myObj, ...);
#   For this reason, always do method dispatching using the name explicitly;
#     generic <- function(...) UseMethod("generic")
# 2002-10-15
# o Created from R.oo Object.R and ideas as described on
#    http://www.maths.lth.se/help/R/
############################################################################
