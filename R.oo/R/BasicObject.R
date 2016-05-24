###########################################################################/**
# @RdocClass BasicObject
#
# @title "A root class like Object but without references"
#
# \description{
#  R.oo\cr
#  \bold{Class BasicObject}\cr
#
#  public class \bold{BasicObject}\cr
# }
#
# @synopsis
#
# \arguments{
#   \item{core}{The core value of the object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
# \keyword{internal}
#*/###########################################################################
setConstructorS3("BasicObject", function(core=NULL) {
  # Create a new environment and wrap it up as a private field of a list.
  if (is.null(core))
    core <- NA;
  this <- core;
  class(this) <- unique(c("BasicObject", class(this)));

  if (getOption("R.oo::BasicObject/instantiationTime", FALSE)) {
    attr(this, "...instantiationTime") <- Sys.time();
  }

  this;
})


###########################################################################/**
# @RdocMethod isReferable
#
# @title "Checks if the object is referable or not"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @logical value, which by default is @TRUE for all
#  @see "BasicObject"'s.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("isReferable", "BasicObject", function(this, ...) {
  TRUE;
}) # isReferable()




###########################################################################/**
# @RdocMethod as.character
#
# @title "Gets a character string representing the object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \value{
#  Returns a @character string representation of the object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("as.character", "BasicObject", function(x, ...) {
  # To please R CMD check
  this <- x;

  paste(class(this)[1L], ": ", getInstantiationTime(this), sep="");
}) # as.character()





###########################################################################/**
# @RdocMethod getInstantiationTime
#
# @title "Gets the time when the object was instantiated"
#
# \description{
#  @get "title" (created) as a POSIXt object.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a POSIXt object, which extends class POSIXct.
# }
#
# \details{
#   The instantiation timestamp is set when the object is created, and
#   only of option \code{"R.oo::BasicObject/instantiationTime"} is @TRUE.
# }
#
# \seealso{
#   For more about time formats and POSIX see @see "base::DateTimeClasses".
#   @seeclass
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("getInstantiationTime", "BasicObject", function(this, ...) {
  time <- attr(this, "...instantiationTime");
  if (!is.null(time)) return(time);

  # Backward compatibility (due to a SPELLING ERROR in an earlier version)
  time <- attr(this, "...instanciationTime");

  NULL;
})




###########################################################################/**
# @RdocMethod hashCode
#
# @title "Gets a hash code for the object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seemethod "equals"
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("hashCode", "BasicObject", function(this, ...) {
  as.integer(getInstantiationTime(this));
})





###########################################################################/**
# @RdocMethod equals
#
# @title "Compares an object with another"
#
# \description{
#  @get "title" and returns @TRUE if they are equal.
#  The equal property must be
#
#  1) \emph{reflexive}, i.e. \code{equals(o1,o1)} should be @TRUE.
#
#  2) \emph{symmetric}, i.e. \code{equals(o1,o2)} is @TRUE if and only
#  if \code{equals(o2,o1)} is @TRUE.
#
#  3) \emph{transitive}, i.e. \code{equals(o1,o2)} is @TRUE and
#  \code{equals(o2,o3)} is @TRUE, then \code{equals(o1,o3)} should
#  be @TRUE.
#
#  5) \emph{consistent}, i.e. \code{equals(o1,o2)} should return the same
#  result on multiple invocations as long as noting has changed.
#
#  6) \code{equals(o1,NULL)} should return @FALSE.
#
#  By default, the method returns @TRUE if and only if the two
#  references compared refer to the same @see "BasicObject", i.e.
#  \code{( !is.null(obj) && (hashCode(this) == hashCode(obj)) )}.
# }
#
# @synopsis
#
# \arguments{
#   \item{other}{The other object this object should be compared to.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if the objects are equal, otherwise @FALSE.
# }
#
# \seealso{
#   @seeclass
# }
#
# @author
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("equals", "BasicObject", function(this, other, ...) {
  ( !is.null(other) && (hashCode(this) == hashCode(other)) );
})




###########################################################################/**
# @RdocMethod print
#
# @title "Prints an BasicObject"
#
# \description{
#  For all objects of class @see "BasicObject", this method will print the
#  value of \code{as.character()} of the object. Note that this function is
#  not called if the argument is not an object of class @see "BasicObject".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# @author
#
# \seealso{
#   @see "base::print.default"
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("print", "BasicObject", function(x, ...) {
  print(as.character(x));
}) # print()




###########################################################################/**
# @RdocMethod objectSize
#
# @title "Gets the size of the BasicObject in bytes"
#
# \description{
#   @get "title" by summing the sizes of all its members. For this reason,
#   the size of memory the BasicObject actually allocates might vary slighty.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer specifying the size of the object in number of bytes.
# }
#
# @author
#
# \seealso{
#   @see "utils::object.size".
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("objectSize", "BasicObject", function(this, ...) {
  object.size(this);
}) # objectSize()




###########################################################################/**
# @RdocMethod getFields
#
# @title "Returns the field names of an BasicObject"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{private}{If @TRUE, private fields will also be returned,
#   otherwise only public fields are returned.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector of field names.
# }
#
# @author
#
# \seealso{
#   To check if a field exists or not, see @seemethod "hasField".
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("getFields", "BasicObject", function(this, private=FALSE, ...) {
  members <- names(attributes(this));
  if (!private) {
    isPrivate <- (regexpr("^[.].*", members) != -1);
    members <- members[!isPrivate];
  }
  members;
}) # getFields()



###########################################################################/**
# @RdocMethod hasField
#
# @title "Checks if a field exists or not"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{field}{@vector of fields to be checked if they exists or not.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @logical @vector indicating for each field if it exists or not.
# }
#
# @author
#
# \seealso{
#   To get the fields of an Object, see @seemethod "getFields".
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("hasField", "BasicObject", function(this, field, ...) {
  !is.na(match(field, getFields(this, private=TRUE)));
}) # hasFields()




###########################################################################/**
# @RdocMethod attach
#
# @title "Attach an BasicObject to the R search path"
#
# \description{
#  Attach the members of an BasicObject to the \R search path.
#
#  If trying to attach the same BasicObject twice without detaching it
#  inbetween, a @warning will be generated and nothing will be done.
# }
#
# @synopsis
#
# \arguments{
#   \item{private}{If @TRUE, private fields will also be attached,
#     otherwise not.}
#   \item{pos}{The position at in search path where the BasicObject should be
#              inserted.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns @TRUE if the @see "BasicObject" was attached, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seemethod "detach" and @see "base::attach", @see "base::detach".
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("attach", "BasicObject", function(this, private=FALSE, pos=2, ...) {
  # To please R CMD check
  attachX <- base::attach;

  attachName <- as.character.BasicObject(this);
  if (is.element(attachName, search())) {
    warning(paste("Object is already attached:", attachName));
    return(invisible(FALSE));
  }

  if (is.list(this)) {
    attachX(unclass(this), name=attachName, pos=pos);
  } else {
    attachX(list(), name=attachName, pos=pos);
  }
  members <- names(attributes(this));

  for (member in members) {
    assign(member, attr(this, member), pos=pos);
  }

  return(invisible(TRUE));
}) # attach()





###########################################################################/**
# @RdocMethod detach
#
# @title "Detach an BasicObject from the R search path"
#
# \description{
#  Detach, from the \R search path, an BasicObject that has previously been
#  attached. If the BasicObject was not attached, a @warning will be
#  generated and nothing will be done.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns @TRUE if the BasicObject was detached, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seemethod "attach" and @see "base::attach", @see "base::detach".
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("detach", "BasicObject", function(this, ...) {
  attachName <- as.character.BasicObject(this);
  if (!is.element(attachName, search())) {
    warning(paste("Object is not attached:", attachName));
    return(invisible(FALSE));
  }

  pos <- which(search() == attachName);
  if (length(pos) == 1L) detach(pos=pos);

  return(invisible(TRUE));
}) # detach()



###########################################################################/**
# @RdocMethod extend
#
# @title "Extends another class"
#
# \description{
#   via a mechanism known as "parasitic inheritance".
#   Simply speaking this method "extends another class". What is actually
#   happening is that it creates an instance of class name \code{...className},
#   by taking another BasicObject instance and add \code{...className} to
#   the class list and also add all the named values in @... as fields to the
#   new instance.
#
#   The method should be used by the constructor of a class and nowhere else.
# }
#
# @synopsis
#
# \arguments{
#   \item{...className}{The name of new class.}
#   \item{...}{Named values representing the fields of the new instance.}
# }
#
# \value{
#  Returns an BasicObject of class \code{className}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("extend", "BasicObject", function(this, ...className, ...) {
  fields <- list(...);
  names <- names(fields);
  for (ii in seq_along(fields)) {
    name <- names[ii];
    if (is.null(name) || nchar(name) == 0) {
      callNames <- names(sys.call());
      callNames <- callNames[nchar(callNames) > 0];
      matchNames <- paste("^", callNames, sep="");
      for (jj in seq_along(matchNames)) {
        if (regexpr(matchNames[jj], "...className") != -1) {
          className <- sys.call()[[3]];
          throw("Could not set field of class (probably called ", className,
                ") because the field name is a prefix to the argument name ",
                "\"...className\": ", callNames[jj]);
        }
      } # for (jj ...)

      throw("Missing name of field #", ii, " in class definition: ", ...className);
    }
    attr(this, name) <- fields[[ii]];
  } # for (ii ...)

  class(this) <- c(...className, class(this));
  this;
}) # extend()




###########################################################################/**
# @RdocMethod newInstance
#
# @title "Creates a new instance of the same class as this object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of the corresponding
#     @see "BasicObject" class.}
# }
#
# \value{
#   Returns a reference to an instance of @see "BasicObject" or a subclass thereof.
# }
#
# @author
#
# \seealso{
#   @see "newInstance.Object".
#   @see "newInstance.Class".
#   @seeclass
# }
#
# @keyword programming
# @keyword methods
#*/###########################################################################
setMethodS3("newInstance", "BasicObject", function(this, ...) {
  # Creates a new instance of the same class
  clazz <- Class$forName(class(this)[1]);
  newInstance(clazz, ...);
}, private=TRUE)



###########################################################################/**
# @RdocMethod $
# @aliasmethod [[
#
# @title "Makes the fields and methods of an BasicObject accessable via the \$ and the [[ operator"
#
# \description{
#   @get "title".
# }
#
# \usage{
#   \method{$}{BasicObject}(this, name)
#   \method{[[}{BasicObject}(this, name)
# }
#
# \arguments{
#   \item{name}{The name of the field or method to be accessed.}
# }
#
# \value{
#  Returns the value of a field or a method (@function).
#  If no such field or method exists, @NULL is returned.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("$", "BasicObject", function(this, name) {
  memberAccessorOrder <- attr(this, ".memberAccessorOrder");
  if (is.null(memberAccessorOrder))
    memberAccessorOrder <- c(1,2,3,4);

  for (memberAccessor in memberAccessorOrder) {
    if (memberAccessor == 1) {
      firstChar <- substr(name, 1,1);
      isPrivate <- identical(firstChar, ".");
      isField <- (regexpr(" ", name) != -1);
      # Do not try to access private fields using a get<Name>() method,
      # because such a functionality means that the user *expects* that
      # there actually is a field called '.<name>', which he or she
      # should not do since it is a private field!
      if (!isField && !isPrivate && is.null(attr(this, "disableGetMethods"))) {
  	# 1. Is it a get<name>() method?
  	getName <- paste(c("get", toupper(firstChar),
  			 substr(name,2,nchar(name))),collapse="");
  	getMethodNames <- paste(getName, class(this), sep=".");
  	for (getMethodName in getMethodNames) {
            # TO DO/FIX ME: This part only works when packages are attached.
            # /HB 2013-10-08
            if (exists(getMethodName, mode="function")) {
  	    ref <- this;
  	    attr(ref, "disableGetMethods") <- TRUE;
  	    return( get(getMethodName, mode="function")(ref) );
  	  }
  	}
      }
    } else if (memberAccessor == 2) {

      # 2. Is it a field?
      value <- attr(this, name);
      if (!is.null(value))
  	return(value);

    } else if (memberAccessor == 3) {

      # 3. Is it a static S3 method?
      methodNames <- paste(name, class(this), sep=".");
      for (methodName in methodNames) {
          # TO DO/FIX ME: This part only works when packages are attached.
          # /HB 2013-10-08
          if (exists(methodName, mode="function")) {
#         # Alt 1. Rather "obfuscated" code
#         method <- get(methodName, mode="function");
#         fcn <- function(...) method(this, ...);
          # Alt 3. Using explicit UseMethod() code
          code <- sprintf("function(...) \"%s\"(this, ...)", name);
          fcn <- eval(base::parse(text=code));
          return(fcn);
  	}
      }
    }
  } # for (memberAccessor in memberAccessorOrder)

  # 5. Otherwise, return NULL.
  NULL;
}) # $()




###########################################################################/**
# @RdocMethod $<-
# @aliasmethod [[<-
#
# @title "Makes the fields and methods of an BasicObject assignable via the \$<- and the [[<- operator"
#
# \description{
#  @get "title".
# }
#
# \usage{
#   \method{$}{BasicObject}(this, name) <- value
#   \method{[[}{BasicObject}(this, name) <- value
# }
#
# \arguments{
#   \item{name}{The name of the \preformatted{set<Name>()} method or the
#     name of the field to be assigned the new value.}
#   \item{value}{The value to be assigned.}
# }
#
# \value{
#  Returns itself, i.e. \code{this}, as all \code{$<-} methods must do.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# \keyword{programming}
# \keyword{methods}
#*/###########################################################################
setMethodS3("$<-", "BasicObject", function(this, name, value) {
  memberAccessorOrder <- attr(this, ".memberAccessorOrder");
  if (is.null(memberAccessorOrder))
    memberAccessorOrder <- c(1,2,3,4);

  for (memberAccessor in memberAccessorOrder) {
    if (memberAccessor == 1) {
      # Do not try to access private fields using a set<Name>() method,
      # because such a functionality means that the user *expects* that
      # there actually is a field called '.<name>', which he or she
      # should not do since it is a private field!
      firstChar <- substr(name, 1,1);
      isPrivate <- identical(firstChar, ".");
      isField <- (regexpr(" ", name) != -1);
      if (!isField && !isPrivate && is.null(attr(this, "disableSetMethods"))) {
  	# 1. Is it a set<name>() method?

  	setName <- paste(c("set", toupper(firstChar),
  			 substr(name,2,nchar(name))),collapse="");
  	setMethodNames <- paste(setName, class(this), sep=".");
  	for (setMethodName in setMethodNames) {
            # TO DO/FIX ME: This part only works when packages are attached.
            # /HB 2013-10-08
            if (exists(setMethodName, mode="function")) {
  	    ref <- this;
  	    attr(ref, "disableSetMethods") <- TRUE;
  	    this <- get(setMethodName, mode="function")(ref, value);
  	    attr(this, "disableSetMethods") <- NULL;
  	    return(this);
  	  }
  	}
      }
    } else if (memberAccessor == 2) {

      # 2. If there exists a field, assign the value to that field.
      if (!is.null(attr(this, name))) {
  	attr(this, name) <- value;
        return(this);
      }
    } else if (memberAccessor == 4) {
      # 4. Otherwise, assign the value to a new field.
      attr(this, name) <- value;
      return(this);
    }
  } # for (memberAccessor in memberAccessorOrder)

  this;
}) # $<-()


setMethodS3("[[", "BasicObject", function(this, name) {
  memberAccessorOrder <- attr(this, ".memberAccessorOrder");
  if (is.null(memberAccessorOrder))
    memberAccessorOrder <- c(1,2,3,4);

  for (memberAccessor in memberAccessorOrder) {
    if (memberAccessor == 1) {
      firstChar <- substr(name, 1,1);
      isPrivate <- identical(firstChar, ".");
      isField <- (regexpr(" ", name) != -1);
      # Do not try to access private fields using a get<Name>() method,
      # because such a functionality means that the user *expects* that
      # there actually is a field called '.<name>', which he or she
      # should not do since it is a private field!
      if (!isField && !isPrivate && is.null(attr(this, "disableGetMethods"))) {
  	# 1. Is it a get<name>() method?
  	getName <- paste(c("get", toupper(firstChar),
  			 substr(name,2,nchar(name))),collapse="");
  	getMethodNames <- paste(getName, class(this), sep=".");
  	for (getMethodName in getMethodNames) {
            # TO DO/FIX ME: This part only works when packages are attached.
            # /HB 2013-10-08
            if (exists(getMethodName, mode="function")) {
  	    ref <- this;
  	    attr(ref, "disableGetMethods") <- TRUE;
  	    return( get(getMethodName, mode="function")(ref) );
  	  }
  	}
      }
    } else if (memberAccessor == 2) {

      # 2. Is it a field?
      value <- attr(this, name);
      if (!is.null(value))
  	return(value);

    } else if (memberAccessor == 3) {

      # 3. Is it a method?
      methodNames <- paste(name, class(this), sep=".");
      for (methodName in methodNames) {
          # TO DO/FIX ME: This part only works when packages are attached.
          # /HB 2013-10-08
          if (exists(methodName, mode="function")) {
  	  method <- get(methodName, mode="function");
  	  return( function(...) method(this, ...) );
  	}
      }
    }
  } # for (memberAccessor in memberAccessorOrder)

  # 5. Otherwise, return NULL.
  NULL;
}) # "[["()


setMethodS3("[[<-", "BasicObject", function(this, name, value) {
  UseMethod("$<-");
}) # "[[<-"()


############################################################################
# HISTORY:
# 2012-12-28
# o Replaced all data.class(obj) with class(obj)[1].
# 2012-12-18
# o R CMD check for R devel no longer gives a NOTE about attach().
# 2012-10-14
# o Now <BasicObject>$<staticFcn>(...) calls <staticFcn>(<BasicObject>, ...).
# 2012-06-22
# o ROBUSTNESS: Now constructor BasicObject() is guaranteed to return
#   an object with non-duplicated class attribute elements.
# o GENERALIZATION: Added newInstance() for BasicObject.
# 2011-04-02
# o Added option "R.oo::Object/instantiationTime", which controls
#   whether the instantiation timestamp should be set when instantiating
#   an Object. Analogoulsy, option "R.oo::BasicObject/instantiationTime"
#   controls ditto for a BasicObject.
# 2008-05-28
# o SPELL CORRECTION: Used '...instanciation' instead of 'instantiation'.
# 2005-02-15
# o Added arguments '...' in order to match any generic functions.
# 2004-10-17
# o Added Rdoc comments.
# 2003-01-18
# o Replaced all occurences of getClass() with data.class(). Will change
#   the use of getClass() in the future to return a Class object.
# 2002-11-05
# o Added class(core) to the class list.
# 2002-11-04
# o Created to be the upcoming root class which does not create a reference
#   object by default.
############################################################################
