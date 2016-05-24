############################################################################
# This source code file contains constructor and function definitions that
# are used for loading this package only.
############################################################################
# To please R CMD check
attachX <- base::attach;

attachX(list(
  Object = function(core=NA, finalize=TRUE) {
    # Create a new environment and wrap it up as a private field of a list.
    this <- core;
    this.env <- new.env();
    attr(this, ".env") <- this.env;
    class(this) <- "Object";

    if (getOption("R.oo::Object/instantiationTime", FALSE)) {
      attr(this, "...instantiationTime", Sys.time());
    }

    finalizer <- function(env) {
      # Nothing needed to be done in this temporary finalizer factory,
      # because it is only utilized by this very package and none of
      # the classes in this package created Object:s that needs to be
      # finalized.
    } # finalizer()

    # Should this Object be finalized?
    if (finalize) {
      onexit <- getOption("R.oo::Object/finalizeOnExit", FALSE);
      reg.finalizer(this.env, finalizer, onexit=onexit);
    }
    assign("...finalize", finalize, envir=this.env, inherits=FALSE);

    this;
  },

  extend = function(this, ...className, ..., ...finalize=TRUE) {
    fields <- list(...);
    names <- names(fields);
    this.env <- attr(this, ".env");
    for (name in names)
      assign(name, fields[[name]], envir=this.env);
    class(this) <- c(...className, class(this));

    # Override this (=unregister finalizer) according to argument
    # '...finalize' of extend()?
    if (!is.na(...finalize) && !isTRUE(...finalize)) {
      # Unregister finalizer (by registering a dummy one)
      reg.finalizer(this.env, f=function(...) {});
    }

    this;
  },

  Class = function(name=NULL, constructor=NULL) {
    if (is.null(name)) {
      constructor <- NA;
    } else if (!is.function(constructor)) {
      throw("Argument 'constructor' must be a function: ", mode(constructor));
    }

    # This is an example where one wants to have the core of an Object to not
    # be NA, but something else.
    this <- extend(Object(constructor), "Class",
      .staticInstance = NULL
    );

    this;
  }
), name="R.oo");

# Cleanup
rm(list="attachX");


############################################################################
# HISTORY:
# 2014-02-22
# o Now the temporary registered finalizer for Object:s returns a
#   finalizer function that does nothing.  This is possible because
#   none of the classes in this package produces Object:s that needs
#   finalization.
# 2014-01-05
# o BUG FIX: The temporary finalizer() registered for Object while
#   loading the R.oo package itself would cause cyclic loading of R.oo.
#   The reason was that it checked whether R.oo was available or not,
#   by only looking at attached namespaces but not loaded ones.
# 2013-10-13
# o Added argument 'finalize' to Object() and '...finalize' to extend()
#   for Object.  The latter override the former.
# 2012-12-18
# o R CMD check for R devel no longer gives a NOTE about attach().
# 2012-11-28
# o LIMITATION: Registered finalizer for pure Object:s (i.e. excluding
#   those which are of a subclass of Object) will no longer be called
#   if the R.oo package has been detached.
# 2011-04-02
# o Added option "R.oo::Object/finalizeOnExit".
# o Added option "R.oo::Object/instantiationTime", which controls
#   whether the instantiation timestamp should be set when instantiating
#   an Object. Analogoulsy, option "R.oo::BasicObject/instantiationTime"
#   controls ditto for a BasicObject.
# 2011-03-11
# o Added explicit 'onexit=FALSE' to all reg.finalizer():s so it is clear
#   that they are not finalized when quitting R.  See 050.Object.R too.
# 2008-05-28
# o SPELL CORRECTION: Used '...instanciation' instead of 'instantiation'.
# 2008-01-10
# o Made the registered finalizer calling finalize() more error prone.
# 2007-08-29
# o BUG FIX: If Object:s are garbage collected after R.oo has been detached,
#   the error 'Error in function (env) : could not find function "finalize"'
#   would be thrown, because the registered finalizer hook tries to call
#   the generic function finalize() in R.oo.  We solve this by trying to
#   reload R.oo (and the unload it again).  Special care was taken so that
#   Object:s allocated by R.oo itself won't cause an endless loop.
# 2005-06-10
# o Added reg.finalizer() to Object() for pure Object:s. However, it must
#   be done in extend.Object() too.
# 2004-10-18
# o Updated the arguments for extend() so that they are identical to the
#   ones in extend.Object.
# 2002-12-15
# o Added reg.finalizer() to the Object class.
# 2002-11-04
# o Added the feature to timestamp each object when it is instanciated.
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
