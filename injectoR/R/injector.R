# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

#' Dependency injection framework
#'
#' @author levk
#' @docType package
#' @name injectoR
NULL;

#' Root binder
.binder <- base::new.env (parent = base::emptyenv ());

#' Binder factory
#' 
#' @param parent of the new binder, injection will propagate up the 
#' parent stack looking for keys; if omitted defaults to root binder
#' @param callback called with the newly created binder and the
#' result is returned; if omitted just the new binder is returned
#' @return result of the injected callback if one is specified,
#' otherwise the new binder
#' @export
#' @examples
#' b <- binder ()
binder <- function (parent = .binder, callback = function (binder) binder)
  callback (base::new.env (parent = parent));

#' Singleton scope, bindings of this scope are provided once, on
#' initial demand
#' 
#' @param provider unscoped delegate, no argument function responsible
#' for provision
#' @export
#' @examples
#' define (three = function () 3, scope = singleton, binder = binder ())
singleton <- function (provider)
  (function (value) function () if (base::is.null (value)) value <<- provider () else value) (NULL);

#' Default scope, bindings are provisioned each time a bean is
#' injected
#' 
#' @param provider unscoped delegate, no argument function responsible
#' for provision
default <- function (provider) provider;

#' Creates a key to factory binding
#' 
#' @param ... injectable bean identifier to factory mappings, the key
#' is the name is matched to a parameter name during injection, the
#' factory responsible for provisioning of the bean, a factory may
#' accept any number of arguments in which case the framework will
#' attempt to inject the argument if a binding to the parameter name
#' exists; if it does not, that argument will not be injected, in
#' which case it is the factory's responsibility to deal with a
#' missing argument
#' @param scope of the bean, wraps the injected factory call
#' specifying provisioning strategy, if omitted a new bean instance
#' will be provisioned each time injection is requested; injectoR also
#' ships with with the singleton scope which will provide once and
#' cache the bean for subsequent calls. Interface allows for custom
#' scoping, the scope parameter must be a function accepting key (name)
#' and the provider - the wrapped injected factory call - a function
#' accepting no parameters responsible for actual provisioning
#' @param binder for this binding, if omitted the new binding is added
#' to the root binder
#' @export
#' @examples
#' define (hello = function () 'world', binder = binder ())
define <- function (..., scope = default, binder = .binder) {
  bindings <- base::list (...);
  base::lapply (base::names (bindings), function (key)
    binder[[ key ]] <- scope (function () inject (bindings[[ key ]], binder)));
  binder;
};

#' Aggregates multiple factories under one key
#' 
#' @param key injectable bean identifier
#' @param scope of the bean, wraps the injected factory call
#' specifying provisioning strategy, if omitted a new bean instance
#' will be provisioned each time injection is requested; injectoR also
#' ships with with the singleton scope which will provide once and
#' cache the bean for subsequent calls. Interface allows for custom
#' scoping, the scope parameter must be a function accepting key (name)
#' and the provider - the wrapped injected factory call - a function
#' accepting no parameters responsible for actual provisioning
#' @param combine aggregation procedure for combination of context
#' and inherited values, a function accepting a list of injectable
#' values from the current binder context and a no argument function
#' to retrieve values of the parent context; if omitted will the binding
#' will aggregate all values
#' @param binder for this binding, if omitted the binding is added to
#' the root binder
#' @return a function accepting one or more factories for adding
#' elements to the binding; naming the factories will result in named
#' values injected; optionally accepts a scope for the bindings, if
#' omitted defaults to provide on injection; please be aware that the
#' scope is called without key for unnamed multibinding
#' @export
#' @examples
#' multibind ('keys', binder = binder ()) (function () 'skeleton')
multibind <- function (key, scope = default,
                       combine = function (this, parent) base::c (this, parent), binder = .binder) 
  if (base::exists (key, envir = binder, inherits = FALSE)) base::attr (binder[[ key ]], 'multibind') else {
    providers <- base::list ();
    binder[[ key ]] <- scope (function () {
      parent <- base::parent.env (binder);
      combine (base::lapply (providers, function (provider) provider ()),
               if (base::exists (key, envir = parent)) base::get (key, envir = parent) () else base::list ())
    });
    base::attr (binder[[ key ]], 'multibind') <- function (..., scope = default)
      providers <<- base::c (providers,
                             base::lapply (base::list (...), function (factory) {
        base::force (factory);
        scope (function () inject (factory, binder));
      }));
  };

#' Shims libraries
#' 
#' @param ... zero or more library names to shim binding each exported
#' variable to the binder; if a library name is specified in a named
#' list format (for example shim(s4='stats4',callback=function(s4.AIC)))
#' all exported variable names from that library will be prepended with
#' that name and a dot (as in the example); if a library cannot be
#' loaded, no bindings are created for that library and no errors are
#' thrown (but there is an error to console as reported by
#' requireNamespace)
#' @param library.paths to use for loading namespace
#' @param callback injected for convenience using the binder specified
#' after shim is completed, if omitted the call returns the binder
#' @param binder for this shim
#' @return result of the callback if specified, binder otherwise
#' @export
#' @examples
#' shim ('injectoR', callback = function (inject) inject, binder = binder ())
shim <- function (..., library.paths = .libPaths (), callback = function () binder, binder = .binder) {
  exports <- base::unlist (base::lapply (base::list (...), function (package)
    if (!base::is.character (package)) base::stop ("Library name list must consist of strings only")
    else if (base::requireNamespace (package, lib.loc = library.paths)) (function (namespace)
      base::lapply (stats::setNames (nm = base::getNamespaceExports (namespace)),
                    function (export)
                      base::getExportedValue (namespace, export))) (base::loadNamespace (package, lib.loc = library.paths))));
  base::lapply (1:base::length (exports),
                function (i) binder[[ base::names (exports)[ i ] ]] <- singleton (function () exports[[ i ]]));
  inject (callback, binder);
};

#' Injects the callback function
#' 
#' @param callback function to inject, a function accepting arguments
#' to be matched to injectable keys; no errors are thrown if no binding
#' is found for a key, this is the intended mechanic for optional
#' injection, if the callback is able to deal with a missing argument
#' the argument becomes optional
#' @param binder containing the injectables, defaults to root binder if
#' omitted
#' @return result of the injected callback evaluation
#' @export
#' @examples
#' inject (function (two) two, define (two = function () 2, binder = binder ()))
#' inject (function (power) power (2, 4), 
#'         define (power = function (power) function (x, n) if (n < 1) 1 else x * power (x, n - 1)))
#' inject (function (fibonacci) fibonacci (8),
#'         define (fibonacci = function (fibonacci)
#'           function (n) if (n < 3) 1
#'                        else fibonacci (n - 1) + fibonacci (n - 2), binder = binder ()))
inject <- function (callback, binder = .binder) {
  args <- base::new.env (parent = base::environment (callback));
  args$missing <- function (x) {
    key <- base::as.character (base::match.call ()$x);
    if (!base::identical (base::parent.frame (), args)) base::missing (x)
    else if (!base::match (key, base::names (base::formals (callback)), nomatch = 0) > 0)
      base::stop ("'missing' can only be used for arguments")
    else !base::exists (key, envir = binder);
  };
  base::lapply (base::names (base::formals (callback)), function (key)
    base::makeActiveBinding (key, (function (value) function (x) 
      if (!base::missing (x)) value <<- x
      else if (base::is.null (value))
        value <<- if (base::exists (key, envir = binder)) base::get (key, envir = binder) ()
        else if (base::formals (callback)[[ key ]] != '') base::formals (callback)[[ key ]]
        else base::stop (base::paste ("Unbound dependency on", key))
      else value) (NULL), args));
  base::eval (base::body (callback), args);
};
