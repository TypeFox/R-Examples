#'
#' @title Groovy Scripting Language Integration
#'
#' @description Functions that integrate the Groovy scripting language with the R Project for Statistical Computing.
#'
#' @details From \href{http://en.wikipedia.org/wiki/Groovy_(programming_language)}{Wikipedia}:
#'
#' "Groovy is an object-oriented programming language for the Java platform. It is a dynamic language with features
#' similar to those of Python, Ruby, Perl, and Smalltalk. It can be used as a scripting language for the Java Platform,
#' is dynamically compiled to Java Virtual Machine (JVM) bytecode, and interoperates with other Java code and
#' libraries."
#'
#' @examples
#'  \dontrun{
#' #
#' # Installation Example
#' #
#' # Since this package does not ship with Groovy the user needs to specify the Groovy jars prior
#' # to using the package -- here's an example how this is accomplished:
#' #
#' groovyJars <- list (
#'     "C:/Temp/groovy.jars/groovy-2.4.5-indy.jar",
#'     # OTHER JAR FILES...
#' )
#'
#' options(GROOVY_JARS=groovyJars)
#'
#' library(rGroovy)
#'
#' Execute (groovyScript="print 'Hello world!'")
#' }
#'
#' @seealso \href{http://groovy.codehaus.org/}{Groovy}
#' @seealso \href{http://www.groovy-lang.org/indy.html}{Invoke Dynamic}
#'
#' @import rJava
#'
#' @docType package
#'
#' @name rGroovy
#'
NULL

#'
#' An environment which is used by this package when managing package-scope
#' variables.
#'
.rGroovyAPI.env <- new.env()

.onLoad <- function (libname, pkgname) {
    #
    # Note that we cannot add documentation to this function -- if we do, when roxygen is executed we will see an error
    # about a missing name and the only way to get past this error is to add @name to this function or leave the
    # documentation off altogether.
    #
    # Function starts the Java virtual machine.
    #
    # @param libname The library name.
    #
    # @param pkgname The package name.
    #

    groovyJars <- getOption("GROOVY_JARS")

    .jpackage(pkgname, lib.loc = libname, morePaths = groovyJars)
}

#' Funtion does nothing.
#'
#' @param libpath The library path.
#'
.onUnload <- function (libpath) {
}

#' Function sets the global instance of GroovyShell that will be used by the Evaluate function whenever it is called
#' with a NULL GroovyShell parameter.
#'
#' @param binding An instance of \href{http://docs.groovy-lang.org/latest/html/api/index.html?groovy/lang/Binding.html}{groovy.lang.Binding}.
#'
#' @examples
#'  \dontrun{
#'  Initialize ()
#'  }
#'
#' @export
#'
Initialize <- function (binding = NULL) {

    if (!is.null (.rGroovyAPI.env$groovyShell)) {
        warning (
            paste (
                "Initialize has been invoked more than once -- are you sure ",
                "this is what you intended to do?"
            )
        )
    }

    if (is.null (binding)) {
        binding <- .jnew("groovy.lang.Binding")
    }

    groovyShell <- .jnew("groovy.lang.GroovyShell", binding)

    assign("groovyShell", groovyShell, envir = .rGroovyAPI.env)
}

#' Function evaluates (executes) the groovy script and returns the result.
#'
#' @param groovyShell The groovyShell with which to execute the specified groovy script. Note that the groovyShell can
#'  be NULL, however if this is NULL then the Initialize function must have been called so that a global groovyShell
#'  instance will be available in the environment otherwise an exception is raised.
#'
#' @param groovyScript The groovy script being executed.
#'
#' @return The result of the script execution.
#'
#' @examples
#'  \dontrun{
#'  Initialize ()
#'  Evaluate (groovyScript="print 'Hello world!'")
#'  }
#'
#' @export
#'
Evaluate <- function (
    groovyShell = NULL,
    groovyScript
) {
    if (is.null (groovyScript)) {
        stop ("The groovyScript parameter cannot be NULL.")
    }

    if (is.null (groovyShell)) {
        groovyShell <- .rGroovyAPI.env$groovyShell
    }

    if (is.null (groovyShell)) {
        stop ("Both the local and global groovyShell(s) are null -- did you forget to call the Initialize function?")
    }

    tryCatch (
        result <- groovyShell$evaluate (groovyScript),
        Throwable = function (exception) {
            stop (
                paste (
                    "An exception was thrown when executing the groovy ",
                    "script -- details follow.",
                    groovyScript, exception$getMessage(), sep=""))
    })

    return (result)
}

#' Function executes the groovy script and returns the result. Execute differs from Evaluate in that references to
#' Groovy objects are not required. The call to Initialize is not required in order to call this function either however
#' keep in mind that a new instance of \href{http://docs.groovy-lang.org/latest/html/api/groovy/lang/GroovyShell.html}{groovy.lang.GroovyShell}
#' will be used every time this function is called.
#'
#' @param groovyScript The groovy script being executed.
#'
#' @param variables The variables that will be passed to the binding that is used when the groovyScript is executed.
#'
#' @examples
#'  \dontrun{
#'  variables <- list ()
#'
#'  variables["name"] = "Tom"
#'  variables["day"]  = "Wednesday"
#'
#'  groovyScript <- "return \"Hello ${name}, how are you doing? Today is ${day}.\""
#'
#'  result <- Execute (groovyScript=groovyScript, variables=variables)
#'  result
#'  }
#'
#' @export
#'
Execute <- function (
    groovyScript,
    variables = list ()
) {
    binding <- .jnew("groovy.lang.Binding")

    groovyShell <- .jnew("groovy.lang.GroovyShell", binding)

    if (is.null(variables)) {
        variables <- list ()
    }

    for (nextName in ls(variables)) {
        nextValue <- variables[[nextName]]
        binding$setVariable (nextName, nextValue)
    }

    return (Evaluate (groovyShell=groovyShell, groovyScript=groovyScript))
}

#' Function prints some information about this package.
#'
#' @export
#'
About <- function () {
    cat (
        " ***********************************************************\n",
        "***                                                     ***\n",
        "***            Welcome to the rGroovy Package           ***\n",
        "***                                                     ***\n",
        "***                    version 1.0.                     ***\n",
        "***                                                     ***\n",
        "***                Follow us on LinkedIn:               ***\n",
        "***                                                     ***\n",
        "***       https://www.linkedin.com/company/229316       ***\n",
        "***                                                     ***\n",
        "***                Follow us on Twitter:                ***\n",
        "***                                                     ***\n",
        "***        https://twitter.com/CoherentLogicCo          ***\n",
        "***                                                     ***\n",
        "***********************************************************\n")
}
