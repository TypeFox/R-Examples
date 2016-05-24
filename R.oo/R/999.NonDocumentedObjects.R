###########################################################################/**
# @RdocDocumentation "Non-documented objects"
#
# % The BasicObject class
# @alias getInstanciationTime.default
# @alias isReferable
#
# % The Class class
# @alias forName
# @alias getDetails
# @alias getKnownSubclasses
# @alias getMethods
# @alias getMethods.default
# @alias getName
# @alias getPackage
# @alias getStaticInstance
# @alias getSuperclasses
# @alias isAbstract
# @alias isBeingCreated
# @alias isDeprecated
# @alias isPrivate
# @alias isProtected
# @alias isPublic
# @alias isStatic
# @alias newInstance
#
# % The Exception class
# @alias getCalls
# @alias getCall
# @alias getLastException
# @alias getMessage
# @alias getStackTrace
# @alias getStackTraceString
# @alias getWhen
# @alias printStackTrace
#
# % The Object class
# @alias attach
# @alias attach.default
# @alias attachLocally
# @alias clone
# @alias clearLookupCache
# @alias clearCache
# @alias detach
# @alias detach.default
# @alias finalize
# @alias gc
# @alias getFields
# @alias getInstanciationTime
# @alias getInstanciationTime.default
# @alias getInstantiationTime
# @alias getInternalAddress
# @alias getFieldModifier
# @alias getFieldModifiers
# @alias hasField
# @alias load
# @alias load.default
# @alias novirtual
# @alias registerFinalizer
# @alias save
# @alias save.default
# @alias staticCode
#
# % The Package class
# @alias getAuthor
# @alias getBundle
# @alias getBundlePackages
# @alias getChangeLog
# @alias getClasses
# @alias getClasses.default
# @alias getContribUrl
# @alias getContents
# @alias getDataPath
# @alias getDate
# @alias getDescription
# @alias getDescriptionFile
# @alias getDevelUrl
# @alias getDocPath
# @alias getEnvironment
# @alias getExamplePath
# @alias getHistory
# @alias getHowToCite
# @alias getLicense
# @alias getMaintainer
# @alias getManPath
# @alias getNews
# @alias getPath
# @alias getPosition
# @alias getTitle
# @alias getUrl
# @alias getVersion
# @alias isLoaded
# @alias isOlderThan
# @alias showChangeLog
# @alias showContents
# @alias showDescriptionFile
# @alias showHistory
# @alias showHowToCite
# @alias showNews
# @alias startupMessage
# @alias unload
#
# % The RccViolationException class
# @alias getRccUrl
#
# % The Rdoc class
# @alias argsToString
# @alias check
# @alias compile
# @alias createManPath
# @alias createName
# @alias declaration
# @alias escapeRdFilename
# @alias getClassS4Usage
# @alias getKeywords
# @alias getNameFormat
# @alias getObject
# @alias getObject.Rdoc
# @alias getPackageNameOf
# @alias getRdDeclaration
# @alias getRdHierarchy
# @alias getRdMethods
# @alias getRdTitle
# @alias getUsage
# @alias hierarchy
# @alias isKeyword
# @alias isVisible
# @alias methodsInheritedFrom
# @alias setManPath
# @alias setNameFormat
#
# % The RdocException class
# @alias getSource
#
# % Trial functions
# @alias gc.default
# @alias callSuperMethodS3
# @alias callSuperMethodS3.default
#
# % Deprecated functions
# @alias setClassS3
# @alias setClassS3.default
# @alias getClass.BasicObject
# @alias getClass.default
#
# \description{
#   This page contains aliases for all "non-documented" objects that 
#   \code{R CMD check} detects in this package. 
#
#   Almost all of them are \emph{generic} functions that have specific 
#   document for the corresponding method coupled to a specific class. 
#   Other functions are re-defined by \code{setMethodS3()} to 
#   \emph{default} methods. Neither of these two classes are non-documented
#   in reality.
#   The rest are deprecated methods.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################

############################################################################
# HISTORY:
# 2012-02-29
# o CLEANUP: Dropped aliases for non-existing environment[.default]().
# 2005-02-10
# o Created to please R CMD check.
############################################################################
