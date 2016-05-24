###########################################################################/**
# @RdocDocumentation "Non-documented objects"
#
# % The HttpDaemon class
# @alias appendRootPaths
# @alias finalize.HttpDaemon
# @alias getConfig
# @alias getCount
# @alias getCount.HttpDaemon
# @alias setCount
# @alias setCount.HttpDaemon
# @alias getDefaultFilenamePattern
# @alias getHttpRequest
# @alias getPort
# @alias getRootPath
# @alias getRootPaths
# @alias insertRootPaths
# @alias isStarted
# @alias openUrl
# @alias processRsp
# @alias restart
# @alias restart.default
# @alias setRootPaths
# @alias sourceTcl
# @alias startHelp
# @alias stop
# @aliaswriteResponse
#
# % The HttpDaemonResponse class
# @alias flush
#
# % The HttpRequest class
# @alias getContentLength
# @alias getContentType
# @alias getContextPath
# @alias getContextPath.HttpRequest
# @alias getDateHeader
# @alias getDateHeader.HttpRequest
# @alias getHeader
# @alias getHeader.HttpRequest
# @alias getParameter
# @alias getParameters
# @alias getProtocol
# @alias getQueryString
# @alias getQueryString.HttpRequest
# @alias getRealPath
# @alias getRemoteAddress
# @alias getRemoteHost
# @alias getRemoteUser
# @alias getRemoteUser.HttpRequest
# @alias getRequestUri
# @alias getRequestUri.HttpRequest
# @alias getRequestUrl
# @alias getRequestUrl.HttpRequest
# @alias getScheme
# @alias getServerName
# @alias getServerPort
# @alias getServletPath
# @alias getServletPath.HttpRequest
# @alias hasParameter
# @alias nbrOfParameters
#
# % The RspResponse class
# @alias getOutput
# @alias import
# @alias write
#
# % The RspLanguage class
# @alias escape
# @alias getComment
# @alias getLanguage
# @alias getNewline
# @alias getVerbatim
#
# % The HtmlRspLanguage class
#
# % The RspPage class
#
# % RRspPackage class
# @alias capabilitiesOf
# @alias isCapableOf
#
# % Miscellanoues stuff
# @alias hexToInt
# @alias includeRsp.default
# @alias stop.default
# @alias write.default
# @alias urlDecode
# @alias compileRsp0
# @alias compileRsp0.default
# @alias epsDev
# @alias dvips
# @alias ps2pdf
# @alias rspCapture
# @alias sourceWithTrim
# @alias sourceWithTrim.default
# @alias sourceRspV2
# @alias sourceRspV2.default
#
# @alias exprToCode
# @alias getAttribute
# @alias getAttributes
# @alias hasAttribute
# @alias setAttribute
# @alias setAttributes
# @alias getCode
# @alias getEcho
# @alias getFile
# @alias getReturn
# @alias getContent
# @alias makeSourceCode
# @alias parseRaw
# @alias toR
# @alias toSourceCode
# @alias flatten
# @alias preprocess
# @alias getType
# @alias findProcessor
# @alias hasProcessor
# @alias process
# @alias tangle
# @alias getCompleteCode
# @alias getSuffixSpecs
# @alias mergeTexts
# @alias getWrap
# @alias asRspString
# @alias dropEmptyText
# @alias trimNonText
# @alias getMetadata
# @alias setMetadata
# @alias escapeRspContent
# @alias escapeRspTags
# @alias extensionToIMT
# @alias parseInternetMediaType
# @alias unescapeRspTags
# @alias nbrOfLines
# @alias requireAttributes
# @alias getItem
# @alias getItem.RspPreprocessingException
# @alias getMessage.RspPreprocessingException
# @alias getNameContentDefaultAttributes
# @alias getNameContentDefaultAttributes.RspDirective
# @alias getInclude
# @alias tidy
# @alias view
#
# @alias parse
# @alias parse.default
#
# @alias withoutGString
# @alias wstring
# @alias wstring.default
#
# @alias getFileSize
# @alias getFileSize.RspFileProduct
#
#
# \description{
# This page contains aliases for all "non-documented" objects that
# \code{R CMD check} detects in this package.
#
# Almost all of them are \emph{generic} functions that have specific
# document for the corresponding method coupled to a specific class.
# Other functions are re-defined by \code{setMethodS3()} to
# \emph{default} methods. Neither of these two classes are non-documented
# in reality.
# The rest are deprecated methods.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################

############################################################################
# HISTORY:
# 2005-02-18
# o Created to please R CMD check.
############################################################################
