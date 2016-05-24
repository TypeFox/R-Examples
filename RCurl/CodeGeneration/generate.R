library(RGCCTranslationUnit)
#library(RAutoGenRunTime)   # for bitlist.

tu = parseTU("enum.c.t00.tu")
edefs = lapply(getEnumerations(tu), resolveType, tu)

curl.h = "/usr/local/include/curl/curl.h"
curl.h = "~/Downloads/curl-7.19.4/include/curl/curl.h"
curl.h = "~/Downloads/curl-7.14.0/include/curl/curl.h"

cpp = getCppDefines(curl.h)
defs = RGCCTranslationUnit:::processDefines(cpp, tu = tu, filter = NULL)
i = defs$macros %in% names(edefs$CURLoption@values)

##############

writeEnumGenerationRCode(edefs[["CURLoption"]]@values, "../inst/enums/Renums.c", includes = '<curl/curl.h>')

##########

ifdef = 
c("KEYPASSWD", "DIRLISTONLY", "APPEND", "KRBLEVEL", "USE_SSL", 
"TIMEOUT_MS", "CONNECTTIMEOUT_MS", "HTTP_TRANSFER_DECODING", 
"HTTP_CONTENT_DECODING", "NEW_FILE_PERMS", "NEW_DIRECTORY_PERMS", 
"POSTREDIR", "OPENSOCKETFUNCTION", 
"OPENSOCKETDATA", "COPYPOSTFIELDS", "PROXY_TRANSFER_MODE", "SEEKFUNCTION", 
"SEEKDATA", "CRLFILE", "ISSUERCERT", "ADDRESS_SCOPE", "CERTINFO", 
"USERNAME", "PASSWORD", "PROXYUSERNAME", "PROXYPASSWORD",
"SSH_HOST_PUBLIC_KEY_MD5", "NOPROXY", "TFTP_BLKSIZE", 
"SOCKS5_GSSAPI_SERVICE", "SOCKS5_GSSAPI_NEC", 
"PROTOCOLS", "REDIR_PROTOCOLS",
"SSH_AUTH_TYPES",
"SSH_PUBLIC_KEYFILE",
"SSH_PRIVATE_KEYFILE",
"FTP_SSL_CCC",
#
"COOKIELIST", "IGNORE_CONTENT_LENGTH", "FTP_SKIP_PASV_IP", 
"FTP_FILEMETHOD", "LOCALPORT", "LOCALPORTRANGE", "CONNECT_ONLY", 
"CONV_FROM_NETWORK_FUNCTION", "CONV_TO_NETWORK_FUNCTION", "CONV_FROM_UTF8_FUNCTION", 
"MAX_SEND_SPEED_LARGE", "MAX_RECV_SPEED_LARGE", "FTP_ALTERNATIVE_TO_USER", 
"SOCKOPTFUNCTION", "SOCKOPTDATA", "SSL_SESSIONID_CACHE"
)

 # For testing in configure.in
cat(sprintf("RCURL_CHECK_ENUM(CURLOPT_%s)\n", ifdef))
cat(sprintf("RCURL_CHECK_ENUM(%s)", names(vals)), sep = "\n")


 # Get the options table.

txt = readLines(curl.h)
cinit = grep("^[[:space:]]*CINIT\\(", txt, val = TRUE)
names(cinit) = gsub(".*CINIT\\(([A-Z0-9_]+),.*", "\\1", cinit)

cinit[ifdef] = sprintf("#ifdef HAVE_CURLOPT_%s\n%s\n#endif\n", ifdef, cinit[ifdef])


cat(cinit,   # grep("CINIT\\(", cinit, val = TRUE)
      sprintf("#ifdef %s\n{ \"%s\", %s}%s\n#endif\n", 
                 names(defs$macros[i]),
                 gsub("^CURLOPT_", "", names(defs$macros[i])),
                 defs$macros[i],
                 c(rep(",", length(defs$macros[i]) - 1), "")),
         sep = "\n", file = "../src/CURLOptTable.h")

              #---------------------

  # CURLINFO enumerations.

vals = edefs$CURLINFO@values
vals = vals[ - match(c("CURLINFO_NONE", "CURLINFO_LASTONE"), names(vals))]
cat(sprintf("#ifdef HAVE_%s\nCURLINFO(%s)%s\n#endif",  
                names(vals), gsub("CURLINFO_", "", names(vals)),
                c(rep(",", length(vals) -1), "")), sep = "\n",
                file = "../src/CURLINFOTable.h")
#cat(paste("CURLINFO(", gsub("CURLINFO_", "", names(vals)), ")", "\n", collapse = ",\n"),
#        file = "../src/CURLINFOTable.h")


               #----------------
con = file("../R/curlEnums.R", "w")
cat("if(!is.null(getClassDef('EnumValue'))) {\n", file = con)
enumNames = c("curl_infotype", "CURLcode", "curl_proxytype", "curl_usessl",
              "curl_ftpccc", "curl_ftpauth", "curl_ftpcreatedir", "curl_ftpmethod",
              "CURL_NETRC_OPTION", 
              "CURLFORMcode",
              "curl_TimeCond",
             "curl_closepolicy",
            )

if(any(!(enumNames %in% names(edefs))))
  stop("can't match ", paste(enumNames[!(enumNames %in% names(edefs))] , collapse = ", "))
# Ignored for now: CURLformoption

# anonymous  CURL_HTTP_VERSION, CURL_SSLVERSION

sapply(edefs[enumNames], writeCode, "r", file = con)


sapply(edefs[c("6368", "6604")], 
           function(x) 
            cat(paste(gsub("CURL_", "", names(x@values)), " = ", x@values, "L", sep = "", collapse = "\n"), "\n\n", file = con))


cat("\n}\n", file = con)
close(con)


 # Generate the documentation for the coercion methods.
sprintf("\\alias{%s-class}", enumNames)
cat(paste(sapply(c("integer", "numeric", "character"), function(x)  sprintf("\\alias{coerce,%s,%s-method}", x, enumNames)), collapse = "\n"))



 # CURLPROTO


i = grep("CURLAUTH", names(defs$macros))
auth = defs$macros[i]

BitwiseValue()
