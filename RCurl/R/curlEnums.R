if(!is.null(getClassDef('EnumValue'))) {
setClass('curl_infotype', contains = 'EnumValue')
`curl_infotypeValues` = EnumDef('curl_infotype', structure(as.integer(c(0,
1,
2,
3,
4,
5,
6,
7)),
names = c("CURLINFO_TEXT",
"CURLINFO_HEADER_IN",
"CURLINFO_HEADER_OUT",
"CURLINFO_DATA_IN",
"CURLINFO_DATA_OUT",
"CURLINFO_SSL_DATA_IN",
"CURLINFO_SSL_DATA_OUT",
"CURLINFO_END")))


setAs('numeric', 'curl_infotype', function(from)  asEnumValue(from, `curl_infotypeValues`, 'curl_infotype', prefix = "CURLINFO_"))
setAs('character', 'curl_infotype', function(from)  asEnumValue(from, `curl_infotypeValues`, 'curl_infotype', prefix = "CURLINFO_"))
setAs('integer', 'curl_infotype', function(from)  asEnumValue(from, `curl_infotypeValues`, 'curl_infotype', prefix = "CURLINFO_"))


`CURLINFO_TEXT` <- GenericEnumValue('CURLINFO_TEXT', 0, 'curl_infotype')
`CURLINFO_HEADER_IN` <- GenericEnumValue('CURLINFO_HEADER_IN', 1, 'curl_infotype')
`CURLINFO_HEADER_OUT` <- GenericEnumValue('CURLINFO_HEADER_OUT', 2, 'curl_infotype')
`CURLINFO_DATA_IN` <- GenericEnumValue('CURLINFO_DATA_IN', 3, 'curl_infotype')
`CURLINFO_DATA_OUT` <- GenericEnumValue('CURLINFO_DATA_OUT', 4, 'curl_infotype')
`CURLINFO_SSL_DATA_IN` <- GenericEnumValue('CURLINFO_SSL_DATA_IN', 5, 'curl_infotype')
`CURLINFO_SSL_DATA_OUT` <- GenericEnumValue('CURLINFO_SSL_DATA_OUT', 6, 'curl_infotype')
`CURLINFO_END` <- GenericEnumValue('CURLINFO_END', 7, 'curl_infotype')

#####################
setClass('CURLcode', contains = 'EnumValue')
`CURLcodeValues` = EnumDef('CURLcode', structure(as.integer(c(0,
1,
2,
3,
4,
5,
6,
7,
8,
9,
10,
11,
12,
13,
14,
15,
16,
17,
18,
19,
20,
21,
22,
23,
24,
25,
26,
27,
28,
29,
30,
31,
32,
33,
34,
35,
36,
37,
38,
39,
40,
41,
42,
43,
44,
45,
46,
47,
48,
49,
50,
51,
52,
53,
54,
55,
56,
57,
58,
59,
60,
61,
62,
63,
64,
65,
66,
67,
68,
69,
70,
71,
72,
73,
74,
75,
76,
77,
78,
79,
80,
81,
82,
83,
84)),
names = c("CURLE_OK",
"CURLE_UNSUPPORTED_PROTOCOL",
"CURLE_FAILED_INIT",
"CURLE_URL_MALFORMAT",
"CURLE_OBSOLETE4",
"CURLE_COULDNT_RESOLVE_PROXY",
"CURLE_COULDNT_RESOLVE_HOST",
"CURLE_COULDNT_CONNECT",
"CURLE_FTP_WEIRD_SERVER_REPLY",
"CURLE_REMOTE_ACCESS_DENIED",
"CURLE_OBSOLETE10",
"CURLE_FTP_WEIRD_PASS_REPLY",
"CURLE_OBSOLETE12",
"CURLE_FTP_WEIRD_PASV_REPLY",
"CURLE_FTP_WEIRD_227_FORMAT",
"CURLE_FTP_CANT_GET_HOST",
"CURLE_OBSOLETE16",
"CURLE_FTP_COULDNT_SET_TYPE",
"CURLE_PARTIAL_FILE",
"CURLE_FTP_COULDNT_RETR_FILE",
"CURLE_OBSOLETE20",
"CURLE_QUOTE_ERROR",
"CURLE_HTTP_RETURNED_ERROR",
"CURLE_WRITE_ERROR",
"CURLE_OBSOLETE24",
"CURLE_UPLOAD_FAILED",
"CURLE_READ_ERROR",
"CURLE_OUT_OF_MEMORY",
"CURLE_OPERATION_TIMEDOUT",
"CURLE_OBSOLETE29",
"CURLE_FTP_PORT_FAILED",
"CURLE_FTP_COULDNT_USE_REST",
"CURLE_OBSOLETE32",
"CURLE_RANGE_ERROR",
"CURLE_HTTP_POST_ERROR",
"CURLE_SSL_CONNECT_ERROR",
"CURLE_BAD_DOWNLOAD_RESUME",
"CURLE_FILE_COULDNT_READ_FILE",
"CURLE_LDAP_CANNOT_BIND",
"CURLE_LDAP_SEARCH_FAILED",
"CURLE_OBSOLETE40",
"CURLE_FUNCTION_NOT_FOUND",
"CURLE_ABORTED_BY_CALLBACK",
"CURLE_BAD_FUNCTION_ARGUMENT",
"CURLE_OBSOLETE44",
"CURLE_INTERFACE_FAILED",
"CURLE_OBSOLETE46",
"CURLE_TOO_MANY_REDIRECTS",
"CURLE_UNKNOWN_TELNET_OPTION",
"CURLE_TELNET_OPTION_SYNTAX",
"CURLE_OBSOLETE50",
"CURLE_PEER_FAILED_VERIFICATION",
"CURLE_GOT_NOTHING",
"CURLE_SSL_ENGINE_NOTFOUND",
"CURLE_SSL_ENGINE_SETFAILED",
"CURLE_SEND_ERROR",
"CURLE_RECV_ERROR",
"CURLE_OBSOLETE57",
"CURLE_SSL_CERTPROBLEM",
"CURLE_SSL_CIPHER",
"CURLE_SSL_CACERT",
"CURLE_BAD_CONTENT_ENCODING",
"CURLE_LDAP_INVALID_URL",
"CURLE_FILESIZE_EXCEEDED",
"CURLE_USE_SSL_FAILED",
"CURLE_SEND_FAIL_REWIND",
"CURLE_SSL_ENGINE_INITFAILED",
"CURLE_LOGIN_DENIED",
"CURLE_TFTP_NOTFOUND",
"CURLE_TFTP_PERM",
"CURLE_REMOTE_DISK_FULL",
"CURLE_TFTP_ILLEGAL",
"CURLE_TFTP_UNKNOWNID",
"CURLE_REMOTE_FILE_EXISTS",
"CURLE_TFTP_NOSUCHUSER",
"CURLE_CONV_FAILED",
"CURLE_CONV_REQD",
"CURLE_SSL_CACERT_BADFILE",
"CURLE_REMOTE_FILE_NOT_FOUND",
"CURLE_SSH",
"CURLE_SSL_SHUTDOWN_FAILED",
"CURLE_AGAIN",
"CURLE_SSL_CRL_BADFILE",
"CURLE_SSL_ISSUER_ERROR",
"CURL_LAST")))


setAs('numeric', 'CURLcode', function(from)  asEnumValue(from, `CURLcodeValues`, 'CURLcode', prefix = "CURL"))
setAs('character', 'CURLcode', function(from)  asEnumValue(from, `CURLcodeValues`, 'CURLcode', prefix = "CURL"))
setAs('integer', 'CURLcode', function(from)  asEnumValue(from, `CURLcodeValues`, 'CURLcode', prefix = "CURL"))


`CURLE_OK` <- GenericEnumValue('CURLE_OK', 0, 'CURLcode')
`CURLE_UNSUPPORTED_PROTOCOL` <- GenericEnumValue('CURLE_UNSUPPORTED_PROTOCOL', 1, 'CURLcode')
`CURLE_FAILED_INIT` <- GenericEnumValue('CURLE_FAILED_INIT', 2, 'CURLcode')
`CURLE_URL_MALFORMAT` <- GenericEnumValue('CURLE_URL_MALFORMAT', 3, 'CURLcode')
`CURLE_OBSOLETE4` <- GenericEnumValue('CURLE_OBSOLETE4', 4, 'CURLcode')
`CURLE_COULDNT_RESOLVE_PROXY` <- GenericEnumValue('CURLE_COULDNT_RESOLVE_PROXY', 5, 'CURLcode')
`CURLE_COULDNT_RESOLVE_HOST` <- GenericEnumValue('CURLE_COULDNT_RESOLVE_HOST', 6, 'CURLcode')
`CURLE_COULDNT_CONNECT` <- GenericEnumValue('CURLE_COULDNT_CONNECT', 7, 'CURLcode')
`CURLE_FTP_WEIRD_SERVER_REPLY` <- GenericEnumValue('CURLE_FTP_WEIRD_SERVER_REPLY', 8, 'CURLcode')
`CURLE_REMOTE_ACCESS_DENIED` <- GenericEnumValue('CURLE_REMOTE_ACCESS_DENIED', 9, 'CURLcode')
`CURLE_OBSOLETE10` <- GenericEnumValue('CURLE_OBSOLETE10', 10, 'CURLcode')
`CURLE_FTP_WEIRD_PASS_REPLY` <- GenericEnumValue('CURLE_FTP_WEIRD_PASS_REPLY', 11, 'CURLcode')
`CURLE_OBSOLETE12` <- GenericEnumValue('CURLE_OBSOLETE12', 12, 'CURLcode')
`CURLE_FTP_WEIRD_PASV_REPLY` <- GenericEnumValue('CURLE_FTP_WEIRD_PASV_REPLY', 13, 'CURLcode')
`CURLE_FTP_WEIRD_227_FORMAT` <- GenericEnumValue('CURLE_FTP_WEIRD_227_FORMAT', 14, 'CURLcode')
`CURLE_FTP_CANT_GET_HOST` <- GenericEnumValue('CURLE_FTP_CANT_GET_HOST', 15, 'CURLcode')
`CURLE_OBSOLETE16` <- GenericEnumValue('CURLE_OBSOLETE16', 16, 'CURLcode')
`CURLE_FTP_COULDNT_SET_TYPE` <- GenericEnumValue('CURLE_FTP_COULDNT_SET_TYPE', 17, 'CURLcode')
`CURLE_PARTIAL_FILE` <- GenericEnumValue('CURLE_PARTIAL_FILE', 18, 'CURLcode')
`CURLE_FTP_COULDNT_RETR_FILE` <- GenericEnumValue('CURLE_FTP_COULDNT_RETR_FILE', 19, 'CURLcode')
`CURLE_OBSOLETE20` <- GenericEnumValue('CURLE_OBSOLETE20', 20, 'CURLcode')
`CURLE_QUOTE_ERROR` <- GenericEnumValue('CURLE_QUOTE_ERROR', 21, 'CURLcode')
`CURLE_HTTP_RETURNED_ERROR` <- GenericEnumValue('CURLE_HTTP_RETURNED_ERROR', 22, 'CURLcode')
`CURLE_WRITE_ERROR` <- GenericEnumValue('CURLE_WRITE_ERROR', 23, 'CURLcode')
`CURLE_OBSOLETE24` <- GenericEnumValue('CURLE_OBSOLETE24', 24, 'CURLcode')
`CURLE_UPLOAD_FAILED` <- GenericEnumValue('CURLE_UPLOAD_FAILED', 25, 'CURLcode')
`CURLE_READ_ERROR` <- GenericEnumValue('CURLE_READ_ERROR', 26, 'CURLcode')
`CURLE_OUT_OF_MEMORY` <- GenericEnumValue('CURLE_OUT_OF_MEMORY', 27, 'CURLcode')
`CURLE_OPERATION_TIMEDOUT` <- GenericEnumValue('CURLE_OPERATION_TIMEDOUT', 28, 'CURLcode')
`CURLE_OBSOLETE29` <- GenericEnumValue('CURLE_OBSOLETE29', 29, 'CURLcode')
`CURLE_FTP_PORT_FAILED` <- GenericEnumValue('CURLE_FTP_PORT_FAILED', 30, 'CURLcode')
`CURLE_FTP_COULDNT_USE_REST` <- GenericEnumValue('CURLE_FTP_COULDNT_USE_REST', 31, 'CURLcode')
`CURLE_OBSOLETE32` <- GenericEnumValue('CURLE_OBSOLETE32', 32, 'CURLcode')
`CURLE_RANGE_ERROR` <- GenericEnumValue('CURLE_RANGE_ERROR', 33, 'CURLcode')
`CURLE_HTTP_POST_ERROR` <- GenericEnumValue('CURLE_HTTP_POST_ERROR', 34, 'CURLcode')
`CURLE_SSL_CONNECT_ERROR` <- GenericEnumValue('CURLE_SSL_CONNECT_ERROR', 35, 'CURLcode')
`CURLE_BAD_DOWNLOAD_RESUME` <- GenericEnumValue('CURLE_BAD_DOWNLOAD_RESUME', 36, 'CURLcode')
`CURLE_FILE_COULDNT_READ_FILE` <- GenericEnumValue('CURLE_FILE_COULDNT_READ_FILE', 37, 'CURLcode')
`CURLE_LDAP_CANNOT_BIND` <- GenericEnumValue('CURLE_LDAP_CANNOT_BIND', 38, 'CURLcode')
`CURLE_LDAP_SEARCH_FAILED` <- GenericEnumValue('CURLE_LDAP_SEARCH_FAILED', 39, 'CURLcode')
`CURLE_OBSOLETE40` <- GenericEnumValue('CURLE_OBSOLETE40', 40, 'CURLcode')
`CURLE_FUNCTION_NOT_FOUND` <- GenericEnumValue('CURLE_FUNCTION_NOT_FOUND', 41, 'CURLcode')
`CURLE_ABORTED_BY_CALLBACK` <- GenericEnumValue('CURLE_ABORTED_BY_CALLBACK', 42, 'CURLcode')
`CURLE_BAD_FUNCTION_ARGUMENT` <- GenericEnumValue('CURLE_BAD_FUNCTION_ARGUMENT', 43, 'CURLcode')
`CURLE_OBSOLETE44` <- GenericEnumValue('CURLE_OBSOLETE44', 44, 'CURLcode')
`CURLE_INTERFACE_FAILED` <- GenericEnumValue('CURLE_INTERFACE_FAILED', 45, 'CURLcode')
`CURLE_OBSOLETE46` <- GenericEnumValue('CURLE_OBSOLETE46', 46, 'CURLcode')
`CURLE_TOO_MANY_REDIRECTS` <- GenericEnumValue('CURLE_TOO_MANY_REDIRECTS', 47, 'CURLcode')
`CURLE_UNKNOWN_TELNET_OPTION` <- GenericEnumValue('CURLE_UNKNOWN_TELNET_OPTION', 48, 'CURLcode')
`CURLE_TELNET_OPTION_SYNTAX` <- GenericEnumValue('CURLE_TELNET_OPTION_SYNTAX', 49, 'CURLcode')
`CURLE_OBSOLETE50` <- GenericEnumValue('CURLE_OBSOLETE50', 50, 'CURLcode')
`CURLE_PEER_FAILED_VERIFICATION` <- GenericEnumValue('CURLE_PEER_FAILED_VERIFICATION', 51, 'CURLcode')
`CURLE_GOT_NOTHING` <- GenericEnumValue('CURLE_GOT_NOTHING', 52, 'CURLcode')
`CURLE_SSL_ENGINE_NOTFOUND` <- GenericEnumValue('CURLE_SSL_ENGINE_NOTFOUND', 53, 'CURLcode')
`CURLE_SSL_ENGINE_SETFAILED` <- GenericEnumValue('CURLE_SSL_ENGINE_SETFAILED', 54, 'CURLcode')
`CURLE_SEND_ERROR` <- GenericEnumValue('CURLE_SEND_ERROR', 55, 'CURLcode')
`CURLE_RECV_ERROR` <- GenericEnumValue('CURLE_RECV_ERROR', 56, 'CURLcode')
`CURLE_OBSOLETE57` <- GenericEnumValue('CURLE_OBSOLETE57', 57, 'CURLcode')
`CURLE_SSL_CERTPROBLEM` <- GenericEnumValue('CURLE_SSL_CERTPROBLEM', 58, 'CURLcode')
`CURLE_SSL_CIPHER` <- GenericEnumValue('CURLE_SSL_CIPHER', 59, 'CURLcode')
`CURLE_SSL_CACERT` <- GenericEnumValue('CURLE_SSL_CACERT', 60, 'CURLcode')
`CURLE_BAD_CONTENT_ENCODING` <- GenericEnumValue('CURLE_BAD_CONTENT_ENCODING', 61, 'CURLcode')
`CURLE_LDAP_INVALID_URL` <- GenericEnumValue('CURLE_LDAP_INVALID_URL', 62, 'CURLcode')
`CURLE_FILESIZE_EXCEEDED` <- GenericEnumValue('CURLE_FILESIZE_EXCEEDED', 63, 'CURLcode')
`CURLE_USE_SSL_FAILED` <- GenericEnumValue('CURLE_USE_SSL_FAILED', 64, 'CURLcode')
`CURLE_SEND_FAIL_REWIND` <- GenericEnumValue('CURLE_SEND_FAIL_REWIND', 65, 'CURLcode')
`CURLE_SSL_ENGINE_INITFAILED` <- GenericEnumValue('CURLE_SSL_ENGINE_INITFAILED', 66, 'CURLcode')
`CURLE_LOGIN_DENIED` <- GenericEnumValue('CURLE_LOGIN_DENIED', 67, 'CURLcode')
`CURLE_TFTP_NOTFOUND` <- GenericEnumValue('CURLE_TFTP_NOTFOUND', 68, 'CURLcode')
`CURLE_TFTP_PERM` <- GenericEnumValue('CURLE_TFTP_PERM', 69, 'CURLcode')
`CURLE_REMOTE_DISK_FULL` <- GenericEnumValue('CURLE_REMOTE_DISK_FULL', 70, 'CURLcode')
`CURLE_TFTP_ILLEGAL` <- GenericEnumValue('CURLE_TFTP_ILLEGAL', 71, 'CURLcode')
`CURLE_TFTP_UNKNOWNID` <- GenericEnumValue('CURLE_TFTP_UNKNOWNID', 72, 'CURLcode')
`CURLE_REMOTE_FILE_EXISTS` <- GenericEnumValue('CURLE_REMOTE_FILE_EXISTS', 73, 'CURLcode')
`CURLE_TFTP_NOSUCHUSER` <- GenericEnumValue('CURLE_TFTP_NOSUCHUSER', 74, 'CURLcode')
`CURLE_CONV_FAILED` <- GenericEnumValue('CURLE_CONV_FAILED', 75, 'CURLcode')
`CURLE_CONV_REQD` <- GenericEnumValue('CURLE_CONV_REQD', 76, 'CURLcode')
`CURLE_SSL_CACERT_BADFILE` <- GenericEnumValue('CURLE_SSL_CACERT_BADFILE', 77, 'CURLcode')
`CURLE_REMOTE_FILE_NOT_FOUND` <- GenericEnumValue('CURLE_REMOTE_FILE_NOT_FOUND', 78, 'CURLcode')
`CURLE_SSH` <- GenericEnumValue('CURLE_SSH', 79, 'CURLcode')
`CURLE_SSL_SHUTDOWN_FAILED` <- GenericEnumValue('CURLE_SSL_SHUTDOWN_FAILED', 80, 'CURLcode')
`CURLE_AGAIN` <- GenericEnumValue('CURLE_AGAIN', 81, 'CURLcode')
`CURLE_SSL_CRL_BADFILE` <- GenericEnumValue('CURLE_SSL_CRL_BADFILE', 82, 'CURLcode')
`CURLE_SSL_ISSUER_ERROR` <- GenericEnumValue('CURLE_SSL_ISSUER_ERROR', 83, 'CURLcode')
`CURL_LAST` <- GenericEnumValue('CURL_LAST', 84, 'CURLcode')

#####################
setClass('curl_proxytype', contains = 'EnumValue')
`curl_proxytypeValues` = EnumDef('curl_proxytype', structure(as.integer(c(0,
1,
4,
5,
6,
7)),
names = c("CURLPROXY_HTTP",
"CURLPROXY_HTTP_1_0",
"CURLPROXY_SOCKS4",
"CURLPROXY_SOCKS5",
"CURLPROXY_SOCKS4A",
"CURLPROXY_SOCKS5_HOSTNAME")))


setAs('numeric', 'curl_proxytype', function(from)  asEnumValue(from, `curl_proxytypeValues`, 'curl_proxytype', prefix = "CURLPROXY_"))
setAs('character', 'curl_proxytype', function(from)  asEnumValue(from, `curl_proxytypeValues`, 'curl_proxytype', prefix = "CURLPROXY_"))
setAs('integer', 'curl_proxytype', function(from)  asEnumValue(from, `curl_proxytypeValues`, 'curl_proxytype', prefix = "CURLPROXY_"))


`CURLPROXY_HTTP` <- GenericEnumValue('CURLPROXY_HTTP', 0, 'curl_proxytype')
`CURLPROXY_HTTP_1_0` <- GenericEnumValue('CURLPROXY_HTTP_1_0', 1, 'curl_proxytype')
`CURLPROXY_SOCKS4` <- GenericEnumValue('CURLPROXY_SOCKS4', 4, 'curl_proxytype')
`CURLPROXY_SOCKS5` <- GenericEnumValue('CURLPROXY_SOCKS5', 5, 'curl_proxytype')
`CURLPROXY_SOCKS4A` <- GenericEnumValue('CURLPROXY_SOCKS4A', 6, 'curl_proxytype')
`CURLPROXY_SOCKS5_HOSTNAME` <- GenericEnumValue('CURLPROXY_SOCKS5_HOSTNAME', 7, 'curl_proxytype')

#####################
setClass('curl_usessl', contains = 'EnumValue')
`curl_usesslValues` = EnumDef('curl_usessl', structure(as.integer(c(0,
1,
2,
3,
4)),
names = c("CURLUSESSL_NONE",
"CURLUSESSL_TRY",
"CURLUSESSL_CONTROL",
"CURLUSESSL_ALL",
"CURLUSESSL_LAST")))


setAs('numeric', 'curl_usessl', function(from)  asEnumValue(from, `curl_usesslValues`, 'curl_usessl', prefix = "CURLUSESSL_"))
setAs('character', 'curl_usessl', function(from)  asEnumValue(from, `curl_usesslValues`, 'curl_usessl', prefix = "CURLUSESSL_"))
setAs('integer', 'curl_usessl', function(from)  asEnumValue(from, `curl_usesslValues`, 'curl_usessl', prefix = "CURLUSESSL_"))


`CURLUSESSL_NONE` <- GenericEnumValue('CURLUSESSL_NONE', 0, 'curl_usessl')
`CURLUSESSL_TRY` <- GenericEnumValue('CURLUSESSL_TRY', 1, 'curl_usessl')
`CURLUSESSL_CONTROL` <- GenericEnumValue('CURLUSESSL_CONTROL', 2, 'curl_usessl')
`CURLUSESSL_ALL` <- GenericEnumValue('CURLUSESSL_ALL', 3, 'curl_usessl')
`CURLUSESSL_LAST` <- GenericEnumValue('CURLUSESSL_LAST', 4, 'curl_usessl')

#####################
setClass('curl_ftpccc', contains = 'EnumValue')
`curl_ftpcccValues` = EnumDef('curl_ftpccc', structure(as.integer(c(0,
1,
2,
3)),
names = c("CURLFTPSSL_CCC_NONE",
"CURLFTPSSL_CCC_PASSIVE",
"CURLFTPSSL_CCC_ACTIVE",
"CURLFTPSSL_CCC_LAST")))


setAs('numeric', 'curl_ftpccc', function(from)  asEnumValue(from, `curl_ftpcccValues`, 'curl_ftpccc', prefix = "CURLFTPSSL_CCC_"))
setAs('character', 'curl_ftpccc', function(from)  asEnumValue(from, `curl_ftpcccValues`, 'curl_ftpccc', prefix = "CURLFTPSSL_CCC_"))
setAs('integer', 'curl_ftpccc', function(from)  asEnumValue(from, `curl_ftpcccValues`, 'curl_ftpccc', prefix = "CURLFTPSSL_CCC_"))


`CURLFTPSSL_CCC_NONE` <- GenericEnumValue('CURLFTPSSL_CCC_NONE', 0, 'curl_ftpccc')
`CURLFTPSSL_CCC_PASSIVE` <- GenericEnumValue('CURLFTPSSL_CCC_PASSIVE', 1, 'curl_ftpccc')
`CURLFTPSSL_CCC_ACTIVE` <- GenericEnumValue('CURLFTPSSL_CCC_ACTIVE', 2, 'curl_ftpccc')
`CURLFTPSSL_CCC_LAST` <- GenericEnumValue('CURLFTPSSL_CCC_LAST', 3, 'curl_ftpccc')

#####################
setClass('curl_ftpauth', contains = 'EnumValue')
`curl_ftpauthValues` = EnumDef('curl_ftpauth', structure(as.integer(c(0,
1,
2,
3)),
names = c("CURLFTPAUTH_DEFAULT",
"CURLFTPAUTH_SSL",
"CURLFTPAUTH_TLS",
"CURLFTPAUTH_LAST")))


setAs('numeric', 'curl_ftpauth', function(from)  asEnumValue(from, `curl_ftpauthValues`, 'curl_ftpauth', prefix = "CURLFTPAUTH_"))
setAs('character', 'curl_ftpauth', function(from)  asEnumValue(from, `curl_ftpauthValues`, 'curl_ftpauth', prefix = "CURLFTPAUTH_"))
setAs('integer', 'curl_ftpauth', function(from)  asEnumValue(from, `curl_ftpauthValues`, 'curl_ftpauth', prefix = "CURLFTPAUTH_"))


`CURLFTPAUTH_DEFAULT` <- GenericEnumValue('CURLFTPAUTH_DEFAULT', 0, 'curl_ftpauth')
`CURLFTPAUTH_SSL` <- GenericEnumValue('CURLFTPAUTH_SSL', 1, 'curl_ftpauth')
`CURLFTPAUTH_TLS` <- GenericEnumValue('CURLFTPAUTH_TLS', 2, 'curl_ftpauth')
`CURLFTPAUTH_LAST` <- GenericEnumValue('CURLFTPAUTH_LAST', 3, 'curl_ftpauth')

#####################
setClass('curl_ftpcreatedir', contains = 'EnumValue')
`curl_ftpcreatedirValues` = EnumDef('curl_ftpcreatedir', structure(as.integer(c(0,
1,
2,
3)),
names = c("CURLFTP_CREATE_DIR_NONE",
"CURLFTP_CREATE_DIR",
"CURLFTP_CREATE_DIR_RETRY",
"CURLFTP_CREATE_DIR_LAST")))


setAs('numeric', 'curl_ftpcreatedir', function(from)  asEnumValue(from, `curl_ftpcreatedirValues`, 'curl_ftpcreatedir', prefix = "CURLFTP_CREATE_DIR"))
setAs('character', 'curl_ftpcreatedir', function(from)  asEnumValue(from, `curl_ftpcreatedirValues`, 'curl_ftpcreatedir', prefix = "CURLFTP_CREATE_DIR"))
setAs('integer', 'curl_ftpcreatedir', function(from)  asEnumValue(from, `curl_ftpcreatedirValues`, 'curl_ftpcreatedir', prefix = "CURLFTP_CREATE_DIR"))


`CURLFTP_CREATE_DIR_NONE` <- GenericEnumValue('CURLFTP_CREATE_DIR_NONE', 0, 'curl_ftpcreatedir')
`CURLFTP_CREATE_DIR` <- GenericEnumValue('CURLFTP_CREATE_DIR', 1, 'curl_ftpcreatedir')
`CURLFTP_CREATE_DIR_RETRY` <- GenericEnumValue('CURLFTP_CREATE_DIR_RETRY', 2, 'curl_ftpcreatedir')
`CURLFTP_CREATE_DIR_LAST` <- GenericEnumValue('CURLFTP_CREATE_DIR_LAST', 3, 'curl_ftpcreatedir')

#####################
setClass('curl_ftpmethod', contains = 'EnumValue')
`curl_ftpmethodValues` = EnumDef('curl_ftpmethod', structure(as.integer(c(0,
1,
2,
3,
4)),
names = c("CURLFTPMETHOD_DEFAULT",
"CURLFTPMETHOD_MULTICWD",
"CURLFTPMETHOD_NOCWD",
"CURLFTPMETHOD_SINGLECWD",
"CURLFTPMETHOD_LAST")))


setAs('numeric', 'curl_ftpmethod', function(from)  asEnumValue(from, `curl_ftpmethodValues`, 'curl_ftpmethod', prefix = "CURLFTPMETHOD_"))
setAs('character', 'curl_ftpmethod', function(from)  asEnumValue(from, `curl_ftpmethodValues`, 'curl_ftpmethod', prefix = "CURLFTPMETHOD_"))
setAs('integer', 'curl_ftpmethod', function(from)  asEnumValue(from, `curl_ftpmethodValues`, 'curl_ftpmethod', prefix = "CURLFTPMETHOD_"))


`CURLFTPMETHOD_DEFAULT` <- GenericEnumValue('CURLFTPMETHOD_DEFAULT', 0, 'curl_ftpmethod')
`CURLFTPMETHOD_MULTICWD` <- GenericEnumValue('CURLFTPMETHOD_MULTICWD', 1, 'curl_ftpmethod')
`CURLFTPMETHOD_NOCWD` <- GenericEnumValue('CURLFTPMETHOD_NOCWD', 2, 'curl_ftpmethod')
`CURLFTPMETHOD_SINGLECWD` <- GenericEnumValue('CURLFTPMETHOD_SINGLECWD', 3, 'curl_ftpmethod')
`CURLFTPMETHOD_LAST` <- GenericEnumValue('CURLFTPMETHOD_LAST', 4, 'curl_ftpmethod')

#####################
setClass('CURL_NETRC_OPTION', contains = 'EnumValue')
`CURL_NETRC_OPTIONValues` = EnumDef('CURL_NETRC_OPTION', structure(as.integer(c(0,
1,
2,
3)),
names = c("CURL_NETRC_IGNORED",
"CURL_NETRC_OPTIONAL",
"CURL_NETRC_REQUIRED",
"CURL_NETRC_LAST")))


setAs('numeric', 'CURL_NETRC_OPTION', function(from)  asEnumValue(from, `CURL_NETRC_OPTIONValues`, 'CURL_NETRC_OPTION', prefix = "CURL_NETRC_"))
setAs('character', 'CURL_NETRC_OPTION', function(from)  asEnumValue(from, `CURL_NETRC_OPTIONValues`, 'CURL_NETRC_OPTION', prefix = "CURL_NETRC_"))
setAs('integer', 'CURL_NETRC_OPTION', function(from)  asEnumValue(from, `CURL_NETRC_OPTIONValues`, 'CURL_NETRC_OPTION', prefix = "CURL_NETRC_"))


`CURL_NETRC_IGNORED` <- GenericEnumValue('CURL_NETRC_IGNORED', 0, 'CURL_NETRC_OPTION')
`CURL_NETRC_OPTIONAL` <- GenericEnumValue('CURL_NETRC_OPTIONAL', 1, 'CURL_NETRC_OPTION')
`CURL_NETRC_REQUIRED` <- GenericEnumValue('CURL_NETRC_REQUIRED', 2, 'CURL_NETRC_OPTION')
`CURL_NETRC_LAST` <- GenericEnumValue('CURL_NETRC_LAST', 3, 'CURL_NETRC_OPTION')

#####################
setClass('CURLFORMcode', contains = 'EnumValue')
`CURLFORMcodeValues` = EnumDef('CURLFORMcode', structure(as.integer(c(0,
1,
2,
3,
4,
5,
6,
7,
8)),
names = c("CURL_FORMADD_OK",
"CURL_FORMADD_MEMORY",
"CURL_FORMADD_OPTION_TWICE",
"CURL_FORMADD_NULL",
"CURL_FORMADD_UNKNOWN_OPTION",
"CURL_FORMADD_INCOMPLETE",
"CURL_FORMADD_ILLEGAL_ARRAY",
"CURL_FORMADD_DISABLED",
"CURL_FORMADD_LAST")))


setAs('numeric', 'CURLFORMcode', function(from)  asEnumValue(from, `CURLFORMcodeValues`, 'CURLFORMcode', prefix = "CURL_FORMADD_"))
setAs('character', 'CURLFORMcode', function(from)  asEnumValue(from, `CURLFORMcodeValues`, 'CURLFORMcode', prefix = "CURL_FORMADD_"))
setAs('integer', 'CURLFORMcode', function(from)  asEnumValue(from, `CURLFORMcodeValues`, 'CURLFORMcode', prefix = "CURL_FORMADD_"))


`CURL_FORMADD_OK` <- GenericEnumValue('CURL_FORMADD_OK', 0, 'CURLFORMcode')
`CURL_FORMADD_MEMORY` <- GenericEnumValue('CURL_FORMADD_MEMORY', 1, 'CURLFORMcode')
`CURL_FORMADD_OPTION_TWICE` <- GenericEnumValue('CURL_FORMADD_OPTION_TWICE', 2, 'CURLFORMcode')
`CURL_FORMADD_NULL` <- GenericEnumValue('CURL_FORMADD_NULL', 3, 'CURLFORMcode')
`CURL_FORMADD_UNKNOWN_OPTION` <- GenericEnumValue('CURL_FORMADD_UNKNOWN_OPTION', 4, 'CURLFORMcode')
`CURL_FORMADD_INCOMPLETE` <- GenericEnumValue('CURL_FORMADD_INCOMPLETE', 5, 'CURLFORMcode')
`CURL_FORMADD_ILLEGAL_ARRAY` <- GenericEnumValue('CURL_FORMADD_ILLEGAL_ARRAY', 6, 'CURLFORMcode')
`CURL_FORMADD_DISABLED` <- GenericEnumValue('CURL_FORMADD_DISABLED', 7, 'CURLFORMcode')
`CURL_FORMADD_LAST` <- GenericEnumValue('CURL_FORMADD_LAST', 8, 'CURLFORMcode')

#####################
setClass('curl_TimeCond', contains = 'EnumValue')
`curl_TimeCondValues` = EnumDef('curl_TimeCond', structure(as.integer(c(0,
1,
2,
3,
4)),
names = c("CURL_TIMECOND_NONE",
"CURL_TIMECOND_IFMODSINCE",
"CURL_TIMECOND_IFUNMODSINCE",
"CURL_TIMECOND_LASTMOD",
"CURL_TIMECOND_LAST")))


setAs('numeric', 'curl_TimeCond', function(from)  asEnumValue(from, `curl_TimeCondValues`, 'curl_TimeCond', prefix = "CURL_TIMECOND_"))
setAs('character', 'curl_TimeCond', function(from)  asEnumValue(from, `curl_TimeCondValues`, 'curl_TimeCond', prefix = "CURL_TIMECOND_"))
setAs('integer', 'curl_TimeCond', function(from)  asEnumValue(from, `curl_TimeCondValues`, 'curl_TimeCond', prefix = "CURL_TIMECOND_"))


`CURL_TIMECOND_NONE` <- GenericEnumValue('CURL_TIMECOND_NONE', 0, 'curl_TimeCond')
`CURL_TIMECOND_IFMODSINCE` <- GenericEnumValue('CURL_TIMECOND_IFMODSINCE', 1, 'curl_TimeCond')
`CURL_TIMECOND_IFUNMODSINCE` <- GenericEnumValue('CURL_TIMECOND_IFUNMODSINCE', 2, 'curl_TimeCond')
`CURL_TIMECOND_LASTMOD` <- GenericEnumValue('CURL_TIMECOND_LASTMOD', 3, 'curl_TimeCond')
`CURL_TIMECOND_LAST` <- GenericEnumValue('CURL_TIMECOND_LAST', 4, 'curl_TimeCond')

#####################
setClass('curl_closepolicy', contains = 'EnumValue')
`curl_closepolicyValues` = EnumDef('curl_closepolicy', structure(as.integer(c(0,
1,
2,
3,
4,
5,
6)),
names = c("CURLCLOSEPOLICY_NONE",
"CURLCLOSEPOLICY_OLDEST",
"CURLCLOSEPOLICY_LEAST_RECENTLY_USED",
"CURLCLOSEPOLICY_LEAST_TRAFFIC",
"CURLCLOSEPOLICY_SLOWEST",
"CURLCLOSEPOLICY_CALLBACK",
"CURLCLOSEPOLICY_LAST")))


setAs('numeric', 'curl_closepolicy', function(from)  asEnumValue(from, `curl_closepolicyValues`, 'curl_closepolicy', prefix = "CURLCLOSEPOLICY_"))
setAs('character', 'curl_closepolicy', function(from)  asEnumValue(from, `curl_closepolicyValues`, 'curl_closepolicy', prefix = "CURLCLOSEPOLICY_"))
setAs('integer', 'curl_closepolicy', function(from)  asEnumValue(from, `curl_closepolicyValues`, 'curl_closepolicy', prefix = "CURLCLOSEPOLICY_"))


`CURLCLOSEPOLICY_NONE` <- GenericEnumValue('CURLCLOSEPOLICY_NONE', 0, 'curl_closepolicy')
`CURLCLOSEPOLICY_OLDEST` <- GenericEnumValue('CURLCLOSEPOLICY_OLDEST', 1, 'curl_closepolicy')
`CURLCLOSEPOLICY_LEAST_RECENTLY_USED` <- GenericEnumValue('CURLCLOSEPOLICY_LEAST_RECENTLY_USED', 2, 'curl_closepolicy')
`CURLCLOSEPOLICY_LEAST_TRAFFIC` <- GenericEnumValue('CURLCLOSEPOLICY_LEAST_TRAFFIC', 3, 'curl_closepolicy')
`CURLCLOSEPOLICY_SLOWEST` <- GenericEnumValue('CURLCLOSEPOLICY_SLOWEST', 4, 'curl_closepolicy')
`CURLCLOSEPOLICY_CALLBACK` <- GenericEnumValue('CURLCLOSEPOLICY_CALLBACK', 5, 'curl_closepolicy')
`CURLCLOSEPOLICY_LAST` <- GenericEnumValue('CURLCLOSEPOLICY_LAST', 6, 'curl_closepolicy')

#####################
HTTP_VERSION_NONE = 0L
HTTP_VERSION_1_0 = 1L
HTTP_VERSION_1_1 = 2L
HTTP_VERSION_LAST = 3L 

SSLVERSION_DEFAULT = 0L
SSLVERSION_TLSv1 = 1L
SSLVERSION_SSLv2 = 2L
SSLVERSION_SSLv3 = 3L
SSLVERSION_LAST = 4L 


}
