# This function is a convenience wrapper to .gettext() which
# passes the domain name and avoid repeating it in the code
.gettext <- function(msg) {
    gettext(msg, domain="R-RcmdrPlugin.temis")
}

.ngettext <- function(n, msg1, msg2) {
    ngettext(n, msg1, msg2, domain="R-RcmdrPlugin.temis")
}

