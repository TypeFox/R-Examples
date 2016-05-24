# This function is a convenience wrapper to .gettext() which
# passes the domain name and avoid repeating it in the code
.Facto_gettext <- function(msg) {
    gettext(msg, domain="R-RcmdrPlugin.FactoMineR")
}

.Facto_ngettext <- function(n, msg1, msg2) {
    ngettext(n, msg1, msg2, domain="R-RcmdrPlugin.FactoMineR")
}

