RS.connect <- function(host=NULL, port=6311L, tls=FALSE, proxy.target=NULL, proxy.wait=TRUE) .Call("RS_connect", host, port, tls, proxy.target, proxy.wait, PACKAGE="RSclient")

RS.close <- function(rsc) .Call("RS_close", rsc)

RS.eval <- function(rsc, x, wait=TRUE, lazy=TRUE) { r <- .Call("RS_eval", rsc, serialize(if (isTRUE(lazy)) substitute(x) else x, NULL, FALSE), wait, PACKAGE="RSclient"); if (is.raw(r)) unserialize(r) else r }

RS.eval.qap <- function(rsc, x, wait=TRUE) .Call("RS_eval_qap", rsc, x, wait, PACKAGE="RSclient")

RS.collect <- function(rsc, timeout = Inf, detail = FALSE) {
  r <- .Call("RS_collect", rsc, timeout, PACKAGE="RSclient")
  if (is.raw(r)) {
    if (length(r)) {
      if (isTRUE(detail))
        list(value = unserialize(r), rsc = attr(r, "rsc"))
      else unserialize(r)
    } else if (isTRUE(detail))
      list(rsc = attr(r, "rsc"))
    else NULL
  } else r
}

RS.server.eval <- function(rsc, text) .Call("RS_ctrl_str", rsc, 0x42L, text, PACKAGE="RSclient")

RS.server.source <- function(rsc, filename) .Call("RS_ctrl_str", rsc, 0x45L, filename, PACKAGE="RSclient")

RS.server.shutdown <- function(rsc) .Call("RS_ctrl_str", rsc, 0x44L, "", PACKAGE="RSclient")

RS.switch <- function(rsc, protocol="TLS") .Call("RS_switch", rsc, protocol, PACKAGE="RSclient")

RS.authkey <- function(rsc, type="rsa-authkey") .Call("RS_authkey", rsc, type, PACKAGE="RSclient")

RS.assign <- function(rsc, name, value, wait = TRUE) {
  if (missing(value)) {
    sym.name <- deparse(substitute(name))
    value <- name
    name <- sym.name
  }
  .Call("RS_assign", rsc, serialize(list(name, value), NULL), wait, PACKAGE="RSclient")
}

RS.login <- function(rsc, user, password, pubkey, authkey) {
  if (missing(user) || missing(password)) stop("user and password must be specified")
  .Call("RS_secauth", rsc, paste(c(user, password, ''), collapse="\n"), authkey, PACKAGE="RSclient")
}

RS.oobCallbacks <- function(rsc, send, msg) {
  if (missing(send) && missing(msg)) return(.Call("RS_oob_cb", rsc, NULL, NULL, TRUE))
  if (missing(send) || missing(msg)) {
    l <- .Call("RS_oob_cb", rsc, NULL, NULL, TRUE, PACKAGE="RSclient")
    if (missing(send)) send <- l$send
    if (missing(msg))  msg <- l$msg
  }
  invisible(.Call("RS_oob_cb", rsc, send, msg, FALSE, PACKAGE="RSclient"))  
}

print.RserveConnection <- function(x, ...) invisible(.Call("RS_print", x))
`==.RserveConnection` <- function(e1, e2) .Call("RS_eq", e1, e2)
`!=.RserveConnection` <- function(e1, e2) !.Call("RS_eq", e1, e2)
