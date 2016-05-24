dbenv_txn_stat_print <- function(dbenv, flags=0L)
{
  .Call("rberkeley_dbenv_txn_stat_print", dbenv, as.integer(flags))
}

dbenv_txn_begin <- function(dbenv, parent=NULL, flags=0L)
{
  .Call("rberkeley_dbenv_txn_begin", dbenv, parent, as.integer(flags))
}

dbtxn_abort <- function(tid)
{
  .Call("rberkeley_dbtxn_abort", tid)
}

dbtxn_commit <- function(tid, flags)
{
  .Call("rberkeley_dbtxn_commit", tid, as.integer(flags))
}

dbtxn_id <- function(tid)
{
  .Call("rberkeley_dbtxn_id", tid)
}
