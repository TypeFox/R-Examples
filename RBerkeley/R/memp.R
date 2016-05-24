db_get_mpf <- function(dbh)
{
  .Call("rberkeley_db_get_mpf", dbh)
}

dbenv_memp_stat_print <- function(dbenv, flags=0L)
{
  .Call("rberkeley_dbenv_stat_print", dbenv, as.integer(flags))
}
