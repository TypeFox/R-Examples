dbenv_mutex_alloc <- function(dbenv, flags)
{
  .Call("rberkeley_dbenv_mutex_alloc", dbenv, flags)
}

dbenv_mutex_free <- function(dbenv, mutex)
{
  .Call("rberkeley_dbenv_mutex_free", dbenv, mutex)
}

dbenv_mutex_stat_print <- function(dbenv, flags)
{
  .Call("rberkeley_dbenv_mutex_stat_print", dbenv, as.integer(flags))
}

dbenv_mutex_lock <- function(dbenv, mutex)
{
  .Call("rberkeley_dbenv_mutex_lock", dbenv, mutex)
}

dbenv_mutex_unlock <- function(dbenv, mutex)
{
  .Call("rberkeley_dbenv_mutex_unlock", dbenv, mutex)
}
