db_env_create <- function(flags=0L) {
  .Call("rberkeley_db_env_create", as.integer(flags))
}

db_get_env <- function(dbh)
{
  .Call("rberkeley_db_get_env", dbh)
}

dbenv_remove <- function(dbenv, db_home, flags=0L)
{
  .Call("rberkeley_dbenv_remove", dbenv, db_home, as.integer(flags))
}

dbenv_dbremove <- function(dbenv, txnid=NULL, file, database, flags)
{
  .Call("rberkeley_dbenv_dbremove", dbenv, txnid, file, database, as.integer(flags))
}

dbenv_dbrename <- function(dbenv, txnid=NULL, file, database, newname, flags)
{
  .Call("rberkeley_dbenv_dbrename", dbenv, txnid, file, database, 
        newname, as.integer(flags))
}

dbenv_set_cachesize <- function(dbenv, gbytes, bytes, ncache)
{
  .Call("rberkeley_dbenv_set_cachesize", 
        dbenv, 
        as.integer(gbytes),
        as.integer(bytes),
        as.integer(ncache))
}

dbenv_get_cachesize <- function(dbenv)
{
  cachesize <- .Call("rberkeley_dbenv_get_cachesize", dbenv)
  names(cachesize) <- c("gbytes","bytes","ncache")
  cachesize
}

dbenv_set_data_dir <- function(dbenv, dir)
{
  .Call("rberkeley_dbenv_set_data_dir", dbenv, dir)
}

dbenv_get_data_dirs <- function(dbenv)
{
  .Call("rberkeley_dbenv_get_data_dirs", dbenv)
}

dbenv_set_tmp_dir <- function(dbenv, dir)
{
  .Call("rberkeley_dbenv_set_tmp_dir", dbenv, dir)
}

dbenv_get_tmp_dir <- function(dbenv)
{
  .Call("rberkeley_dbenv_get_tmp_dir", dbenv)
}

dbenv_set_flags <- function(dbenv, flags, onoff)
{
  .Call("rberkeley_dbenv_set_flags", dbenv, flags, onoff)
}

dbenv_get_flags <- function(dbenv)
{
  .Call("rberkeley_dbenv_get_flags", dbenv)
}

dbenv_set_intermediate_dir_mode <- function(dbenv, mode)
{
  .Call("rberkeley_dbenv_set_intermediate_dir_mode", dbenv, mode)
}

dbenv_get_intermediate_dir_mode <- function(dbenv)
{
  .Call("rberkeley_dbenv_get_intermediate_dir_mode", dbenv)
}

dbenv_set_shm_key <- function(dbenv, shm_key)
{
  .Call("rberkeley_dbenv_set_shm_key", dbenv, as.integer(shm_key))
}

dbenv_get_shm_key <- function(dbenv)
{
  .Call("rberkeley_dbenv_get_shm_key", dbenv)
}

dbenv_open <- function(dbenv, db_home, flags, mode)
{
  .Call("rberkeley_dbenv_open", dbenv, db_home, flags, mode)
}

dbenv_get_home <- function(dbenv)
{
  if(check_pointer(dbenv))
  .Call("rberkeley_dbenv_get_home", dbenv)
}

dbenv_get_open_flags <- function(dbenv)
{
  if(check_pointer(dbenv))
  .Call("rberkeley_dbenv_get_open_flags", dbenv)
}

dbenv_close <- function(dbenv, flags=0L)
{
  if(check_pointer(dbenv))
    .Call("rberkeley_dbenv_close", dbenv, flags)
}

dbenv_stat_print <- function(dbenv, flags)
{
  if(check_pointer(dbenv))
    .Call("rberkeley_dbenv_stat_print", dbenv, flags)
}

dbenv_set_verbose <- function(dbenv, which, onoff=1L)
{
  .Call("rberkeley_dbenv_set_verbose", dbenv, as.integer(which), as.integer(onoff))
}

dbenv_get_verbose <- function(dbenv, which)
{
  if(missing(which))
    stop("'which' must be a valid flag set with mkFlags")
  .Call("rberkeley_dbenv_get_verbose", dbenv, as.integer(which))
}
