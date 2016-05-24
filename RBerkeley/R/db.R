DBT <- function(data=NULL,size=NULL,ulen=NULL,dlen=NULL,doff=NULL,flags=NULL) {
  structure(list(data,size,ulen,dlen,doff,flags),class="DBT")
}

db_create <- function(dbenv=NULL, flags=0L) {
  structure(0, conn_id=.Call("rberkeley_db_create", dbenv, flags), class="DB")
}

db_strerror <- function(error)
{
  .Call("rberkeley_db_strerror", as.integer(error))
} 

db_del <- function(dbh, txnid=NULL, key, flags=0L)
{
  if(!is.raw(key))
    key <- serialize(key, NULL)
  .Call("rberkeley_db_del", as.DB(dbh), txnid, key, as.integer(flags))
}

db_open <- function(dbh, txnid=NULL, file="access.db",
                    database=NULL, type="BTREE", flags=0L) {
  if(!is.numeric(type)) {
  type <- switch(gsub("DB_","",toupper(type)), 
                       BTREE  =1L,
                       HASH   =2L,
                       RECNO  =3L,
                       QUEUE  =4L,
                       UNKNOWN=5L)
  if(is.null(type)) stop("'type' must be a one of BTREE, HASH, RECNO, QUEUE, or UNKNOWN")
  }
  .Call("rberkeley_db_open", as.DB(dbh), txnid, file, database, as.integer(type), flags)
}

db_close <- function(dbh) {
  .Call("rberkeley_db_close", as.DB(dbh))
}

db_upgrade <- function(dbh, file, flags)
{
  .Call("rberkeley_db_upgrade", as.DB(dbh), file, as.integer(flags))
}

db_compact <- function(dbh, txnid=NULL, start=NULL, stop=NULL,
                       c_data=NULL, flags=0L)
{
  .Call("rberkeley_db_compact", as.DB(dbh), txnid, start, stop, c_data, as.integer(flags))
}

db_set_priority <- function(dbh, priority)
{
  if(is.character(priority))
    priority <- mkFlags(priority)
  .Call("rberkeley_db_set_priority", as.DB(dbh), as.integer(priority))
}

db_get_priority <- function(dbh)
{
  .Call("rberkeley_db_get_priority", as.DB(dbh))
}

db_set_flags <- function(dbh, flags=0L)
{
  .Call("rberkeley_db_set_flags", as.DB(dbh), as.integer(flags))
}

db_get_flags <- function(dbh)
{
  .Call("rberkeley_db_get_flags", as.DB(dbh))
}

db_put <- function(dbh, txnid=NULL, key, data, flags=0L)
{
  if(!is.raw(key))
    key <- serialize(key, NULL)
  if(!is.raw(data))
    data <- serialize(data, NULL)

  invisible(.Call("rberkeley_db_put", as.DB(dbh), txnid, key, data, as.integer(flags)))
}

db_get <- function(dbh, txnid=NULL, key, data=NULL, flags=0L)
{
  if(!is.raw(key))
    key <- serialize(key, NULL)
  if(!is.null(data) && !is.raw(data))
    data <- serialize(data, NULL)

  # should add error checking here... 
  .Call("rberkeley_db_get", as.DB(dbh), txnid, key, data, as.integer(flags))
}

db_getP <- function(dbh, txnid=NULL, key, data=NULL, flags=0L)
{
  if(!is.raw(key))
    key <- serialize(key, NULL)

  if(!inherits(data,"DBT")) {
    if(!is.null(data) && !is.raw(data))
      data <- serialize(data,NULL)
    data <- DBT(data)
  }

  # should add error checking here... 
  .Call("rberkeley_db_getP", as.DB(dbh), txnid, key, data, as.integer(flags))
}

db_key_range <- function(dbh, txnid=NULL, key, flags=0L)
{
  if(!is.raw(key))
    key <- serialize(key, NULL)
  .Call("rberkeley_db_key_range", as.DB(dbh), txnid, key, as.integer(flags))
}

db_exists <- function(dbh, txnid=NULL, key, flags=0L)
{
  if(!is.raw(key))
    key <- serialize(key, NULL)
  .Call("rberkeley_db_exists", as.DB(dbh), txnid, key, flags)
}

db_truncate <- function(dhh)
{

}
db_get_byteswapped <- function(dbh)
{
  .Call("rberkeley_db_get_byteswapped", dbh)
}

db_set_cachesize <- function(dbh, gbytes, bytes, ncache)
{
  .Call("rberkeley_db_set_cachesize", 
        as.DB(dbh), 
        as.integer(gbytes),
        as.integer(bytes),
        as.integer(ncache))
}

db_get_cachesize <- function(dbh)
{
  cachesize <- .Call("rberkeley_db_get_cachesize", as.DB(dbh))
  names(cachesize) <- c("gbytes","bytes","ncache")
  cachesize
}

db_set_pagesize <- function(dbh, pagesize)
{
  .Call("rberkeley_db_set_pagesize", as.DB(dbh), as.integer(pagesize))
}

db_get_pagesize <- function(dbh)
{
  pagesize <- .Call("rberkeley_db_get_pagesize", as.DB(dbh))
  pagesize
}

db_set_encrypt <- function(dbh, passwd, flags)
{
  if(missing(flags))
    flags = mkFlags("DB_ENCRYPT_AES")
  
  .Call("rberkeley_db_set_encrypt", as.DB(dbh), passwd, flags)
}

db_get_encrypt_flags <- function(dbh)
{
  .Call("rberkeley_db_get_encrypt_flags", as.DB(dbh))
}
 
db_set_lorder <- function(dbh, lorder)
{
  lorder <- as.integer(lorder)
  if(lorder != 1234L || lorder != 4321L)
    stop("incorrect 'lorder' value")

  .Call("rberkeley_db_set_lorder", as.DB(dbh), lorder)
}

db_get_lorder <- function(dbh)
{
  lorder <- .Call("rberkeley_db_get_lorder", as.DB(dbh))
  if(lorder == 1234L)
   return("little.endian")
  if(lorder == 4321L)
   return("big.endian")
}

db_stat <- function(dbh, txnid=NULL, flags=mkFlags("DB_FAST_STAT"))
{
  .Call("rberkeley_db_stat", as.DB(dbh), txnid, as.integer(flags)) 
}

db_stat_print <- function(dbh, flags=0L)
{
  .Call("rberkeley_db_stat_print", as.DB(dbh), as.integer(flags)) 
}

db_remove <- function(dbh, file, database)
{
  .Call("rberkeley_db_remove", as.DB(dbh), file, database)
}

db_rename <- function(dbh, file, database, newname)
{
  .Call("rberkeley_db_rename", as.DB(dbh), file, database, newname)
}

db_version <- function() 
{
  .Call("rberkeley_db_version")
}

db_get_dbname <- function(dbh)
{
  .Call("rberkeley_db_get_dbname", as.DB(dbh))
}

db_sync <- function(dbh)
{
  .Call("rberkeley_db_sync", as.DB(dbh))
}

db_get_type <- function(dbh)
{
  .Call("rberkeley_db_get_type", as.DB(dbh))
}

db_set_errpfx <- function(dbh, errpfx)
{
  .Call("rberkeley_db_set_errfx", as.DB(dbh), errpfx)
}

db_get_errpfx <- function(dbh)
{
  .Call("rberkeley_db_get_errfx", as.DB(dbh))
}

db_set_errfile <- function(dbh, errfile)
{
  if(!is.null(errfile) && !file.exists(errfile)) {
    if(!file.create(errfile)) stop(paste("could not create file"))
  }
  .Call("rberkeley_db_set_errfile", dbh, errfile)
}

db_set_msgfile <- function(dbh, msgfile)
{
  if(!is.null(msgfile) && !file.exists(msgfile)) {
    if(!file.create(msgfile)) stop(paste("could not create file"))
  }
  .Call("rberkeley_db_set_msgfile", dbh, msgfile)
}

db_set_re_source <- function(dbh, source)
{
	if(!file.exists(source)) {
		stop(paste("could not set_re_source source to file",source))
	}
	.Call("rberkeley_db_set_re_source", 
			as.DB(dbh), 
			as.character(source))
}

