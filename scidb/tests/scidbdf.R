# scidb array object tests
# Set the SCIDB_TEST_HOST system environment variable to the hostname or IP
# address of SciDB to run these tests. If SCIDB_TEST_HOST is not set the tests
# are not run.

check = function(a,b)
{
  print(match.call())
  stopifnot(all.equal(a,b,check.attributes=FALSE,check.names=FALSE))
}

library("scidb")
host = Sys.getenv("SCIDB_TEST_HOST")
if(nchar(host)>0)
{
  scidbconnect(host)
  options(scidb.debug=TRUE)

# Upload
  data("iris")
  x = as.scidb(iris)
# Factor levels are lost
  i = iris
  i[,5] = as.character(i[,5])
  check(x[],i)

# cast
  cast(x,x)

# Selection by name
  check(x[,"Petal_Length"][], iris[,"Petal.Length"])
  check(x$Petal_Length[], iris$Petal.Length)
# Selection along rows
  check(x[1:5,"Petal_Length"][], iris[1:5,"Petal.Length"])

# Unique of a single-attribute array
  u = unique(x$Species)
  check(count(u), 3)

# Aggregation by a non-integer attribute with a project thrown in
  check(aggregate(iris$Petal.Length,by=list(iris$Species),FUN=mean)[,2],
        aggregate(project(x,c('Petal_Length','Species')), by = 'Species', FUN='avg(Petal_Length)')[][,1])

# Conversion tests
 # bool          logical
 # char          character
 # datetime      double (aka real, numeric)
 # datetimetz    double
 # float         double
 # double        double
 # int64         double
 # uint64        double
 # uint32        double
 # int8          integer
 # uint8         integer
 # int16         integer
 # uint16        integer
 # int32         integer
 # string        character
  x = scidb("apply(build(<b:bool>[i=0:0,1,0],true),
             c, char('c'),
             dt, datetime(i),
             dz, datetime(i),
             f, float(i),
             d, double(i),
             i64, int64(i),
             ui64, uint64(i),
             ui32, uint32(i),
             i8, int8(i),
             i16, int16(i),
             ui16, uint16(i),
             x, string(i))")
  y = x[]


# $ indexing
 x = build("random()%100",20,type="double",start=1)
 x = bind(x,"w",2)
 f = function(a)  unique(a$val)
 u = f(x)

}
gc()
