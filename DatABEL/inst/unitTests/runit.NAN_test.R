### --- Test setup ---

if(FALSE) {
  ## Not really needed, but can be handy when writing tests
  library("RUnit")
  library("DatABEL")
}

#test.empty <- function(){}
### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

source(paste(path,"/shared_functions.R",sep=""))

quiet <- TRUE

### --- Test functions ---

NA_export_import <- function(TYPE="DOUBLE",tolmult=10) 
{
    unlink("*.fv*");
    unlink("tmp*")

    v <- make_matrix_of_type_with_NA(TYPE,0.1)
    
    f <- matrix2databel(v,"tmp",64,TYPE);
    databel2text(f,"tmp_text","NAN");
    f2 <-  text2databel("tmp_text","tmp3",R_matrix=TRUE,type=TYPE,readonly=TRUE, naString="NAN");
    m2 <- as(f2,"matrix");
    checkEqualsNumeric(m2,v,tolerance=tolmult*sqrt(.Machine$double.eps));
    
    rm(f,f2);gc()
    
    unlink("*.fv*");
    unlink("tmp*")    
}

# make matrix, convert to fv, back and compare
NA_test <- function(TYPE)
{
    unlink("*.fv*");
    unlink("tmp*")
    
    v <- make_matrix_of_type_with_NA(TYPE,0.1)
    f <- matrix2databel(v, "tmp", 64, "DOUBLE");
    v1 <- as(f,"matrix")
    checkEqualsNumeric(v,v1);
    
    rm(f);gc()
	
    unlink("*.fv*");
    unlink("tmp*");
}

test.NA_string <- function(){
  unlink("*.fv*");
  unlink("tmp*")

  v=matrix(c(1,2,4,123),2,2);
  f=matrix2databel(v,"tmp",64,"INT");
  databel2text(f,"tmp2","NAN");
  f2 = text2databel("tmp2","tmp3",skipcols=1,skiprows=1,
      R_matrix=TRUE,type="INT",readonly=TRUE, naString="123");
  v2 <- v
  v2[v2=="123"] <- NA
  m2=as(f2,"matrix");
  checkEquals(m2,v2)
  checkEquals(m2[1,1],1);
  checkEquals(m2[2,1],2);
  checkEquals(m2[1,2],4);
  checkTrue(is.na(m2[2,2]));

  rm(f,f2);gc()
  
  unlink("*.fv*");
  unlink("tmp*")
}

test.NA_test <- function()
{
	NA_test("DOUBLE");
	NA_test("FLOAT");
	NA_test("INT");
	NA_test("UNSIGNED_INT");
	NA_test("SHORT_INT");
	NA_test("UNSIGNED_SHORT_INT");
	NA_test("CHAR");
	NA_test("UNSIGNED_CHAR");
}

test.NA_export_import <- function()
{
    NA_export_import("DOUBLE");
    NA_export_import("FLOAT");
    NA_export_import("INT");
    NA_export_import("UNSIGNED_INT");
    NA_export_import("SHORT_INT");
    NA_export_import("UNSIGNED_SHORT_INT");
    NA_export_import("CHAR");
    NA_export_import("UNSIGNED_CHAR");
}