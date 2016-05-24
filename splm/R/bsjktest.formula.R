`bsjktest.formula` <-
function(x, data, index=NULL, listw, test=c("C.1","C.2","C.3","J"), ...){

  ## transform listw if needed
  if("listw" %in% class(listw)) {
    w <- listw2mat(listw)
  } else {
    w <- listw
  }

  ## transform data if needed
  if(!is.null(index)) {
    #require(plm)
    data <- plm.data(data, index)
    }

  gindex <- dimnames(data)[[2]][1]
  tindex <- dimnames(data)[[2]][2]

  switch(match.arg(test), C.1 = {

    bsjk = pbsjkSDtest(formula=x, data=data, w=w, index=index, ...)

  }, C.2 = {

    bsjk = pbsjkARtest(formula=x, data=data, w=w, index=index, ...)

  }, C.3 = {

    stop("C.3 test not yet available")
    #bsjk = pbsjkREtest(formula=x, data=data, w=w, index=index, ...)

  }, J = {

    bsjk = pbsjkJtest(formula=x, data=data, w=w, index=index, ...)

  })

  return(bsjk)

}

