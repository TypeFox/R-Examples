`CEF` <-
function (typ=NA, fun=NA, dis=NA, a2=NA, c=NA, p1=NA, p2=p1) {
  if (!is.na(typ)) {
    if (!is.na(a2) && is.na(c) && is.na(p1) && is.na(p2)) {
      result <- a2.cef(typ,a2)
      }
    else {
      if (is.na(a2) && !is.na(c) && is.na(p1) && is.na(p2)) {
        result <- c.cef(typ,c)
        }
      else {
        if (is.na(a2) && is.na(c) && !is.na(p1) && !is.na(p2)) {
          result <- p1p2.cef(typ,p1,p2)
          }
          else {
            result <- NA
            }
        }
      }
    }
  else {
    if (!is.na(dis)) {
      if (!is.na(a2)) {
        result <- distort.a2(dis, fun, a2)
        }
      else {
        result <- distort(dis, fun, p1, p2)
        }
      }
    else {
      result <- fun
      }
    }
  return(result)
  }

