#' Cost of goods sold and ending inventory under three methods (FIFO,LIFO,Weighted average)
#'
#' @param uinv units of beginning inventory
#' @param pinv prince of beginning inventory
#' @param units nx1 vector of inventory units. inventory purchased ordered by time (from first to last)
#' @param price nx1 vector of inventory price. same order as units
#' @param sinv units of sold inventory
#' @param method inventory methods: FIFO (first in first out, permitted under both US and IFRS), LIFO (late in first out, US only), WAC (weighted average cost,US and IFRS)
#' @export
#' @examples
#' cogs(uinv=2,pinv=2,units=c(3,5),price=c(3,5),sinv=7,method="FIFO")
#'
#' cogs(uinv=2,pinv=2,units=c(3,5),price=c(3,5),sinv=7,method="LIFO")
#'
#' cogs(uinv=2,pinv=2,units=c(3,5),price=c(3,5),sinv=7,method="WAC")
cogs <- function(uinv,pinv,units,price,sinv,method="FIFO"){
  n=length(units)
  m=length(price)
  costOfGoods=0
  endingInventory=0
  if(m!=n){
    stop("length of units and price are not the same\n")
  }else{
    if(method=="FIFO"){
      if(sinv <= uinv){
        costOfGoods = sinv * pinv
        endingInventory = (uinv - sinv) * pinv
        for(i in 1:n){
          endingInventory = endingInventory + units[i] * price[i]
        }
      }else{
        costOfGoods = uinv * pinv
        sinv = sinv - uinv
        for(i in 1:n){
          if(sinv <= units[i]){
            costOfGoods = costOfGoods +sinv * price[i]
            endingInventory = (units[i] - sinv) * price[i]
            if(i < n){
              temp=i + 1
              for(j in temp:n){
                endingInventory = endingInventory + units[j] * price[j]
              }
            }
            sinv=0
            next
          }else{
            costOfGoods = costOfGoods +units[i] * price[i]
            sinv = sinv -units[i]
          }
        }
        if(sinv >0){
          stop("Inventory is not enough to sell\n")
        }
      }
    }else if(method=="WAC"){
      endingInventory = uinv * pinv
      tu = uinv
      for(i in 1:n){
        endingInventory = endingInventory + units[i] * price[i]
        tu = tu + units[i]
      }
      if(tu >= sinv){
        costOfGoods = endingInventory/tu*sinv
        endingInventory = endingInventory/tu*(tu-sinv)
      }else{
        stop("Inventory is not enough to sell\n")
      }
      
    }else if(method=="LIFO"){
      for(i in n:1){
        if(sinv <= units[i]){
          costOfGoods = costOfGoods +sinv * price[i]
          endingInventory = (units[i] - sinv) * price[i]
          if(i > 1){
            temp=i - 1
            for(j in temp:1){
              endingInventory = endingInventory + units[j] * price[j]
            }
          }
          endingInventory = endingInventory + uinv * pinv
          sinv=0
          next
        }else{
          costOfGoods = costOfGoods +units[i] * price[i]
          sinv = sinv -units[i]
        }
      }
      if(sinv > 0){
        if(sinv <= uinv){
          costOfGoods = costOfGoods + sinv * pinv
          endingInventory = endingInventory + (uinv - sinv) * pinv
        }
        else{
          stop("Inventory is not enough to sell\n")
        }
      }
    }
    
  }
  return(list("costOfGoods"=costOfGoods,"endingInventory"=endingInventory))
}
