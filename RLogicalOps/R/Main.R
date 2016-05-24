#' @title Processing Logical Operations
#' @description Dynamically process the logical operations over a list of string and get the result
#' @param logicalOperation
#' A String containing the logical operation to be performed
#' @param stringList
#' A list of String over which the logical operation is to be performed
#' @return A list of the output values(0/1)
#' @examples
#' logicalOperation <- "AND(Cake,OR(Birthday,Anniversary))"
#' stringList <- c("The cake at the birthday was awesome."
#'                 ,"Their anniversary was last week.")
#'
#' processLogiOper(logicalOperation, stringList)
#' @usage processLogiOper(logicalOperation, stringList)
#' @export
processLogiOper <- function(logicalOperation, stringList){
  res <- unlist(lapply(1:length(stringList), function(i){
    if(i %% 1000 == 0) print(i)
    testStack <- rstackdeque::rstack()
    out <- ""
    logOper <- logicalOperation

    while(stringr::str_detect(logOper,"\\(|\\)|,")){       #Encounter opening bracket/ closing bracket or comma
      end <- stringr::str_locate(logOper,"\\(|\\)|,")[1]   #Get first index of opening bracket/ closing bracket or comma

      if(stringr::str_sub(logOper,end,end) == "("){

        testStack <- rstackdeque::insert_top(testStack, stringr::str_sub(logOper, 1, end-1))    #Stack inserted upto first opening bracket
        testStack <- rstackdeque::insert_top(testStack, stringr::str_sub(logOper, end, end))    #Stack inserted with '('
        logOper <- stringr::str_sub(logOper, end+1, nchar(logOper))                #logOper updated with remaining text

      } else if(stringr::str_sub(logOper,end,end) == ")"){
        word1 <- stringr::str_sub(logOper, 1, end-1)               #Word till closing bracket
        if(stringr::str_detect(word1,"@")){                        #If @ encountered in word
          foundVarList <- word1
          while(rstackdeque::peek_top(testStack) != "("){              #Pop top element of stack till opening bracket is encountered
            foundVarList <- c(rstackdeque::peek_top(testStack), foundVarList)    #Element between opening and closing brackets
            testStack <- rstackdeque::without_top(testStack)
          }
          foundVarList <- gsub("@","",foundVarList)       #Supress @
        }else{                                            #If @ not encountered
          foundVarList <- ifelse(stringr::str_detect(stringList[i], stringr::fixed(word1, ignore_case = T)),1,0)    #Check if word is present in String
          while(rstackdeque::peek_top(testStack) != "("){                        #Text between opening and closing brackets with occurence in String
            foundVarList <- c(rstackdeque::peek_top(testStack), foundVarList)
            testStack <- rstackdeque::without_top(testStack)
          }
        }
        testStack <- rstackdeque::without_top(testStack)           #Remove top element from stack
        if(rstackdeque::peek_top(testStack) == "OR"){              #If top element is OR
          out <- 0
          foundVarList <- gsub("@","",foundVarList)
          for(j in 1:length(foundVarList)){
            out <- ifelse(out|as.numeric(foundVarList[j]),1,0)
          }
        }else if(rstackdeque::peek_top(testStack) == "NOT"){       #If top element is NOT
          foundVarList <- gsub("@","",foundVarList)
          out <- ifelse(as.numeric(foundVarList),0,1)
        }else if(rstackdeque::peek_top(testStack) == "AND"){       #If top element is AND
          out <- 1
          foundVarList <- gsub("@","",foundVarList)
          for(j in 1:length(foundVarList)){
            out <- ifelse(out&as.numeric(foundVarList[j]),1,0)
          }
        }
        testStack <- rstackdeque::without_top(testStack)
        logOper <- stringr::str_replace(logOper,paste0(stringr::str_sub(logOper, 1, end-1),"\\)"),paste0("@",out))
        rm(word1, foundVarList)
      }else {
        word2 <- stringr::str_sub(logOper, 1, end-1)
        if(stringr::str_detect(word2,"@")){
          testStack <- rstackdeque::insert_top(testStack,word2)
          logOper <- stringr::str_sub(logOper, end+1, nchar(logOper))
        }else{
          foundVar <- ifelse(stringr::str_detect(stringList[i], stringr::fixed(word2, ignore_case = T)),1,0)
          testStack <- rstackdeque::insert_top(testStack,foundVar)
          logOper <- stringr::str_sub(logOper, end+1, nchar(logOper))
        }
        rm(word2)
      }
    }
    return(out)
  }))
  return(res)
}
