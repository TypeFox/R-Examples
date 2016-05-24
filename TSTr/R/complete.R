complete <- function(tree,string,node=1,match = "") {
  b<-tree
  guardo<-c()
  if (nchar(string) > 0) {
    stav <- str_sub(string,1,1)
    if (stav < b$ch[node]) {
      if (is.na(b$L[node])) {
        print("No match found")
      }
      else {
        guardo<-append(guardo,complete(b,string,b$L[node],match))
      }
    }
    else if (stav > b$ch[node]) {
      if (is.na(b$R[node])) {
        print("No match found")
      }
      else {
        guardo<-append(guardo,complete(b,string,b$R[node],match))
      }
    }
    else {
      if (nchar(string) == 1) {
        if (!is.na(node) && b$flag[node] == 1) {
          guardo<-append(guardo,paste0(match,b$ch[node]))
          guardo<-append(guardo,spdfs(b,b$C[node],paste0(match,b$ch[node],b$ch[b$C[node]])))
        }
        else if (!is.na(b$C[node])) {
          guardo<-append(guardo,spdfs(b,b$C[node],paste0(match,b$ch[node],b$ch[b$C[node]])))
        }
        else {
          print(1)
        }
      }
      else {
        guardo<-append(guardo,complete(b,str_sub(string,2,-1),b$C[node],paste0(match,stav)))
      }
    }
  }
  else {
    print("Invalid string")
  }
  return(guardo)
}