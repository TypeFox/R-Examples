get.x.est <- function(z,model,p,par2,group){
  x.est <- z
  numGroups <- length(table(group))
  for (i in 1:numGroups){
    # verschiedene x.est je nach gruppe
    pt <- p.truth(model,p=p,par2=par2,group=i)
    y.mu <- y.mean(model,p=p,par2=par2,group=i)
    sel <- group == i
    x.est[sel] <- (z[sel] - (1-pt)*y.mu)/pt  
  }
  x.est
}

## in coming call for selection: get.quotient(X[sel,i],models[i],p.list[[i]], par2, group=1)
# => response z within group 1 etc
get.quotient <- function(z,model,p,par2, group, pi, piVar){
  pt <- p.truth(model,p,par2,group)
  # second, irrelevant distribution y
  y.mu <- y.mean(model,p,par2,group)
  y.var <- y.variance(y.mu,model,p,par2)
  
  # sensitive attribute x
#   if (model =="FR"){
#     x.est <- get.x.est (z,model,p,par2,rep(group,length(z)) )
#     x.mu <- mean(x.est)
#     x.var <- (var(z)-pt*(1-pt)*(x.mu-y.mu)^2 - (1-pt)*y.var)/pt
#     print("in cor.models")
#     print(x.mu)
#     print(x.var)
#   }else{
    # sensitive attribute x: take mean and variance from both groups!!!
    x.mu <- pi
    x.var <- piVar
    # NOT VARIANCE OF x.est !!!!! => VARIANCE OF TRUE VALUES NEEDED !!! 
#   }
   
  # noise (in paper u): mask.mu = 0
  # mask=noise term u (random error)
  mask.var <- (1-pt)/pt*(x.var+y.var/pt+(x.mu - y.mu)^2)
  quotient <- mask.var/x.var
}
  
# probability of answering to sensitive question
p.truth <- function(model,p,par2,group){
  switch(model,
         "direct" = pt <- 1 ,
         "Warner" =    pt <- 2*p-1, # ifelse(2*p-1<0, 1-2*p,2*p-1),
         "Crosswise" = pt <- 2*p-1, # ifelse(2*p-1<0, 1-2*p,2*p-1),
         "Mangat" = pt <- p,
         "Kuk" =  pt <- p[1]-p[2],
         "FR" = pt <- 1-sum(p) ,
         
         "UQTknown" = pt <- p[1],
         "SLD" = pt <-  par2 -1 +p[group], 
#          "CDM" = pt <- 1-p[group],
         "UQTunknown" = pt <- p[group],
         "mix.norm" = pt <- p[1],
         "mix.exp" = pt <- p[1],
         "mix.unknown" = pt <- p[group]
#          "CDMsym" = {
#            idx <- 2*group-1   # idx=1 für G1 ; idx=3 für G2
#            pt <- 1- p[idx]-p[idx+1]  # 1-p1-p2 für gruppe 1
#          }
         )
  pt
}

# mean of unrelated question
y.mean <- function(model,p,par2,group){
  switch(model,
         "SLD" = y.mu <- (1-p[group])/(2-par2-p[group]),
         "UQTunknown" = y.mu <- par2,
#          "CDM" = y.mu <- 1-par2,
#          "CDMsym" = {
#            idx <- 2*group-1   # idx=1 für G1 ; idx=3 für G2
#            y.mu <- (p[idx]*(1-par2)) / (p[idx]+p[idx+1])
#          },
         "direct" = y.mu <- 0 ,
         "Warner" =    y.mu <-  1/2,
         "Crosswise" = y.mu <-  1/2,
         "Mangat" = y.mu <- 1,
         "Kuk" = y.mu <- p[2]/(1-p[1]+p[2]),
         "FR" = y.mu <- ifelse(all(p==0),0, 1/sum(p)*p %*% 0:(length(p)-1)),
         "UQTknown" = y.mu <- p[2],
         "mix.norm" = y.mu <- p[2],
         "mix.exp" = y.mu <- p[2],
         "mix.unknown" = y.mu <- par2[1]
         )
  y.mu
}

# variance of unrelated question
y.variance <- function(y.mu,model,p,par2=NULL){
    if (model =="FR"){
      y.var <- ifelse(all(p==0),0,1/sum(p)*p %*% (0:(length(p)-1))^2 - 
                          (1/sum(p)*p %*% (0:(length(p)-1)))^2)
    }else if (model == "mix.norm"){
      y.var <- p[3]^2
    }else if (model == "mix.exp"){
      y.var <- p[2]^2
    }else if (model == "mix.unknown"){
      y.var <- par2[2]
    }else{
      y.var <- y.mu* (1-y.mu)
    }
    y.var
}

# reverse sign of correlation for certain combinatins of model and p
# if done for both variables, nothing happens
r.sign <- function(r,mod1,p1,mod2,p2){
  reverse1 <- ifelse (mod1 %in% c("Warner","Crosswise") && p1<.5,TRUE,FALSE)
  reverse2 <- ifelse (mod2 %in% c("Warner","Crosswise") && p2<.5,TRUE,FALSE) #|| (mod2=="Kuk" && p2[1]<p2[2]),
  r <- ifelse(reverse1+reverse2 == 1, -r, r)
  return(r)
}
