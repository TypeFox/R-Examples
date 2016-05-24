`association.fit` <-
function (var, dep, adj, quantitative, type, level, nIndiv, genotypingRate=0, ...) 
{

    if (!quantitative) {
      if (length(unique(dep))==1)
      {
       res <- "Genot error"
      } 
      else 
       {
        co<-dom<-co<-dom<-rec<-over<-Ad<-NULL
        dep <- as.factor(dep)
        controlGeno <- ifelse(is.null(levels(var)),0,(length(var)/nIndiv)*100)
        if (genotypingRate >= controlGeno)
          {
            res <- c(paste("Genot ", round(controlGeno, 1), "\\%", sep = ""))
          }

        else if (length(table(as.character(var)))==1) {
            res <- "Monomorphic"
        }
        else {
         if (length(table(as.character(var))) == 3) {
                var.co <- codominant(var)
              if (any(type%in%6) | any(type%in%2))
                var.dom <- dominant(var)
              if (any(type%in%6) | any(type%in%3)) 
                var.rec <- recessive(var)
              if (any(type%in%6) | any(type%in%4))
                var.over <- overdominant(var)
            if (is.null(adj)) {
                  m.co <- glm(dep ~ var.co, family = binomial, ...)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ NULL, subset = subset, family = binomial, ...)
              if (any(type%in%6) | any(type%in%2))
                  m.dom <- glm(dep ~ var.dom, subset = subset, 
                    family = binomial, ...)
              if (any(type%in%6) | any(type%in%3))
                  m.rec <- glm(dep ~ var.rec, subset = subset, 
                    family = binomial, ...)
              if (any(type%in%6) | any(type%in%4))
                  m.over <- glm(dep ~ var.over, subset = subset, 
                    family = binomial, ...)
              if (any(type%in%6) | any(type%in%5))
                  m.ad <- glm(dep ~ as.numeric(var.co), subset = subset, 
                    family = binomial, ...)
                }
            else {
                  m.co <- glm(dep ~ . + var.co, family = binomial, 
                    data = adj, ...)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ ., subset = subset, family = binomial, 
                    data = adj, ...)
              if (any(type%in%6) | any(type%in%2))
                  m.dom <- glm(dep ~ . + var.dom, subset = subset, 
                    family = binomial, data = adj, ...)
              if (any(type%in%6) | any(type%in%3))
                  m.rec <- glm(dep ~ . + var.rec, subset = subset, 
                    family = binomial, data = adj, ...)
              if (any(type%in%6) | any(type%in%4))
                  m.over <- glm(dep ~ . + var.over, subset = subset, 
                    family = binomial, data = adj, ...)
              if (any(type%in%6) | any(type%in%5))
                  m.ad <- glm(dep ~ . + as.numeric(var.co), subset = subset, 
                    family = binomial, data = adj, ...)
                }


              if (any(type%in%6) | any(type%in%1))
               {
                temptp<-Table.N.Per(var.co, dep, subset)$tp
                co <- cbind(temptp, 
                  intervals.or(m.co, level, m.b, var)$or.ic, 
                  c(round(AIC(m.co), 1), NA, NA))
                if (any(temptp == 0) & is.null(adj)) 
                    {

                      pp<-fisher.test(dep,var.co)$p
                      co[1, 8] <- pp
                    }
               }
              if (any(type%in%6) | any(type%in%2))
               {
                temptp<-Table.N.Per(var.dom, dep, subset)$tp
                dom <- cbind(temptp, 
                  intervals.or(m.dom, level, m.b, var.dom)$or.ic, 
                  c(round(AIC(m.dom), 1), NA))
                if (any(temptp == 0) & is.null(adj)) 
                    {
                      pp<-fisher.test(dep,var.dom)$p
                      dom[1, 8] <- pp
                    }
               } 
              if (any(type%in%6) | any(type%in%3))
               {
                temptp<-Table.N.Per(var.rec, dep, subset)$tp
                rec <- cbind(temptp, 
                  intervals.or(m.rec, level, m.b, var.rec)$or.ic, 
                  c(round(AIC(m.rec), 1), NA))
                if (any(temptp == 0) & is.null(adj)) 
                    {
                      pp<-fisher.test(dep,var.rec)$p
                      rec[1, 8] <- pp
                    }
               }
              if (any(type%in%6) | any(type%in%4))
               {
                temptp<-Table.N.Per(var.over, dep, subset)$tp
                over <- cbind(temptp, 
                  intervals.or(m.over, level, m.b, var.over)$or.ic, 
                  c(round(AIC(m.over), 1), NA))
                if (any(temptp == 0) & is.null(adj)) 
                    {
                      pp<-fisher.test(dep,var.over)$p
                      over[1, 8] <- pp
                    }
               }
              if (any(type%in%6) | any(type%in%5))
               {
                temptp<-Table.N.Per(var.co, dep, subset)$tp
                totals<-round(table(dep),1)
                prop.totals<-round(100*prop.table(totals),1)
                ansTot<-c(totals[1],prop.totals[1], totals[2],prop.totals[2])
                Ad <- c(ansTot, intervals.or(m.ad, 
                  level, m.b)$or.ic, round(AIC(m.ad), 1))
                if (any(temptp == 0) & is.null(adj)) 
                    {
                      pp<-fisher.test(dep,var.co)$p
                      Ad[8] <- pp
                    }
               } 
             res<-NULL 
             if(!is.null(co))
               res<-rbind(Codominant=rep(NA,9),co)
             if(!is.null(dom))
               res<-rbind(res,Dominant=rep(NA,9),dom)  
             if(!is.null(rec))
               res<-rbind(res,Recessive=rep(NA,9),rec)
             if(!is.null(over))
               res<-rbind(res,Overdominant=rep(NA,9),over)
             if(!is.null(Ad))
               res<-rbind(res,"log-Additive"=rep(NA,9),"0,1,2"=Ad)

             dimnames(res)[[2]][5:9] <- c("OR","lower","upper","p-value","AIC")
               
            }
            else if (length(table(as.character(var))) == 2) {
                var.co <- codominant(var)
                if (is.null(adj)) {
                  m.co <- glm(dep ~ var.co, family = binomial, ...)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ NULL, subset = subset, family = binomial, ...)
                }
                else {
                  m.co <- glm(dep ~ . + var.co, family = binomial, 
                    data = adj, ...)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ ., subset = subset, family = binomial, 
                    data = adj, ...)
                }
                co <- cbind(Table.N.Per(var.co, dep, subset)$tp, 
                  intervals.or(m.co, level, m.b, var.co)$or.ic, 
                  c(round(AIC(m.co), 1), NA))
#                Ad <- c(rep(NA, times = 4), intervals.or(m.ad, level, m.b)$or.ic, round(AIC(m.ad), 1))

                totals<-table(dep)
                prop.totals<-round(100*prop.table(totals),1)
                ansTot<-c(totals[1],prop.totals[1], totals[2],prop.totals[2])
                Ad <- c(ansTot, intervals.or(m.co, level, m.b)$or.ic, round(AIC(m.co), 1))
                  
                Ad[8]<-NA 

               if(any(Table.N.Per(var.co, dep, subset)$tp==0) & is.null(adj))
                {
                 pp<-fisher.test(dep,var.co)$p
                 Ad[8]<-pp
                 co[1,8]<-pp 
                }

                res <- rbind(Codominant=rep(NA,9),co,"log-Additive"=rep(NA,9), "0,1,2"=Ad)
                dimnames(res)[[2]][5:9] <- c("OR","lower","upper","p-value","AIC")
            }
        }
       }
    }
    else {    # quantitative trait
        co<-dom<-co<-dom<-rec<-over<-Ad<-NULL 
        controlGeno <- ifelse(is.null(levels(var)),0,(length(var)/nIndiv)*100)
        if (genotypingRate >= controlGeno)
          {
            res <- c(paste("Genot ", round(controlGeno, 1), "\\%", sep = ""))
          }

        else if (length(table(as.character(var)))==1) {
            res <- "Monomorphic"
        }
        else {
            if (length(table(as.character(var))) == 3) {
              var.co <- codominant(var)
              if (any(type%in%6) | any(type%in%2))
                var.dom <- dominant(var)
              if (any(type%in%6) | any(type%in%3))
                var.rec <- recessive(var)
              if (any(type%in%6) | any(type%in%4))
                var.over <- overdominant(var)
                if (is.null(adj)) {
                  m.co <- glm(dep ~ var.co, family = gaussian, ...)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ NULL, subset = subset, family = gaussian, ...)
              if (any(type%in%6) | any(type%in%2))
                  m.dom <- glm(dep ~ var.dom, subset = subset, 
                    family = gaussian, ...)
              if (any(type%in%6) | any(type%in%3))
                  m.rec <- glm(dep ~ var.rec, subset = subset, 
                    family = gaussian, ...)
              if (any(type%in%6) | any(type%in%4))
                  m.over <- glm(dep ~ var.over, subset = subset, 
                    family = gaussian, ...)
              if (any(type%in%6) | any(type%in%5))
                  m.ad <- glm(dep ~ as.numeric(var.co), subset = subset, 
                    family = gaussian, ...)
                }
                else {
                  m.co <- glm(dep ~ . + var.co, family = gaussian, 
                    data = adj, ...)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ ., subset = subset, family = gaussian, 
                    data = adj, ...)
              if (any(type%in%6) | any(type%in%2))
                  m.dom <- glm(dep ~ . + var.dom, subset = subset, 
                    family = gaussian, data = adj, ...)
              if (any(type%in%6) | any(type%in%3))
                  m.rec <- glm(dep ~ . + var.rec, subset = subset, 
                    family = gaussian, data = adj, ...)
              if (any(type%in%6) | any(type%in%4))
                  m.over <- glm(dep ~ . + var.over, subset = subset, 
                    family = gaussian, data = adj, ...)
              if (any(type%in%6) | any(type%in%5))
                  m.ad <- glm(dep ~ . + as.numeric(var.co), subset = subset, 
                    family = gaussian, data = adj, ...)
                }
              if (any(type%in%6) | any(type%in%1))
                co <- cbind(Table.mean.se(var.co, dep, subset)$tp, 
                  intervals.dif(m.co, level, m.b, var)$m, AIC = c(round(AIC(m.co), 
                    1), NA, NA))
              if (any(type%in%6) | any(type%in%2))
                dom <- cbind(Table.mean.se(var.dom, dep, subset)$tp, 
                  intervals.dif(m.dom, level, m.b, var.dom)$m, 
                  AIC = c(round(AIC(m.dom), 1), NA))
              if (any(type%in%6) | any(type%in%3))
                rec <- cbind(Table.mean.se(var.rec, dep, subset)$tp, 
                  intervals.dif(m.rec, level, m.b, var.rec)$m, 
                  AIC = c(round(AIC(m.rec), 1), NA))
              if (any(type%in%6) | any(type%in%4))
                over <- cbind(Table.mean.se(var.over, dep, subset)$tp, 
                  intervals.dif(m.over, level, m.b, var.over)$m, 
                  AIC = c(round(AIC(m.over), 1), NA))
              if (any(type%in%6) | any(type%in%5))
                Ad <- c(rep(NA, 3), intervals.dif(m.ad, level, 
                  m.b)$m, AIC(m.ad))

             res<-NULL 
             if(!is.null(co))
               res<-rbind(Codominant=rep(NA,8),co)
             if(!is.null(dom))
               res<-rbind(res,Dominant=rep(NA,8),dom)  
             if(!is.null(rec))
               res<-rbind(res,Recessive=rep(NA,8),rec)
             if(!is.null(over))
               res<-rbind(res,Overdominant=rep(NA,8),over)
             if(!is.null(Ad))
               res<-rbind(res,"log-Additive"=rep(NA,8),"0,1,2"=Ad)
              dimnames(res)[[2]][4:8] <- c("dif","lower","upper","p-value","AIC")
            }


           else if (length(table(as.character(var))) == 2) {
                var.co <- codominant(var)
                if (is.null(adj)) {
                  m.co <- glm(dep ~ var.co, family = gaussian, ...)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ NULL, subset = subset, family = gaussian, ...)
                }
                else {
                  m.co <- glm(dep ~ . + var.co, family = gaussian, 
                    data = adj, ...)
                  subset <- 1:length(var) %in% as.numeric(rownames(m.co$model))
                  m.b <- glm(dep ~ ., subset = subset, family = gaussian, 
                    data = adj, ...)
                }
                co <- cbind(Table.mean.se(var.co, dep, subset)$tp, 
                  intervals.dif(m.co, level, m.b, var.co)$m, AIC = c(AIC(m.co), 
                    NA))
#                Ad <- c(rep(NA, 3), intervals.dif(m.ad, level, m.b)$m, round(AIC(m.ad), 1))
                Ad <- c(rep(NA, 3), intervals.dif(m.co, level,m.b)$m, round(AIC(m.co), 1)) 

                Ad[7]<-NA 

                if(any(Table.mean.se(var.co, dep, subset)$tp==0) & is.null(adj))
                 {
                  pp<-fisher.test(dep,var.co)$p
                  Ad[7]<-pp
                  co[1,7]<-pp 
                 }

                res <- rbind(Codominant=rep(NA,8),co,"log-Additive"=rep(NA,8), "0,1,2"=Ad)
                dimnames(res)[[2]][4:8] <- c("dif","lower","upper","p-value","AIC")
            }
        }
    }
    res
}


