#####################################################################################
#                      Functions for Bilder & Loughin (2007)                        #
#####################################################################################

# data.format()
# A function that reformats the raw data or bootstrap resample data into model form

data.format<-function(data, I, J, K, nvars, add.constant = .5, model.vars = NULL,
                      predict.func = FALSE) {
  nrows<-(2^nvars)*I*J*max(1,K)
  if (predict.func) {
    if (MRCV_globals$print.status) {
      MRCV_globals$pb.counter<-MRCV_globals$pb.counter+1.7
      setTxtProgressBar(MRCV_globals$pb, MRCV_globals$pb.counter)
    }
  }
  if (is.null(model.vars)) {
    model.data.unsorted<-data.frame(matrix(data = NA, nrow = nrows, ncol = (2*nvars+1)))
    counter<-0
    if (nvars==2) {
      for (i in 1:I) {
        for (j in 1:J) {
          counter<-counter+2^nvars
          table.count<-table(data[,i],data[,(I+j)])
          model.data.unsorted[(counter-3),]<-c(names(data)[i], names(data)[(I+j)],0,0, 
                                               as.numeric(table.count[1,1]))
          model.data.unsorted[(counter-2),]<-c(names(data)[i], names(data)[(I+j)],0,1, 
                                               as.numeric(table.count[1,2]))
          model.data.unsorted[(counter-1),]<-c(names(data)[i], names(data)[(I+j)],1,0, 
                                               as.numeric(table.count[2,1]))
          model.data.unsorted[(counter),]  <-c(names(data)[i], names(data)[(I+j)],1,1, 
                                               as.numeric(table.count[2,2]))
        }
      }
      colnames(model.data.unsorted)<-c("W", "Y", "wi", "yj", "count")
    }
    if (nvars==3) {
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K) {
            counter<-counter+2^nvars
            table.count<-table(data[,i],data[,(I+j)],data[,(I+J+k)])
            model.data.unsorted[(counter-7),]<-c(names(data)[i], names(data)[(I+j)],names(data)[(I+J+k)],0,0,0, 
                                                 as.numeric(table.count[1,1,1]))
            model.data.unsorted[(counter-6),]<-c(names(data)[i], names(data)[(I+j)],names(data)[(I+J+k)],0,0,1, 
                                                 as.numeric(table.count[1,1,2]))
            model.data.unsorted[(counter-5),]<-c(names(data)[i], names(data)[(I+j)],names(data)[(I+J+k)],0,1,0, 
                                                 as.numeric(table.count[1,2,1]))
            model.data.unsorted[(counter-4),]<-c(names(data)[i], names(data)[(I+j)],names(data)[(I+J+k)],0,1,1, 
                                                 as.numeric(table.count[1,2,2]))
            model.data.unsorted[(counter-3),]<-c(names(data)[i], names(data)[(I+j)],names(data)[(I+J+k)],1,0,0, 
                                                 as.numeric(table.count[2,1,1]))
            model.data.unsorted[(counter-2),]<-c(names(data)[i], names(data)[(I+j)],names(data)[(I+J+k)],1,0,1, 
                                                 as.numeric(table.count[2,1,2]))
            model.data.unsorted[(counter-1),]<-c(names(data)[i], names(data)[(I+j)],names(data)[(I+J+k)],1,1,0, 
                                                 as.numeric(table.count[2,2,1]))
            model.data.unsorted[(counter),]  <-c(names(data)[i], names(data)[(I+j)],names(data)[(I+J+k)],1,1,1, 
                                                 as.numeric(table.count[2,2,2]))
          }
        }
      }
      colnames(model.data.unsorted)<-c("W", "Y", "Z", "wi", "yj", "zk", "count")       
    }
    model.data.unsorted[,1:nvars]<-lapply(model.data.unsorted[,1:nvars], as.factor)
    model.data.unsorted[,((nvars+1):(ncol(model.data.unsorted)))]<-(lapply(model.data.unsorted[,((nvars+1):(ncol(model.data.unsorted)))], as.numeric))
    model.data.unsorted[,ncol(model.data.unsorted)]<-(apply(X = as.matrix(model.data.unsorted[,ncol(model.data.unsorted)]), MARGIN = 1, 
                                                            FUN = check.zero, add.constant = add.constant))
  }
  if (!is.null(model.vars)) {
    if (MRCV_globals$print.status) {
      MRCV_globals$pb.counter<-MRCV_globals$pb.counter+MRCV_globals$weight.B.max
      setTxtProgressBar(MRCV_globals$pb, MRCV_globals$pb.counter)
    }
    boot.sample<-cbind(model.vars, data)
    model.data.unsorted<-data.frame(matrix(data = NA, nrow = nrows, ncol = 1))
    counter<-0
    if (nvars==2) {
      for (i in 1:I) {
        for (j in 1:J) {
          counter<-counter+2^nvars
          model.data.unsorted[(counter-3),1]<-sum(boot.sample[((boot.sample[,i]==0)
                                                               &(boot.sample[,(I+j)]==0)),(I+J+1)])
          model.data.unsorted[(counter-2),1]<-sum(boot.sample[((boot.sample[,i]==0)
                                                               &(boot.sample[,(I+j)]==1)),(I+J+1)])
          model.data.unsorted[(counter-1),1]<-sum(boot.sample[((boot.sample[,i]==1)
                                                               &(boot.sample[,(I+j)]==0)),(I+J+1)])
          model.data.unsorted[(counter),1]  <-sum(boot.sample[((boot.sample[,i]==1)
                                                               &(boot.sample[,(I+j)]==1)),(I+J+1)])
        }
      }
    }
    if (nvars==3) {
      for (i in 1:I) {
        for (j in 1:J) {
          for (k in 1:K) {
            counter<-counter+2^nvars
            model.data.unsorted[(counter-7),1]<-sum(boot.sample[((boot.sample[,i]==0)
                                                                 &(boot.sample[,(I+j)]==0)
                                                                 &(boot.sample[,(I+J+k)]==0)),(I+J+K+1)])
            model.data.unsorted[(counter-6),1]<-sum(boot.sample[((boot.sample[,i]==0)
                                                                 &(boot.sample[,(I+j)]==0)
                                                                 &(boot.sample[,(I+J+k)]==1)),(I+J+K+1)])
            model.data.unsorted[(counter-5),1]<-sum(boot.sample[((boot.sample[,i]==0)
                                                                 &(boot.sample[,(I+j)]==1)
                                                                 &(boot.sample[,(I+J+k)]==0)),(I+J+K+1)])
            model.data.unsorted[(counter-4),1]<-sum(boot.sample[((boot.sample[,i]==0)
                                                                 &(boot.sample[,(I+j)]==1)
                                                                 &(boot.sample[,(I+J+k)]==1)),(I+J+K+1)])
            model.data.unsorted[(counter-3),1]<-sum(boot.sample[((boot.sample[,i]==1)
                                                                 &(boot.sample[,(I+j)]==0)
                                                                 &(boot.sample[,(I+J+k)]==0)),(I+J+K+1)])
            model.data.unsorted[(counter-2),1]<-sum(boot.sample[((boot.sample[,i]==1)
                                                                 &(boot.sample[,(I+j)]==0)
                                                                 &(boot.sample[,(I+J+k)]==1)),(I+J+K+1)])
            model.data.unsorted[(counter-1),1]<-sum(boot.sample[((boot.sample[,i]==1)
                                                                 &(boot.sample[,(I+j)]==1)
                                                                 &(boot.sample[,(I+J+k)]==0)),(I+J+K+1)])
            model.data.unsorted[(counter),1]  <-sum(boot.sample[((boot.sample[,i]==1)
                                                                 &(boot.sample[,(I+j)]==1)
                                                                 &(boot.sample[,(I+J+k)]==1)),(I+J+K+1)])
          }
        }
      }
    }
    model.data.unsorted[,1]<-apply(X = as.matrix(model.data.unsorted[,1]), MARGIN = 1, 
                                   FUN = check.zero, add.constant = add.constant)
    model.data.unsorted<-data.matrix(model.data.unsorted)
  }
  model.data.unsorted
}

#####################################################################################

# genloglin.fit()
# A function that estimates the model of interest
#limit.output argument needed when function is used with apply
#model.vars argument is needed for estimation involving bootstrap resamples

genloglin.fit<-function(data, model, nvars, limit.output = FALSE, model.vars = NULL) {
  op<-options()
  on.exit(options(op))
  model.data<-data
  if (!is.null(model.vars)){
    if (MRCV_globals$print.status) {
      MRCV_globals$pb.counter<-MRCV_globals$pb.counter+MRCV_globals$B/MRCV_globals$B.use
      setTxtProgressBar(MRCV_globals$pb, MRCV_globals$pb.counter)
    }
    model.data<-data.frame(model.vars, data)
    colnames(model.data)[ncol(model.data)]<-"count"
  }
  if (limit.output & is.null(model.vars)){
    if (MRCV_globals$print.status) {
      MRCV_globals$pb.counter<-MRCV_globals$pb.counter+1
      setTxtProgressBar(MRCV_globals$pb, MRCV_globals$pb.counter)
    }
  }
  options(warn = -1)
  if (model == "spmi") {
    mod.fit<-glm(formula = count ~ -1 + W:Y + wi%in%W:Y + yj%in%W:Y, 
                 data = model.data, family = poisson(link = log))
  }
  if (model == "homogeneous") {
    mod.fit<-glm(formula = count ~ -1 + W:Y + wi%in%W:Y + yj%in%W:Y + wi:yj, 
                 data = model.data, family = poisson(link = log)) 
  }
  if (model == "w.main") {
    mod.fit<-glm(formula = count ~ -1 + W:Y + wi%in%W:Y + yj%in%W:Y + wi:yj 
                 + wi:yj%in%W, data = model.data, family = poisson(link = log)) 
  }
  if (model == "y.main") {
    mod.fit<-glm(formula = count ~ -1 + W:Y + wi%in%W:Y + yj%in%W:Y + wi:yj 
                 + wi:yj%in%Y, data = model.data, family = poisson(link = log)) 
  }
  if (model == "wy.main") {
    mod.fit<-glm(formula = count ~ -1 + W:Y + wi%in%W:Y + yj%in%W:Y + wi:yj 
                 + wi:yj%in%W + wi:yj%in%Y, data = model.data, 
                 family = poisson(link = log)) 
  }
  if (model == "saturated") {
    mod.fit<-glm(formula = count ~ -1 + W:Y + wi%in%W:Y + yj%in%W:Y + wi:yj 
                 + wi:yj%in%W + wi:yj%in%Y + wi:yj%in%W:Y, 
                 data = model.data, family = poisson(link = log))  
  }  
  if (class(model)=="formula") {
    mod.fit<-glm(formula = model, data = model.data, family = poisson(link = log))
    formula.char<-Reduce(paste, deparse(model, width.cutoff=499L))
    mod.fit$call<-paste("glm(formula =", formula.char, ", family = poisson(link = log), data = model.data)", 
                        collapse = " ")
  }
  options(warn = 0)
  output<-mod.fit 
  if (limit.output) {
    output<-c(mod.fit$fitted.values, mod.fit$deviance)
    if (is.null(model.vars)){
      output<-mod.fit$fitted.values
    }
  }
  output
}

#####################################################################################

# ipf.genloglin()
# A function that uses the iterative proportional fitting algorithm to obtain
#   multinomial probabilities under the null hypothesis

ipf.genloglin<-function(data, I, J, K, nvars, p, p.theta.2, p.theta.3, x.theta.2, 
                        x.theta.3) {
  p.theta.2.new<-p.theta.2
  p.theta.3.new<-p.theta.3
  if (nvars==2) {
    for (i in 1:(I+J-1)) {
      for (j in (i+1):(I+J)) {
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j]) 
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==0)),5]<-sum(p[((p[,i]==0)&(p[,j]==0)), 
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==1)),5]<-sum(p[((p[,i]==0)&(p[,j]==1)),
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==0)),5]<-sum(p[((p[,i]==1)&(p[,j]==0)),
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==1)),5]<-sum(p[((p[,i]==1)&(p[,j]==1)),
                                                                     ncol(p)])
        
        p[((p[,i]==0)&(p[,j]==0)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==0)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==0)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==0)&(x.theta.2[,4]==0)),5])
        p[((p[,i]==0)&(p[,j]==1)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==1)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==1)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==0)&(x.theta.2[,4]==1)),5])
        p[((p[,i]==1)&(p[,j]==0)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==0)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==0)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==1)&(x.theta.2[,4]==0)),5])
        p[((p[,i]==1)&(p[,j]==1)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==1)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==1)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==1)&(x.theta.2[,4]==1)),5])
        
        p.theta.2<-p.theta.2.new
      }
    }
  }
  if (nvars==3) {
    for (i in 1:I) {
      for (j in (I+1):(I+J)) {
        for (k in (I+J+1):(I+J+K)) {
          p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==0)&(p.theta.3[,5]==0)
                     &(p.theta.3[,6]==0)),7]<-sum(p[((p[,i]==0)&(p[,j]==0)&(p[,k]==0)),ncol(p)])
          p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==0)&(p.theta.3[,5]==0)
                     &(p.theta.3[,6]==1)),7]<-sum(p[((p[,i]==0)&(p[,j]==0)&(p[,k]==1)),ncol(p)])
          p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==0)&(p.theta.3[,5]==1)
                     &(p.theta.3[,6]==0)),7]<-sum(p[((p[,i]==0)&(p[,j]==1)&(p[,k]==0)),ncol(p)])
          p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==0)&(p.theta.3[,5]==1)
                     &(p.theta.3[,6]==1)),7]<-sum(p[((p[,i]==0)&(p[,j]==1)&(p[,k]==1)),ncol(p)])
          p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==1)&(p.theta.3[,5]==0)
                     &(p.theta.3[,6]==0)),7]<-sum(p[((p[,i]==1)&(p[,j]==0)&(p[,k]==0)),ncol(p)])
          p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==1)&(p.theta.3[,5]==0)
                     &(p.theta.3[,6]==1)),7]<-sum(p[((p[,i]==1)&(p[,j]==0)&(p[,k]==1)),ncol(p)])
          p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==1)&(p.theta.3[,5]==1)
                     &(p.theta.3[,6]==0)),7]<-sum(p[((p[,i]==1)&(p[,j]==1)&(p[,k]==0)),ncol(p)])
          p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k]) &(p.theta.3[,4]==1)&(p.theta.3[,5]==1)
                     &(p.theta.3[,6]==1)),7]<-sum(p[((p[,i]==1)&(p[,j]==1)&(p[,k]==1)),ncol(p)])
          
          p[((p[,i]==0)&(p[,j]==0)&(p[,k]==0)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==0)&(p[,k]==0)),ncol(p)]
                     /p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==0)&(p.theta.3[,5]==0)
                     &(p.theta.3[,6]==0)),7]*x.theta.3[((x.theta.3[,1]==names(data)[i])
                     &(x.theta.3[,2]==names(data)[j])&(x.theta.3[,3]==names(data)[k])
                     &(x.theta.3[,4]==0)&(x.theta.3[,5]==0)&(x.theta.3[,6]==0)),7])
          p[((p[,i]==0)&(p[,j]==0)&(p[,k]==1)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==0)&(p[,k]==1)),ncol(p)]
                     /p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==0)&(p.theta.3[,5]==0)
                     &(p.theta.3[,6]==1)),7]*x.theta.3[((x.theta.3[,1]==names(data)[i])
                     &(x.theta.3[,2]==names(data)[j])&(x.theta.3[,3]==names(data)[k])
                     &(x.theta.3[,4]==0)&(x.theta.3[,5]==0)&(x.theta.3[,6]==1)),7])
          p[((p[,i]==0)&(p[,j]==1)&(p[,k]==0)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==1)&(p[,k]==0)),ncol(p)]
                     /p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==0)&(p.theta.3[,5]==1)
                     &(p.theta.3[,6]==0)),7]*x.theta.3[((x.theta.3[,1]==names(data)[i])
                     &(x.theta.3[,2]==names(data)[j])&(x.theta.3[,3]==names(data)[k])
                     &(x.theta.3[,4]==0)&(x.theta.3[,5]==1)&(x.theta.3[,6]==0)),7])
          p[((p[,i]==0)&(p[,j]==1)&(p[,k]==1)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==1)&(p[,k]==1)),ncol(p)]
                     /p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==0)&(p.theta.3[,5]==1)
                     &(p.theta.3[,6]==1)),7]*x.theta.3[((x.theta.3[,1]==names(data)[i])
                     &(x.theta.3[,2]==names(data)[j])&(x.theta.3[,3]==names(data)[k])
                     &(x.theta.3[,4]==0)&(x.theta.3[,5]==1)&(x.theta.3[,6]==1)),7])
          p[((p[,i]==1)&(p[,j]==0)&(p[,k]==0)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==0)&(p[,k]==0)),ncol(p)]
                     /p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==1)&(p.theta.3[,5]==0)
                     &(p.theta.3[,6]==0)),7]*x.theta.3[((x.theta.3[,1]==names(data)[i])
                     &(x.theta.3[,2]==names(data)[j])&(x.theta.3[,3]==names(data)[k])
                     &(x.theta.3[,4]==1)&(x.theta.3[,5]==0)&(x.theta.3[,6]==0)),7])
          p[((p[,i]==1)&(p[,j]==0)&(p[,k]==1)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==0)&(p[,k]==1)),ncol(p)]
                     /p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==1)&(p.theta.3[,5]==0)
                     &(p.theta.3[,6]==1)),7]*x.theta.3[((x.theta.3[,1]==names(data)[i])
                     &(x.theta.3[,2]==names(data)[j])&(x.theta.3[,3]==names(data)[k])
                     &(x.theta.3[,4]==1)&(x.theta.3[,5]==0)&(x.theta.3[,6]==1)),7])
          p[((p[,i]==1)&(p[,j]==1)&(p[,k]==0)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==1)&(p[,k]==0)),ncol(p)]
                     /p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==1)&(p.theta.3[,5]==1)
                     &(p.theta.3[,6]==0)),7]*x.theta.3[((x.theta.3[,1]==names(data)[i])
                     &(x.theta.3[,2]==names(data)[j])&(x.theta.3[,3]==names(data)[k])
                     &(x.theta.3[,4]==1)&(x.theta.3[,5]==1)&(x.theta.3[,6]==0)),7])
          p[((p[,i]==1)&(p[,j]==1)&(p[,k]==1)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==1)&(p[,k]==1)),ncol(p)]
                     /p.theta.3[((p.theta.3[,1]==names(data)[i])&(p.theta.3[,2]==names(data)[j])
                     &(p.theta.3[,3]==names(data)[k])&(p.theta.3[,4]==1)&(p.theta.3[,5]==1)
                     &(p.theta.3[,6]==1)),7]*x.theta.3[((x.theta.3[,1]==names(data)[i])
                     &(x.theta.3[,2]==names(data)[j])&(x.theta.3[,3]==names(data)[k])
                     &(x.theta.3[,4]==1)&(x.theta.3[,5]==1)&(x.theta.3[,6]==1)),7])
          
          p.theta.3<-p.theta.3.new
        }
      }
    }
    for (i in 1:(I-1)) {
      for (j in (i+1):I) {
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j]) 
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==0)),5]<-sum(p[((p[,i]==0)&(p[,j]==0)), 
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==1)),5]<-sum(p[((p[,i]==0)&(p[,j]==1)),
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==0)),5]<-sum(p[((p[,i]==1)&(p[,j]==0)),
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==1)),5]<-sum(p[((p[,i]==1)&(p[,j]==1)),
                                                                     ncol(p)])
        
        p[((p[,i]==0)&(p[,j]==0)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==0)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==0)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==0)&(x.theta.2[,4]==0)),5])
        p[((p[,i]==0)&(p[,j]==1)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==1)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==1)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==0)&(x.theta.2[,4]==1)),5])
        p[((p[,i]==1)&(p[,j]==0)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==0)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==0)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==1)&(x.theta.2[,4]==0)),5])
        p[((p[,i]==1)&(p[,j]==1)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==1)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==1)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==1)&(x.theta.2[,4]==1)),5])
        
        p.theta.2<-p.theta.2.new
      }
    }
    for (i in (I+1):(I+J-1)) {
      for (j in (i+1):(I+J)) {
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j]) 
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==0)),5]<-sum(p[((p[,i]==0)&(p[,j]==0)), 
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==1)),5]<-sum(p[((p[,i]==0)&(p[,j]==1)),
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==0)),5]<-sum(p[((p[,i]==1)&(p[,j]==0)),
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==1)),5]<-sum(p[((p[,i]==1)&(p[,j]==1)),
                                                                     ncol(p)])
        
        p[((p[,i]==0)&(p[,j]==0)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==0)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==0)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==0)&(x.theta.2[,4]==0)),5])
        p[((p[,i]==0)&(p[,j]==1)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==1)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==1)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==0)&(x.theta.2[,4]==1)),5])
        p[((p[,i]==1)&(p[,j]==0)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==0)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==0)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==1)&(x.theta.2[,4]==0)),5])
        p[((p[,i]==1)&(p[,j]==1)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==1)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==1)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==1)&(x.theta.2[,4]==1)),5])
        
        p.theta.2<-p.theta.2.new
      }
    }
    for (i in (I+J+1):(I+J+K-1)) {
      for (j in (i+1):(I+J+K)) {
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j]) 
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==0)),5]<-sum(p[((p[,i]==0)&(p[,j]==0)), 
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==1)),5]<-sum(p[((p[,i]==0)&(p[,j]==1)),
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==0)),5]<-sum(p[((p[,i]==1)&(p[,j]==0)),
                                                                     ncol(p)])
        p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==1)),5]<-sum(p[((p[,i]==1)&(p[,j]==1)),
                                                                     ncol(p)])
        
        p[((p[,i]==0)&(p[,j]==0)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==0)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==0)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==0)&(x.theta.2[,4]==0)),5])
        p[((p[,i]==0)&(p[,j]==1)),ncol(p)]<-(p[((p[,i]==0)&(p[,j]==1)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==0)&(p.theta.2[,4]==1)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==0)&(x.theta.2[,4]==1)),5])
        p[((p[,i]==1)&(p[,j]==0)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==0)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==0)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==1)&(x.theta.2[,4]==0)),5])
        p[((p[,i]==1)&(p[,j]==1)),ncol(p)]<-(p[((p[,i]==1)&(p[,j]==1)),ncol(p)]
                   /p.theta.2[((p.theta.2[,1]==names(data)[i])&(p.theta.2[,2]==names(data)[j])
                   &(p.theta.2[,3]==1)&(p.theta.2[,4]==1)),5]
                   *x.theta.2[((x.theta.2[,1]==names(data)[i])&(x.theta.2[,2]==names(data)[j])
                   &(x.theta.2[,3]==1)&(x.theta.2[,4]==1)),5])
        
        p.theta.2<-p.theta.2.new
      }
    }
  }
  p
}

#####################################################################################

# check.margins()
# A function that checks whether resamples have all positive or negative responses 
# for an item

check.margins<-function(data, I, J, K, nvars, model.vars, item.names) {
  nrows<-(2^nvars)*I*J*max(1,K)
  model.data<-data.frame(matrix(data = NA, nrow = nrows, ncol = (2*nvars+1)))
  model.data[,1:(2*nvars)]<-model.vars
  model.data[,(2*nvars+1)]<-data
  if (nvars ==2) {
    pos.count<-numeric(I+J)
    neg.count<-numeric(I+J)
    for (i in 1:I) {
      pos.count[i]<-sum(model.data[((model.data[,1]==item.names[i])
                                    &(model.data[,3]==1)),5])
      neg.count[i]<-sum(model.data[((model.data[,1]==item.names[i])
                                    &(model.data[,3]==0)),5])
    }
    for (j in 1:J) {
      pos.count[(I+j)]<-sum(model.data[((model.data[,2]==item.names[(I+j)])
                                        &(model.data[,4]==1)),5])     
      neg.count[(I+j)]<-sum(model.data[((model.data[,2]==item.names[(I+j)])
                                        &(model.data[,4]==0)),5])
    }
  }
  if (nvars ==3) {
    pos.count<-numeric(I+J+K)
    neg.count<-numeric(I+J+K)
    for (i in 1:I) {
      pos.count[i]<-sum(model.data[((model.data[,1]==item.names[i])
                                    &(model.data[,4]==1)),7])
      neg.count[i]<-sum(model.data[((model.data[,1]==item.names[i])
                                    &(model.data[,4]==0)),7])
    }
    for (j in 1:J) {
      pos.count[(I+j)]<-sum(model.data[((model.data[,2]==item.names[(I+j)])
                                        &(model.data[,5]==1)),7])     
      neg.count[(I+j)]<-sum(model.data[((model.data[,2]==item.names[(I+j)])
                                        &(model.data[,5]==0)),7])
    }
    for (k in 1:K) {
      pos.count[(I+J+k)]<-sum(model.data[((model.data[,3]==item.names[(I+J+k)])
                                          &(model.data[,6]==1)),7])     
      neg.count[(I+J+k)]<-sum(model.data[((model.data[,3]==item.names[(I+J+k)])
                                          &(model.data[,6]==0)),7])
    }
  }  
  ((min(pos.count)>0)&(min(neg.count)>0))
}

#####################################################################################

# print.genloglin() 
# A method function used to control display of genloglin() output

print.genloglin<-function(x, digits = max(3, getOption("digits") - 3), ...) {
  op<-options()
  on.exit(options(op))
  options(scipen=5)
  x<-x$mod.fit
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", 
      sep = "")
  if (length(coef(x))) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts))
      cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n\n")
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  cat("Null Deviance:\t   ", format(signif(x$null.deviance, digits)), 
      "\nResidual Deviance:", format(signif(x$deviance, digits)), "\n") 
  invisible(x)
}

#####################################################################################

# genloglin() 
# A function that estimates a model for two MRCVs using a marginal estimation approach

genloglin<-function(data, I, J, K = NULL, model, add.constant = .5, boot = TRUE, 
                    B = 1999, B.max = B, print.status = TRUE) {  
  op<-options()
  on.exit(options(op))
  options(warn=1)
  if((class(data)!="data.frame")&(class(data)!="matrix")) {
    stop("The \"data\" argument requires an object of class \"data.frame\".")
  }
  data<-as.data.frame(data)
  if ((I==0)|(J==0)) {
    stop("\"I\" and \"J\" must be greater than 0.")
  }
  if (!is.numeric(I)|!is.numeric(J)) {
    stop("\"I\" and \"J\" must be numeric.")
  }
  if ((I%%1!=0)|(J%%1!=0)) {
    stop("\"I\" and \"J\" must be integers.")
  }
  I<-as.integer(I)
  J<-as.integer(J)
  if (!is.null(K)) {
    if (!is.numeric(K)) {
      stop("\"K\" must be numeric or NULL.")
    }
    if (K%%1!=0) {
      stop("\"K\" must be an integer.")
    }
    K<-as.integer(K)
    if (K==0) {
      K<-NULL
    }
  }
  if (!is.null(K)) {
    if (class(model)!="formula") {
      stop("For the 3 MRCV case, only user-supplied formulas are accepted by the \"model\" argument.")
    }
  }
  if ((class(model)!="formula")&(model!="spmi")&(model!="homogeneous")&(model!="w.main")&(model!="y.main")&(model!="wy.main")&(model!="saturated")) {
    stop("The \"model\" argument requires a formula or one of 6 recognized character strings. \n  See help(genloglin) for details.")    
  }
  if (!is.numeric(B)|!is.numeric(B.max)) {
    warning("\"B\" and \"B.max\" must be numeric. \n  The input values have been changed to the default value of 1999.")
    B<-1999
    B.max=B
  }
  if (B.max<B) {
    warning("\"B.max\" must be greater than or equal to \"B\". \n  \"B.max\" has been set equal to \"B\".")
    B.max=B
  }
  if ((B%%1!=0)|(B.max%%1!=0)) {
    warning("\"B\" and \"B.max\" must be integers. \n  The input values have been rounded up to the nearest integer.")
  }
  B<-as.integer(ceiling(B))
  B.max<-as.integer(ceiling(B.max))
  if (!is.numeric(add.constant)) {
    if (add.constant==FALSE) {
      add.constant<-0
    }
    if (add.constant!=FALSE) {
      warning("The \"add.constant\" argument only accepts numeric values. \n  The input value has been changed to the default value of .5.")
      add.constant<-.5
    }
  }
  if (add.constant<0) {
    warning("The \"add.constant\" argument cannot be negative. \n  The input value has been changed to the default value of .5.")
    add.constant<-.5
  }
  if ((class(boot)!="logical")&(boot!=1)&(boot!=0)) {
    warning("The \"boot\" argument requires an object of class \"logical\". \n  The input value has been changed to the default value of TRUE.")
    boot<-TRUE
  }
  if ((class(print.status)!="logical")&(print.status!=1)&(print.status!=0)) {
    warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    print.status<-TRUE
  }
  boot<-as.logical(boot)
  print.status<-as.logical(print.status)
  MRCV_globals$B<-B
  MRCV_globals$print.status<-print.status
  n<-nrow(data)
  if (is.numeric(K)) {
    if (K==0) {
      K<-NULL
    }
  }
  nvars<-2+is.numeric(K)
  #Reformat raw data
  model.data.unsorted<-data.format(data = data, I = I, J = J, K = K, nvars = nvars, 
                                   add.constant = add.constant)
  #Need to sort data for later calculations
  if (nvars==2) {
    model.data<-model.data.unsorted[order(-model.data.unsorted$wi,-model.data.unsorted$yj),]
  }
  if (nvars==3) {
    model.data<-model.data.unsorted[order(-model.data.unsorted$wi,-model.data.unsorted$yj,
                                          -model.data.unsorted$zk),]
  }
  #Add additional variables specified in model formula
  if (class(model)=="formula") {
    for (i in 1:I) {
      parm<-paste("W", i, sep = "")
      if (length(agrep(parm, model, max.distance=0)) > 0) {
        model.data<-data.frame(model.data, 
                               as.numeric((model.data[,1]==names(data)[i])))
        colnames(model.data)[ncol(model.data)]<-parm
      } 
    }
    for (j in 1:J) {
      parm<-paste("Y", j, sep = "")
      if (length(agrep(parm, model, max.distance=0)) > 0) {
        model.data<-data.frame(model.data, 
                               as.numeric((model.data[,2]==names(data)[(I+j)])))
        colnames(model.data)[ncol(model.data)]<-parm
      } 
    }
    if (nvars==3) {
      for (k in 1:K) {
        parm<-paste("Z", k, sep = "")
        if (length(agrep(parm, model, max.distance=0)) > 0) {
          model.data<-data.frame(model.data, 
                                 as.numeric((model.data[,3]==names(data)[(I+J+k)])))
          colnames(model.data)[ncol(model.data)]<-parm
        } 
      }
    }
  }
  #Estimate model of interest
  mod.fit<-genloglin.fit(data = model.data, model = model, nvars = nvars)
  model.data<-model.data[,1:(2*nvars+1)]
  mu.hat<-mod.fit$fitted.values
  sum.fit<-summary(mod.fit)
  
  #RS Calculations (Appendix A)
  W.counts<-as.data.frame(table(data[,1:I])) #Get all possible combos of W's
  cols <- c(1:I) #Need to order W's in ascending order starting with W1
  W.counts<-W.counts[do.call("order", as.data.frame(W.counts[,cols])),]
  Y.counts<-as.data.frame(table(data[,(I+1):(I+J)])) #All possible Y's
  cols <- c(1:J) #Need to order Y's in ascending order starting with Y1
  Y.counts<-Y.counts[do.call("order", as.data.frame(Y.counts[,cols])),]
  if (nvars==3) {
    Z.counts<-as.data.frame(table(data[,(I+J+1):(I+J+K)])) #All possible Z's
    cols <- c(1:K) #Need to order Z's in ascending order starting with Z1
    Z.counts<-Z.counts[do.call("order", as.data.frame(Z.counts[,cols])),]
  }
  n.counts<-as.data.frame(table(data)) #Get all possible combos of W's, Y's, and Z's
  cols <- c(1:(ncol(data))) #Need to order W's, Y's, Z's in ascending order W1, W2, ...
  n.counts<-n.counts[do.call("order", as.data.frame(n.counts[,cols])),]
  G<-t(data.matrix(W.counts[,1:I])-1) #rx2^r matrix
  H<-t(data.matrix(Y.counts[,1:J])-1) #cx2^c matrix
  if (nvars==3) {
    L<-t(data.matrix(Z.counts[,1:K])-1) 
  }
  tau<-n.counts[,ncol(n.counts)]/n  #Vector of multinomial probabilities
  Jr<-matrix(data = 1, nrow = I, ncol = 2^I)
  Jc<-matrix(data = 1, nrow = J, ncol = 2^J)
  if (nvars==3) {
    Jq<-matrix(data = 1, nrow = K, ncol = 2^K)
  }
  if (nvars==2) {
    B.matrix<-rbind(kronecker(G,H), kronecker(G,(Jc-H)), kronecker((Jr-G),H), 
                    kronecker((Jr-G),(Jc-H)))                
  }
  if (nvars==3) {
    B.matrix<-rbind(kronecker(kronecker(G,H),L), kronecker(kronecker(G,H),(Jq-L)), 
                    kronecker(kronecker(G,(Jc-H)),L), kronecker(kronecker(G,(Jc-H)),(Jq-L)), 
                    kronecker(kronecker((Jr-G),H),L), kronecker(kronecker((Jr-G),H),(Jq-L)), 
                    kronecker(kronecker((Jr-G),(Jc-H)),L), kronecker(kronecker((Jr-G),(Jc-H)),(Jq-L)))                
  }
  V<-n*B.matrix%*%tcrossprod((diag(tau) - tcrossprod(tau)),B.matrix) #Asymptotic variance for m
  X<-model.matrix(mod.fit)  #Design matrix
  #Covariance matrix for B.hat
  d.mu.hat<-diag(mu.hat)
  XP.mu.X.inv<-solve(crossprod(X,d.mu.hat)%*%X)
  sigma<-tcrossprod(XP.mu.X.inv,X)%*%V%*%X%*%XP.mu.X.inv
  rs.se<-sqrt(diag(sigma)) #RS2 standard errors for B.hats
  cov.mu<-d.mu.hat%*%X%*%tcrossprod(sigma,X)%*%d.mu.hat #Cov matrix for mu.hat
  i.matrix<-diag(nrow(mod.fit$data))
  #Covariance matrix for the residuals (m - mu.hat)
  E.part<-i.matrix-d.mu.hat%*%X%*%tcrossprod(XP.mu.X.inv,X)
  E<-E.part%*%tcrossprod(V,E.part)
  gamma<-Re(eigen(diag(1/mu.hat)%*%E)$values) #Eigenvalues (only use real part)
  #Swap default output with RS2 results
  sum.fit$coefficients[,2]<-rs.se
  sum.fit$coefficients[,3]<-sum.fit$coefficients[,1]/sum.fit$coefficients[,2]
  sum.fit$coefficients[,4]<-2*(1-pnorm(abs(sum.fit$coefficients[,3])))
  rnames<-names(mod.fit$coefficients)
  cnames<-c("Estimate", "RS SE", "z value", "Pr(>|z|)") 
  dimnames(sum.fit$coefficients)<-list(rnames,cnames)
  sum.fit$cov.unscaled<-sigma
  sum.fit$cov.scaled<-sigma
  
  if (boot) {
    #Use the "Gange bootstrap" algorithm
    #Get observed pairwise counts for W's
    counter<-0
    w.m<-data.frame(matrix(data = NA, nrow = 4*choose(I,2), ncol = 5))
    for (i in 1:(I-1)) {
      for (j in (i+1):I) {
        counter<-counter+4
        w.m[(counter-3),]<-c(names(data)[(i)],names(data)[(j)],0,0,
                             table(data[,i],data[,j])[1,1])
        w.m[(counter-2),]<-c(names(data)[(i)],names(data)[(j)],0,1,
                             table(data[,i],data[,j])[1,2])
        w.m[(counter-1),]<-c(names(data)[(i)],names(data)[(j)],1,0,
                             table(data[,i],data[,j])[2,1])
        w.m[(counter),]  <-c(names(data)[(i)],names(data)[(j)],1,1,
                             table(data[,i],data[,j])[2,2])
      }
    }
    #Get observed pairwise counts for Y's
    counter<-0
    y.m<-data.frame(matrix(data = NA, nrow = 4*choose(J,2), ncol = 5))
    for (i in 1:(J-1)) {
      for (j in (i+1):J) {
        counter<-counter+4
        y.m[(counter-3),]<-c(names(data)[(i+I)],names(data)[(j+I)],0,0,
                             table(data[,(i+I)],data[,(j+I)])[1,1])
        y.m[(counter-2),]<-c(names(data)[(i+I)],names(data)[(j+I)],0,1,
                             table(data[,(i+I)],data[,(j+I)])[1,2])
        y.m[(counter-1),]<-c(names(data)[(i+I)],names(data)[(j+I)],1,0,
                             table(data[,(i+I)],data[,(j+I)])[2,1])
        y.m[(counter),]  <-c(names(data)[(i+I)],names(data)[(j+I)],1,1,
                             table(data[,(i+I)],data[,(j+I)])[2,2])
      }
    }
    if (nvars==3) {
      #Get observed pairwise counts for Z's
      counter<-0
      z.m<-data.frame(matrix(data = NA, nrow = 4*choose(K,2), ncol = 5))
      for (i in 1:(K-1)) {
        for (j in (i+1):K) {
          counter<-counter+4
          z.m[(counter-3),]<-c(names(data)[(i+I+J)],names(data)[(j+I+J)],0,0,
                               table(data[,(i+I+J)],data[,(j+I+J)])[1,1])
          z.m[(counter-2),]<-c(names(data)[(i+I+J)],names(data)[(j+I+J)],0,1,
                               table(data[,(i+I+J)],data[,(j+I+J)])[1,2])
          z.m[(counter-1),]<-c(names(data)[(i+I+J)],names(data)[(j+I+J)],1,0,
                               table(data[,(i+I+J)],data[,(j+I+J)])[2,1])
          z.m[(counter),]  <-c(names(data)[(i+I+J)],names(data)[(j+I+J)],1,1,
                               table(data[,(i+I+J)],data[,(j+I+J)])[2,2])
        }
      }
    }
    #Get model estimated counts for each pair
    nrows<-(2^nvars)*I*J*max(1,K)
    est.m<-data.frame(matrix(data = NA, nrow = nrows, ncol = (2*nvars+1)))
    est.m[,1:(ncol(est.m)-1)]<-model.data[,1:(ncol(est.m)-1)]
    est.m[,ncol(est.m)]<-mu.hat
    #Initialize multinomial probability matrix
    p<-n.counts
    p[,ncol(p)]<-.5
    colnames(p)<-c(names(data), "p")
    if (nvars==2) {
      x.theta.2<-rbind(w.m,y.m,est.m)
    }
    if (nvars==3) {
      x.theta.2<-rbind(w.m,y.m,z.m)
    }
    x.theta.2[,1:2]<-lapply(x.theta.2[,1:2], as.factor)
    x.theta.2[,3:5]<-lapply(x.theta.2[,3:5], as.numeric)
    x.theta.2[,5]<-apply(X = as.matrix(x.theta.2[,5]), MARGIN = 1, FUN = check.zero, 
                         add.constant = add.constant)
    x.theta.2[,5]<-x.theta.2[,5]/n
    p.theta.2<-x.theta.2
    p.theta.2[,5]<-0
    p.theta.3<-NULL
    x.theta.3<-NULL
    if (nvars==3) {
      est.m[,7]<-apply(X = as.matrix(est.m[,7]), MARGIN = 1, FUN = check.zero, 
                       add.constant = add.constant)
      est.m[,7]<-est.m[,7]/n
      p.theta.3<-est.m
      p.theta.3[,7]<-0
      x.theta.3<-est.m
    }
    tol<-0.00000001
    save.p<-1
    counter<-1
    #Use the iterative proportional fitting algorithm
    while(max(abs(save.p-p[,ncol(p)]))>tol){ 
      save.p<-p[,ncol(p)]
      p<-ipf.genloglin(data = data, I = I, J = J, K = K, nvars = nvars, p = p, 
                       p.theta.2 = p.theta.2, p.theta.3 = p.theta.3, x.theta.2 = x.theta.2, 
                       x.theta.3 = x.theta.3)
      if (print.status) {
        if (((counter%%5)==0)|(counter==1)) {
          cat(counter, "iterations of the iterative proportional fitting algorithm completed", "\n")
        }
      }
      counter<-counter+1
    }  
    if (print.status) {
      cat(counter, "total iterations required", "\n")
    }
    #Generate bootstrap resamples under the null hypothesis
    boot.sample<-cbind(p[,1:(ncol(p)-1)], rmultinom(n = B.max, size = n, 
                                                    prob = p[,ncol(p)]))
    model.data.star.unsorted<-data.frame(matrix(data = NA, nrow = nrows, 
                                                ncol = (2*nvars+B.max)))
    model.data.star.unsorted[,1:(2*nvars)]<-model.data.unsorted[,1:(2*nvars)]
    #Create progress bar for bootstrapping
    if (print.status) {
      cat("Bootstrap Progress:", "\n")
      weight.B.max<-((5+2/3)*B)/B.max
      MRCV_globals$weight.B.max<-weight.B.max
      MRCV_globals$pb<-txtProgressBar(min = 0, max = (B.max*weight.B.max+B), style = 3)
      MRCV_globals$pb.counter<-0
    }
    #Reformat the bootstrap resamples into model form
    model.data.star.unsorted[,(2*nvars+1):(2*nvars+B.max)]<-apply(X = 
          as.matrix(boot.sample[,((ncol(p)):(ncol(boot.sample)))]), MARGIN = 2, 
          FUN = data.format, I = I, J = J, K = K, nvars = nvars, 
          add.constant = add.constant, model.vars = boot.sample[,(1:(ncol(p)-1))])
    if (nvars==2) {
      colnames(model.data.star.unsorted)<-c("W", "Y", "wi", "yj", rep("count", B.max))  
      model.data.star<-model.data.star.unsorted[order(-model.data.star.unsorted$wi, 
                                                      -model.data.star.unsorted$yj),]
    }
    if (nvars==3) {
      colnames(model.data.star.unsorted)<-c("W", "Y", "Z", "wi", "yj", "zk", 
                                            rep("count", B.max))  
      model.data.star<-model.data.star.unsorted[order(-model.data.star.unsorted$wi, 
                                                      -model.data.star.unsorted$yj,
                                                      -model.data.star.unsorted$zk),]
    }
    #Only keep valid resamples (no items with all positive or negative counts)
    keep<-apply(X = as.matrix(model.data.star[,(2*nvars+1):(2*nvars+B.max)]), MARGIN = 2, 
                FUN = check.margins, I=I, J=J, K=K, nvars = nvars, 
                model.vars=model.data.star[,1:(2*nvars)], item.names=names(data))
    model.data.star<-model.data.star[,c(rep(TRUE, 2*nvars),keep)]
    B.discard<-B.max-ncol(model.data.star)+2*nvars
    #Only keep B resamples (or all valid resamples if value < B)
    model.data.star<-model.data.star[,1:(min((B+2*nvars),(ncol(model.data.star)+2*nvars)))]
    B.use<-ncol(model.data.star)-2*nvars
    MRCV_globals$B.use<-B.use
    #Add additional variables specified in model formula
    if (class(model)=="formula") {
      for (i in 1:I) {
        parm<-paste("W", i, sep = "")
        if (length(agrep(parm, model, max.distance=0)) > 0) {
          model.data.star<-data.frame(model.data.star, 
                                      as.numeric((model.data.star[,1]==names(data)[i])))
          colnames(model.data.star)[ncol(model.data.star)]<-parm
        } 
      }
      for (j in 1:J) {
        parm<-paste("Y", j, sep = "")
        if (length(agrep(parm, model, max.distance=0)) > 0) {
          model.data.star<-data.frame(model.data.star, 
                                      as.numeric((model.data.star[,2]==names(data)[(I+j)])))
          colnames(model.data.star)[ncol(model.data.star)]<-parm
        } 
      }
      if (nvars==3) {
        for (k in 1:K) {
          parm<-paste("Z", k, sep = "")
          if (length(agrep(parm, model, max.distance=0)) > 0) {
            model.data.star<-data.frame(model.data.star, 
                                        as.numeric((model.data.star[,3]==names(data)[(I+J+k)])))
            colnames(model.data.star)[ncol(model.data.star)]<-parm
          } 
        }
      }
    }
    var.cols<-c(1:(2*nvars))
    if (ncol(model.data.star) > (B.use+2*nvars)) {
      var.cols<-c(1:(2*nvars),(B.use+2*nvars+1):ncol(model.data.star))
    }
    #Estimate null model for each resample and keep mu.hat and deviance output
    mod.fit.star<-t(apply(X = as.matrix(model.data.star[,(2*nvars+1):(B.use+2*nvars)]), 
                          MARGIN = 2, FUN = genloglin.fit, model = model, nvars = nvars,
                          limit.output = TRUE, model.vars = model.data.star[,var.cols]))
    if (print.status) {
      setTxtProgressBar(MRCV_globals$pb, (B.max*weight.B.max+B))
      close(MRCV_globals$pb)
    }
    model.data.star<-model.data.star[,1:(B.use+2*nvars)]
    mu.hat.star<-mod.fit.star[,(1:nrows)]
    deviance.star<-mod.fit.star[,ncol(mod.fit.star)]
    #Assume alternative model is saturated (user can specify a different alternative
    #   model in anova.genloglin())
    mu.hat.star.HA<-t(model.data.star[,(2*nvars+1):(B.use+2*nvars)])
    #Get t*'s
    chisq.star<-rowSums(((mu.hat.star.HA-mu.hat.star)^2)/mu.hat.star)
    lrt.star<-deviance.star
    residual.star<-mu.hat.star.HA-mu.hat.star
    
    boot.results<-list(B.use = B.use, B.discard = B.discard, model.data.star = 
                         model.data.star, mod.fit.star = mod.fit.star, chisq.star = 
                         chisq.star, lrt.star = lrt.star, residual.star = 
                         residual.star)
  }
  mod.fit<-mod.fit[-11]
  sum.fit<-sum.fit[-5]
  original.arg<-list(data = data, I = I, J = J, K = K, nvars = nvars, model = model, 
                     add.constant = add.constant, boot = boot)
  rs.results<-list(cov.mu = cov.mu, E = E, gamma = gamma)
  if (!boot) {
    output<-list(original.arg = original.arg, mod.fit = mod.fit, sum.fit = sum.fit, 
                 rs.results = rs.results)
  }
  if (boot) {
    output<-list(original.arg = original.arg, mod.fit = mod.fit, sum.fit = sum.fit, 
                 rs.results = rs.results, boot.results = boot.results) 
  }
  class(output)<-"genloglin"
  output
}

#####################################################################################

# print.summary.genloglin
# A method function that controls the display of summary.genloglin()

print.summary.genloglin<-function(x, digits = max(3, getOption("digits") - 3), 
                                  symbolic.cor = x$symbolic.cor, signif.stars = 
                                  getOption("show.signif.stars"), ...) {
  op<-options()
  on.exit(options(op))
  options(scipen=5)
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", 
      sep = "")
  cat("Deviance Residuals: \n")
  if (x$df.residual > 5) {
    x$deviance.resid<-quantile(x$deviance.resid, na.rm = TRUE)
    names(x$deviance.resid)<-c("Min", "1Q", "Median", "3Q", "Max")
  }
  xx<-zapsmall(x$deviance.resid, digits + 1)
  print.default(xx, digits = digits, na.print = "", print.gap = 2)
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    df<-if ("df" %in% names(x)) 
      x[["df"]]
    else NULL
    if (!is.null(df) && (nsingular <- df[3L] - df[1L])) 
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
  }
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
      format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null", 
      "    Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance", 
      "deviance")]), digits = max(5, digits + 1))), 1L, paste, collapse = " "), sep = "")
  if (nzchar(mess <- naprint(x$na.action))) 
    cat("  (", mess, ")\n", sep = "")
  cat("\n","        Number of Fisher Scoring iterations: ", x$iter, "\n", sep = "")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2, digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}

#####################################################################################

# summary.genloglin() 
# A method function that summarizes results given by genloglin()

summary.genloglin<-function(object, ...) {
  op<-options()
  on.exit(options(op))
  options(scipen=5)
  sum.fit<-object$sum.fit
  class(sum.fit)<-"summary.genloglin"
  sum.fit    
}

#####################################################################################

# residuals.genloglin() 
# A method function that calculates standardized residuals

residuals.genloglin<-function(object, ...) {
  data<-object$original.arg$data
  boot<-object$original.arg$boot
  model<-object$original.arg$model
  nvars<-object$original.arg$nvars
  I<-object$original.arg$I
  J<-object$original.arg$J
  K<-object$original.arg$K
  model.data<-object$mod.fit$data[,1:(2*nvars+1)]
  mu.hat<-object$mod.fit$fitted.values
  mu.hat.sat<-model.data[,(2*nvars+1)]
  E<-object$rs.results$E
  
  #Calculate standardized pearson residuals using the estimated asymptotic variance
  resid.num<-mu.hat.sat-mu.hat
  for (i in 1:length(resid.num)) {
    if (abs(resid.num[i]) < .000000001) {
      resid.num[i]<-0
    }
  }
  std.pearson.res.asymp.var<-resid.num/sqrt(diag(E))
  std.pearson.res.asymp.var<-data.frame(model.data[,1:(2*nvars)], 
                                        res = round(std.pearson.res.asymp.var,2))
  if (nvars==2) {
    res.table.asymp<-tabular(Heading()*res*Heading()*(mean)*W*Heading()
                             *Factor(wi, wi, c(1,0))~Heading()*Y*Heading()
                             *Factor(yj,yj,c(1,0)), data = std.pearson.res.asymp.var)
    output<-list(std.pearson.res.asymp.var = res.table.asymp)
  }
  if (nvars==3) {
    output<-list(std.pearson.res.asymp.var = std.pearson.res.asymp.var)
  }
  
  if (boot) {
    B.use<-object$boot.results$B.use
    B.discard<-object$boot.results$B.discard
    residual.star<-object$boot.results$residual.star
    resid.denom<-apply(X = residual.star, MARGIN = 2, FUN = sd)
    #Calculate standardized pearson residuals using the bootstrap variance
    std.pearson.res.boot.var<-resid.num/resid.denom
    std.pearson.res.boot.var<-data.frame(model.data[,1:(2*nvars)], 
                                         res = round(std.pearson.res.boot.var,2))
    if (nvars==2) {
      res.table.boot<-tabular(Heading()*res*Heading()*(mean)*W*Heading()
                              *Factor(wi, wi, c(1,0))~Heading()*Y*Heading()
                              *Factor(yj,yj,c(1,0)), data = std.pearson.res.boot.var)
      output<-list(std.pearson.res.asymp.var = res.table.asymp, B.use = B.use, 
                   B.discard = B.discard, std.pearson.res.boot.var = res.table.boot)
    }
    if (nvars==3) {
      output<-list(std.pearson.res.asymp.var = std.pearson.res.asymp.var, 
                   B.use = B.use, B.discard = B.discard, std.pearson.res.boot.var = 
                   std.pearson.res.boot.var)
    }
  }
  output
}

#####################################################################################

# print.anova.genloglin()
# A function used to control the display of output provided by anova.genloglin()

print.anova.genloglin<-function(x, ...) {
  op<-options()
  on.exit(options(op))
  options(scipen=10)
  type<-names(x)
  model<-x$original.arg$model
  if (class(model)=="formula") {
    model<-Reduce(paste, deparse(model, width.cutoff=499L))
  }
  model.HA<-x$original.arg$model.HA
  if (class(model.HA)=="formula") {
    model.HA<-Reduce(paste, deparse(model.HA, width.cutoff=499L))
  }
  gof<-x$original.arg$gof
  Pearson.chisq<-round(x$test.statistics$Pearson.chisq, 2)
  lrt<-round(x$test.statistics$lrt, 2)
  if (model=="saturated") {
    chisq.p.rs<-x$rs.results$Pearson.chisq.rs$p.value
    lrt.p.rs<-x$rs.results$lrt.rs$p.value
    chisq.p.boot<-x$boot.results$Pearson.chisq.boot$p.value
    lrt.p.boot<-x$boot.results$lrt.boot$p.value
    cat("\n")
    cat("Model comparison statistics for", "\n")
    cat("H0 =", model, "\n")
    cat("HA =", model.HA, "\n", "\n")
    cat("Pearson chi-square statistic =", Pearson.chisq, "\n")
    cat("LRT statistic =", lrt, "\n", "\n")
    cat("Second-Order Rao-Scott Adjusted Results:", "\n")
    cat("Pearson chi-square p-value =", chisq.p.rs, "\n")
    cat("LRT p-value =", lrt.p.rs, "\n", "\n")
    cat("Bootstrap Results:", "\n")
    cat("Pearson chi-square p-value =", chisq.p.boot, "\n")
    cat("LRT p-value =", lrt.p.boot, "\n", "\n")
  }
  if (model!="saturated") {
    cat("\n")
    cat("Model comparison statistics for", "\n") 
    cat("H0 =", model, "\n")
    cat("HA =", model.HA, "\n", "\n")
    cat("Pearson chi-square statistic =", Pearson.chisq, "\n")
    cat("LRT statistic =", lrt, "\n")
    if (any(type=="rs.results")){
      Pearson.chisq.rs<-round(x$rs.results$Pearson.chisq.rs$Pearson.chisq.rs, 2)
      df.rs<-round(x$rs.results$Pearson.chisq.rs$df, 2)
      chisq.p.rs<-as.name(paste("=", round(x$rs.results$Pearson.chisq.rs$p.value, 4)))
      if (round(x$rs.results$Pearson.chisq.rs$p.value, 4)<.0001) {
        chisq.p.rs<-as.name(paste("<", .0001))
      }
      lrt.rs<-round(x$rs.results$lrt.rs$lrt.rs, 2)
      lrt.p.rs<-as.name(paste("=", round(x$rs.results$lrt.rs$p.value, 4)))
      if (round(x$rs.results$lrt.rs$p.value, 4)<.0001) {
        lrt.p.rs<-as.name(paste("<", .0001))
      }
      cat("\n")
      cat("Second-Order Rao-Scott Adjusted Results:", "\n")
      cat("Rao-Scott Pearson chi-square statistic =", paste(Pearson.chisq.rs, ",", sep = ""), "df =", paste(df.rs, ",", sep = ""), "p", chisq.p.rs, "\n")
      cat("Rao-Scott LRT statistic =", paste(lrt.rs, ",", sep = ""), "df =", paste(df.rs, ",", sep = ""), "p", lrt.p.rs, "\n")
      if (any(type=="boot.results")){
        B.use<-x$boot.results$B.use
        B.discard<-x$boot.results$B.discard
        p.chisq.boot<-as.name(paste("=", round(x$boot.results$p.chisq.boot, 4)))
        if (x$boot.results$p.chisq.boot==0) {
          p.chisq.boot<-as.name(paste("<", round(1/B.use, 4)))
        }
        p.lrt.boot<-as.name(paste("=", round(x$boot.results$p.lrt.boot, 4)))
        if (x$boot.results$p.lrt.boot==0) {
          p.lrt.boot<-as.name(paste("<", round(1/B.use, 4)))
        }
        cat("\n")
        cat("Bootstrap Results:", "\n")
        if (B.discard > 0) {
          cat(B.discard, "resamples were removed from the analysis due to")  
          cat(" not having all rows or columns represented in a 2x2 table", "\n")
        }
        cat("Final results based on", B.use, "resamples", "\n")
        cat("Pearson chi-square p-value", p.chisq.boot, "\n")
        cat("LRT p-value", p.lrt.boot, "\n")
      }
    }
    if (gof&(model.HA!="saturated")) {
      cat("\n")
      cat("-------------------------------------------------------------------------------------")
      cat("\n", "\n")
      Pearson.chisq.gof = round(x$test.statistics$Pearson.chisq.gof, 2)
      lrt.gof = round(x$test.statistics$lrt.gof, 2)
      cat("Goodness of fit statistics for", "\n")
      cat("H0 =", model, "\n", "\n")
      cat("Pearson chi-square GOF statistic =", Pearson.chisq.gof, "\n")
      cat("LRT GOF statistic =", lrt.gof, "\n")
      if (any(type=="rs.results")){
        Pearson.chisq.gof.rs<-(round(x$rs.results$Pearson.chisq.gof.rs$Pearson.chisq.gof.rs, 2))
        df.rs<-round(x$rs.results$Pearson.chisq.gof.rs$df, 2)
        chisq.gof.p.rs<-(as.name(paste("=", round(x$rs.results$Pearson.chisq.gof.rs$p.value, 4))))
        if (round(x$rs.results$Pearson.chisq.gof.rs$p.value, 4)<.0001) {
          chisq.gof.p.rs<-as.name(paste("<", .0001))
        }
        lrt.gof.rs<-round(x$rs.results$lrt.gof.rs$lrt.gof.rs, 2)
        lrt.gof.p.rs<-as.name(paste("=", round(x$rs.results$lrt.gof.rs$p.value, 4)))
        if (round(x$rs.results$lrt.gof.rs$p.value, 4)<.0001) {
          lrt.gof.p.rs<-as.name(paste("<", .0001))
        }
        cat("\n")
        cat("Second-Order Rao-Scott Adjusted Results:", "\n")
        cat("Rao-Scott Pearson chi-square GOF statistic =", paste(Pearson.chisq.gof.rs, ",", sep = ""), "df =", paste(df.rs, ",", sep = ""), "p", chisq.gof.p.rs, "\n")
        cat("Rao-Scott LRT GOF statistic =", paste(lrt.gof.rs, ",", sep = ""), "df =", paste(df.rs, ",", sep = ""), "p", lrt.gof.p.rs, "\n")
      }
      if (any(type=="boot.results")){
        p.chisq.gof.boot<-as.name(paste("=", round(x$boot.results$p.chisq.gof.boot, 4)))
        if (x$boot.results$p.chisq.gof.boot==0) {
          p.chisq.gof.boot<-as.name(paste("<", round(1/B.use, 4)))
        }
        p.lrt.gof.boot<-as.name(paste("=", round(x$boot.results$p.lrt.gof.boot, 4)))
        if (x$boot.results$p.lrt.gof.boot==0) {
          p.lrt.gof.boot<-as.name(paste("<", round(1/B.use, 4)))
        }
        cat("\n")
        cat("Bootstrap Results:", "\n")
        cat("Pearson chi-square GOF p-value", p.chisq.gof.boot, "\n")
        cat("LRT GOF p-value", p.lrt.gof.boot, "\n")
      }
    }
    cat("\n")
  }
  invisible(x)
}

#####################################################################################

# anova.genloglin() 
# A method function for comparing two models through Pearson and lrt statistics

anova.genloglin<-function(object, model.HA="saturated", type = "all", gof = TRUE, 
                          print.status = TRUE, ...) {
  op<-options()
  on.exit(options(op))
  options(warn=1)
  if ((class(print.status)!="logical")&(print.status!=1)&(print.status!=0)) {
    warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    print.status<-TRUE
  }
  if ((class(gof)!="logical")&(gof!=1)&(gof!=0)) {
    warning("The \"gof\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    gof<-TRUE
  }
  print.status<-as.logical(print.status)
  gof<-as.logical(gof)
  if (length(type)>1) {
    warning("The \"type\" argument requires an object of length 1. \n  Only the first element is used.")
    type<-type[1]
  }
  if (!(type%in%c("all","boot","rs2"))) {
    warning("The \"type\" argument can only take on values of \"boot\", \"rs2\", and \"all\". \n  The input value has been changed to the default value of \"all\".")
    type<-"all"
  }
  boot<-object$original.arg$boot
  if(!boot) {
    if (type!="rs2") {
      warning("You must specify the boot option in genloglin() in order to obtain bootstrap results. \n  The \"type\" argument has been changed to \"rs2\".")
      type = "rs2"
    }
  }
  data<-object$original.arg$data
  nvars<-object$original.arg$nvars
  I<-object$original.arg$I
  J<-object$original.arg$J
  K<-object$original.arg$K
  if (!is.null(K)) {
    if ((class(model.HA)!="formula")&(model.HA!="saturated")) {
      warning("For the 3 MRCV case, only \"saturated\" or user-supplied formulas are accepted by the \"model.HA\" argument. \n  The input value has been changed to the default value of \"saturated\".")
      model.HA<-"saturated"
    }
  }
  if ((class(model.HA)!="formula")&(model.HA!="homogeneous")&(model.HA!="w.main")&(model.HA!="y.main")&(model.HA!="wy.main")&(model.HA!="saturated")) {
    warning("The \"model.HA\" argument requires a formula or one of 5 recognized character strings. \n  See help(anova.genloglin) for details. \n  The input value has been changed to the default value of \"saturated\".")    
    model.HA<-"saturated"
  }
  model<-object$original.arg$model
  model.data<-object$mod.fit$data[,1:(2*nvars+1)]
  mu.hat<-object$mod.fit$fitted.values
  mu.hat.HA<-model.data[,(2*nvars+1)]
  deviance<-object$mod.fit$deviance
  gamma<-object$rs.results$gamma
  MRCV_globals$print.status<-print.status
  #Compute observed test statistics
  chisq.obs<-sum(((mu.hat.HA-mu.hat)^2)/mu.hat)
  lrt.obs<-deviance
  chisq.gof.obs<-chisq.obs
  lrt.gof.obs<-lrt.obs
  #If alternative model is not saturated model then need to estimate model
  if (model.HA!="saturated") {    
    if (class(model.HA)=="formula") {
      for (i in 1:I) {
        parm<-paste("W", i, sep = "")
        if (length(agrep(parm, model.HA, max.distance=0)) > 0) {
          model.data<-data.frame(model.data, 
                                 as.numeric((model.data[,1]==names(data)[i])))
          colnames(model.data)[ncol(model.data)]<-parm
        } 
      }
      for (j in 1:J) {
        parm<-paste("Y", j, sep = "")
        if (length(agrep(parm, model.HA, max.distance=0)) > 0) {
          model.data<-data.frame(model.data, 
                                 as.numeric((model.data[,2]==names(data)[(I+j)])))
          colnames(model.data)[ncol(model.data)]<-parm
        } 
      }
      if (nvars==3) {
        for (k in 1:K) {
          parm<-paste("Z", k, sep = "")
          if (length(agrep(parm, model.HA, max.distance=0)) > 0) {
            model.data<-data.frame(model.data, 
                                   as.numeric((model.data[,3]==names(data)[(I+J+k)])))
            colnames(model.data)[ncol(model.data)]<-parm
          } 
        }
      }
    }
    mod.fit.HA<-genloglin.fit(data = model.data, model = model.HA, nvars = nvars)
    model.data<-model.data[,1:(2*nvars+1)]
    mu.hat.HA<-mod.fit.HA$fitted.values
    deviance.HA<-mod.fit.HA$deviance
    chisq.obs<-sum(((mu.hat.HA-mu.hat)^2)/mu.hat)
    lrt.obs<-deviance-deviance.HA
  }
  original.arg<-list(model = model, model.HA = model.HA, gof = gof)
  test.statistics<-list(Pearson.chisq = chisq.obs, lrt = lrt.obs)
  #If gof=TRUE and alternative model is not saturated then provide gof statistics
  if (gof&(model.HA!="saturated")) {
    test.statistics<-list(Pearson.chisq = chisq.obs, lrt = lrt.obs, 
                          Pearson.chisq.gof = chisq.gof.obs, 
                          lrt.gof = lrt.gof.obs)
  }
  
  if (any(type == "rs2" | type == "all")) {
    sum.gamma<-sum(gamma)
    sum.gamma.sq<-sum(gamma^2)
    chisq.rs<-(sum.gamma*chisq.obs)/sum.gamma.sq
    lrt.rs<-(sum.gamma*lrt.obs)/sum.gamma.sq
    chisq.gof.rs<-(sum.gamma*chisq.gof.obs)/sum.gamma.sq
    lrt.gof.rs<-(sum.gamma*lrt.gof.obs)/sum.gamma.sq
    df.rs<-(sum.gamma^2)/sum.gamma.sq
    p.chisq.rs<-1-pchisq(q = chisq.rs, df = df.rs)
    p.lrt.rs<-1-pchisq(q = lrt.rs, df = df.rs)
    p.chisq.gof.rs<-1-pchisq(q = chisq.gof.rs, df = df.rs)
    p.lrt.gof.rs<-1-pchisq(q = lrt.gof.rs, df = df.rs)
    rs.results<-list(Pearson.chisq.rs = list(Pearson.chisq.rs = chisq.rs, 
                     df = df.rs, p.value = p.chisq.rs), lrt.rs = 
                     list(lrt.rs = lrt.rs, df = df.rs, p.value = p.lrt.rs))
    if (gof&(model.HA!="saturated")) {
      rs.results<-list(Pearson.chisq.rs = list(Pearson.chisq.rs = chisq.rs, 
                     df = df.rs, p.value = p.chisq.rs), lrt.rs = list(lrt.rs = 
                     lrt.rs, df = df.rs, p.value = p.lrt.rs), 
                     Pearson.chisq.gof.rs = list(Pearson.chisq.gof.rs = 
                     chisq.gof.rs, df = df.rs, p.value = p.chisq.gof.rs), 
                     lrt.gof.rs = list(lrt.gof.rs = lrt.gof.rs, df = df.rs, 
                     p.value = p.lrt.gof.rs))
    }
    output<-list(original.arg = original.arg, test.statistics = test.statistics, 
                 rs.results = rs.results)
  }
  
  if (any(type == "boot" | type == "all")) {
    nrows<-(2^nvars)*I*J*max(1,K)
    B.use<-object$boot.results$B.use
    B.discard<-object$boot.results$B.discard
    model.data.star<-object$boot.results$model.data.star
    mu.hat.star<-object$boot.results$mod.fit.star[,1:nrows]
    deviance.star<-object$boot.results$mod.fit.star[,(nrows+1)]
    chisq.star<-object$boot.results$chisq.star
    lrt.star<-object$boot.results$lrt.star
    chisq.gof.star<-chisq.star
    lrt.gof.star<-lrt.star
    #If alternative model is not saturated then need to estimate t*'s for HA
    if (model.HA!="saturated") {
      if (class(model.HA)=="formula") {
        for (i in 1:I) {
          parm<-paste("W", i, sep = "")
          if (length(agrep(parm, model.HA, max.distance=0)) > 0) {
            model.data.star<-data.frame(model.data.star, 
                                        as.numeric((model.data.star[,1]==names(data)[i])))
            colnames(model.data.star)[ncol(model.data.star)]<-parm
          } 
        }
        for (j in 1:J) {
          parm<-paste("Y", j, sep = "")
          if (length(agrep(parm, model.HA, max.distance=0)) > 0) {
            model.data.star<-data.frame(model.data.star, 
                                        as.numeric((model.data.star[,2]==names(data)[(I+j)])))
            colnames(model.data.star)[ncol(model.data.star)]<-parm
          } 
        }
        if (nvars==3) {
          for (k in 1:K) {
            parm<-paste("Z", k, sep = "")
            if (length(agrep(parm, model.HA, max.distance=0)) > 0) {
              model.data.star<-data.frame(model.data.star, 
                                          as.numeric((model.data.star[,3]==names(data)[(I+J+k)])))
              colnames(model.data.star)[ncol(model.data.star)]<-parm
            } 
          }
        }
      }
      var.cols<-c(1:(2*nvars))
      if (ncol(model.data.star) > (B.use+(2*nvars))) {
        var.cols<-c(1:(2*nvars),(B.use+(2*nvars+1)):ncol(model.data.star))
      }
      #Create progress bar for bootstrapping
      if (print.status) {
        cat("Bootstrap Progress:", "\n")
        MRCV_globals$pb <- txtProgressBar(min = 0, max = B.use, style = 3)
        MRCV_globals$pb.counter<-0
      }
      mod.fit.star.HA<-matrix(NA, nrow = B.use, ncol = (nrows+1))
      mod.fit.star.HA<-t(apply(X = 
                               as.matrix(model.data.star[,(2*nvars+1):(B.use+(2*nvars))]), 
                               MARGIN = 2, FUN = genloglin.fit, model = model.HA, 
                               nvars = nvars, limit.output = TRUE, model.vars = 
                               model.data.star[,var.cols]))
      if (print.status) {
        setTxtProgressBar(MRCV_globals$pb, B.use)
        close(MRCV_globals$pb)
      }
      model.data.star<-model.data.star[,1:(B.use+2*nvars)]
      mu.hat.star.HA<-mod.fit.star.HA[,(1:nrows)]
      deviance.star.HA<-mod.fit.star.HA[,ncol(mod.fit.star.HA)]
      chisq.star<-rowSums(((mu.hat.star.HA-mu.hat.star)^2)/mu.hat.star)
      lrt.star<-deviance.star-deviance.star.HA
    }
    p.chisq.boot<-(1/B.use)*sum(chisq.star>=chisq.obs)
    p.lrt.boot<-(1/B.use)*sum(lrt.star>=lrt.obs)
    p.chisq.gof.boot<-(1/B.use)*sum(chisq.gof.star>=chisq.gof.obs)
    p.lrt.gof.boot<-(1/B.use)*sum(lrt.gof.star>=lrt.gof.obs)
    boot.results<-list(B.use = B.use, B.discard = B.discard, p.chisq.boot = 
                       p.chisq.boot, p.lrt.boot = p.lrt.boot)
    if (gof&(model.HA!="saturated")) {
      boot.results<-list(B.use = B.use, B.discard = B.discard, p.chisq.boot = 
                         p.chisq.boot, p.lrt.boot = p.lrt.boot, p.chisq.gof.boot = 
                         p.chisq.gof.boot, p.lrt.gof.boot = p.lrt.gof.boot)
    } 
    output<-list(original.arg = original.arg, test.statistics = test.statistics, 
                 boot.results = boot.results)
  }
  if (type == "all") {
    output<-list(original.arg = original.arg, test.statistics = test.statistics, 
                 rs.results = rs.results, boot.results = boot.results)
  }
  if (model=="saturated") {
    test.statistics<-list(Pearson.chisq = 0, lrt = 0)
    rs.results<-list(Pearson.chisq.rs = list(p.value = 1), 
                     lrt.rs = list(p.value = 1))
    boot.results<-list(Pearson.chisq.boot = list(p.value = 1), lrt.boot = 
                         list(p.value = 1))
    output<-list(original.arg = original.arg, test.statistics = test.statistics, 
                 rs.results = rs.results, boot.results = boot.results)
  }
  class(output)<-"anova.genloglin"
  output
}

#####################################################################################

# est.jack()
# A function that calculates model estimated odds ratios based on n-1 observations

est.jack<-function(mu.hat, i, I, J, K) {
  nrows<-I*J*max(1,K)
  output<-exp(log(mu.hat[i]) - log(mu.hat[(nrows+i)]) - log(mu.hat[(2*nrows+i)]) 
              + log(mu.hat[(3*nrows+i)]))
}

#####################################################################################

# print.predict.genloglin()
# A function used to control the display of output provided by predict.genloglin()

print.predict.genloglin<-function(x, ...) {
  op<-options()
  on.exit(options(op))
  options(scipen=5)
  type<-names(x)
  data<-x$original.arg$data
  nvars<-x$original.arg$nvars
  I<-x$original.arg$I
  J<-x$original.arg$J
  K<-x$original.arg$K
  alpha<-x$original.arg$alpha
  coverage<-paste(((1-alpha)*100), "%", sep = "")
  obs<-round(x$OR.obs,2)
  model.asymp<-round(x$OR.model.asymp,2)
  if (nvars==2) {
    OR.obs<-matrix(data = NA, nrow = I, ncol = J)
    OR.model.asymp<-matrix(data = NA, nrow = I, ncol = J)
    for (i in 1:I) {
      for (j in 1:J) {
        OR.obs[i,j]<-paste(obs[((i-1)*J+j),1], paste("(", paste(obs[((i-1)*J+j),2], 
                                                                obs[((i-1)*J+j),3], sep = ", "),")", sep = ""))
        OR.model.asymp[i,j]<-paste(model.asymp[((i-1)*J+j),1], 
                                   paste("(", paste(model.asymp[((i-1)*J+j),2], 
                                         model.asymp[((i-1)*J+j),3], sep = ", "),")", 
                                         sep = ""))
      }
    }
    rownames(OR.obs)<-names(data)[1:I]
    colnames(OR.obs)<-names(data)[(I+1):(I+J)]
    rownames(OR.model.asymp)<-names(data)[1:I]
    colnames(OR.model.asymp)<-names(data)[(I+1):(I+J)]
    cat("Observed odds ratios with", coverage, "asymptotic confidence intervals", "\n") 
    print.default(OR.obs, quote = FALSE)
    cat("\n")
    cat("Model-predicted odds ratios with", coverage, "asymptotic confidence intervals", "\n") 
    print.default(OR.model.asymp, quote = FALSE)
    cat("\n")
  }
  if (nvars==3) {
    cat("Observed odds ratios with", coverage, "asymptotic confidence intervals", "\n") 
    print.default(obs, quote = FALSE)
    cat("\n")
    cat("Model-predicted odds ratios with", coverage, "asymptotic confidence intervals", "\n") 
    print.default(model.asymp, quote = FALSE)
    cat("\n") 
  }
  
  if (any(type=="boot.results")){
    B.use<-x$boot.results$B.use
    B.discard<-x$boot.results$B.discard
    model.BCa<-round(x$boot.results$OR.model.BCa, 2)
    if (nvars==2) {
      OR.model.BCa<-matrix(data = NA, nrow = I, ncol = J)
      for (i in 1:I) {
        for (j in 1:J) {
          OR.model.BCa[i,j]<-paste(model.BCa[((i-1)*J+j),1], 
                                   paste("(", paste(model.BCa[((i-1)*J+j),2], 
                                         model.BCa[((i-1)*J+j),3], sep = ", "),")", 
                                         sep = ""))
        }
      }
      rownames(OR.model.BCa)<-names(data)[1:I]
      colnames(OR.model.BCa)<-names(data)[(I+1):(I+J)]
      cat("Bootstrap Results:", "\n")
      if (B.discard > 0) {
        cat(B.discard, "resamples were removed from the analysis due to")  
        cat(" not having all rows or columns represented in a 2x2 table", "\n")
      }
      cat("Final results based on", B.use, "resamples", "\n")
      cat("Model-predicted odds ratios with", coverage, "bootstrap BCa confidence intervals", "\n") 
      print.default(OR.model.BCa, quote = FALSE)
      cat("\n")
    }
    if (nvars==3) {
      cat("Bootstrap Results:", "\n")
      if (B.discard > 0) {
        cat(B.discard, "resamples were removed from the analysis due to")  
        cat(" not having all rows or columns represented in a 2x2 table", "\n")
      }
      cat("Final results based on", B.use, "resamples", "\n")
      cat("Model-predicted odds ratios with", coverage, "bootstrap BCa confidence intervals", "\n") 
      print.default(model.BCa, quote = FALSE)
      cat("\n")
    }
  }
  invisible(x)
}

#####################################################################################

# predict.genloglin() 
# A function that calculates observed and model-estimated odds ratios and their 
#   corresponding confidence intervals

predict.genloglin<-function(object, alpha = .05, pair = "WY", print.status = TRUE, ...) {
  op<-options()
  on.exit(options(op))
  options(warn=1)
  if (!is.numeric(alpha)) {
    warning("The \"alpha\" argument only accepts numeric values. \n  The input value has been changed to the default value of .05.")
    alpha<-.05
  }
  if ((alpha<=0)|(alpha>=1)) {
    warning("The \"alpha\" argument must be between 0 and 1. \n  The input value has been changed to the default value of .05.")
    alpha<-.05
  }
  if (length(pair)>1) {
    warning("The \"pair\" argument requires an object of length 1. \n  Only the first element is used.")
    pair<-pair[1]
  }
  if (!(pair%in%c("WY","WZ","YZ","YW","ZW","ZY","wy","wz","yz","yw","zw","zy"))) {
    warning("The \"pair\" argument can only take on values of \"WY\", \"WZ\", and \"YZ\". \n  The input value has been changed to the default value of \"WY\".")
    pair<-"WY"
  }
  if ((class(print.status)!="logical")&(print.status!=1)&(print.status!=0)) {
    warning("The \"print.status\" argument requires a logical value. \n  The input value has been changed to the default value of TRUE.")
    print.status<-TRUE
  }
  print.status<-as.logical(print.status)
  boot<-object$original.arg$boot
  data<-object$original.arg$data
  n<-nrow(data)
  nvars<-object$original.arg$nvars
  I<-object$original.arg$I
  J<-object$original.arg$J
  K<-object$original.arg$K
  add.constant<-object$original.arg$add.constant
  model<-object$original.arg$model
  model.data<-object$mod.fit$data[,1:(2*nvars+1)]
  mu.hat<-object$mod.fit$fitted.values
  cov.mu<-object$rs.results$cov.mu
  nrows<-I*J*max(1,K)
  MRCV_globals$print.status<-print.status
  if (nvars==3) {
    if (any(pair == "WZ" | pair == "ZW"| pair == "wz"| pair == "zw")) {
      cov.mu<-cbind(cov.mu,model.data[,1:6])
      cov.mu<-cov.mu[order(-cov.mu$yj,-cov.mu$wi,-cov.mu$zk,cov.mu$Y,cov.mu$W,cov.mu$Z),]
      cov.mu<-cbind(t(cov.mu[,1:(ncol(cov.mu)-6)]),model.data[,1:6])
      cov.mu<-cov.mu[order(-cov.mu$yj,-cov.mu$wi,-cov.mu$zk,cov.mu$Y,cov.mu$W,cov.mu$Z),]
      cov.mu<-t(cov.mu[,1:(ncol(cov.mu)-6)])
      mu.hat<-cbind(as.data.frame(mu.hat),model.data[,1:6])
      mu.hat<-mu.hat[order(-mu.hat$yj,-mu.hat$wi,-mu.hat$zk,mu.hat$Y,mu.hat$W,mu.hat$Z),]
      mu.hat<-as.matrix(mu.hat[,1])
      model.data<-model.data[order(-model.data$yj,-model.data$wi,-model.data$zk,
                                   model.data$Y,model.data$W,model.data$Z),]
    }
    if (any(pair == "WY" | pair == "YW"| pair == "wy"| pair == "yw")) {
      cov.mu<-cbind(cov.mu,model.data[,1:6])
      cov.mu<-cov.mu[order(-cov.mu$zk,-cov.mu$wi,-cov.mu$yj,cov.mu$Z,cov.mu$W,cov.mu$Y),]
      cov.mu<-cbind(t(cov.mu[,1:(ncol(cov.mu)-6)]),model.data[,1:6])
      cov.mu<-cov.mu[order(-cov.mu$zk,-cov.mu$wi,-cov.mu$yj,cov.mu$Z,cov.mu$W,cov.mu$Y),]
      cov.mu<-t(cov.mu[,1:(ncol(cov.mu)-6)])
      mu.hat<-cbind(as.data.frame(mu.hat),model.data[,1:6])
      mu.hat<-mu.hat[order(-mu.hat$zk,-mu.hat$wi,-mu.hat$yj,mu.hat$Z,mu.hat$W,mu.hat$Y),]
      mu.hat<-as.matrix(mu.hat[,1])
      model.data<-model.data[order(-model.data$zk,-model.data$wi,-model.data$yj,
                                   model.data$Z,model.data$W,model.data$Y),]
    }
  }
  #Get q (used for computing asymptotic variance of model estimated ORs)
  if (nvars==2) {
    q<-matrix(data = NA, nrow = nrows, ncol = ((2^nvars)*nrows))
    for (i in 1:nrows) {
      delta<-as.matrix(c(rep(0,(i-1)),(1/mu.hat[i]),rep(0,(nrows-1)),
                         (-1/mu.hat[(nrows+i)]),rep(0,(nrows-1)),
                         (-1/mu.hat[(2*nrows+i)]),rep(0,(nrows-1)),
                         (1/mu.hat[(3*nrows+i)]),rep(0,(nrows-i))))
      q[i,]<-delta
    }
    logOR.obs<-matrix(data = NA, nrow = nrows, ncol = 3)
    var.logOR.obs<-matrix(data = NA, nrow = nrows, ncol = 1)
    logOR.est<-matrix(data = NA, nrow = nrows, ncol = 3)
  }
  if (nvars==3) {
    q<-matrix(data = NA, nrow = 2*nrows, ncol = ((2^nvars)*nrows))
    for (i in 1:nrows) {
      delta<-as.matrix(c(rep(0,(i-1)),(1/mu.hat[i]),rep(0,(nrows-1)),
                         (-1/mu.hat[(nrows+i)]),rep(0,(nrows-1)),
                         (-1/mu.hat[(2*nrows+i)]),rep(0,(nrows-1)),
                         (1/mu.hat[(3*nrows+i)]),rep(0,(nrows-i)),
                         rep(0,(4*nrows))))
      q[i,]<-delta
    }
    for (i in 1:nrows) {
      delta<-as.matrix(c(rep(0,(4*nrows)),rep(0,(i-1)),(1/mu.hat[(i+4*nrows)]),rep(0,(nrows-1)),
                         (-1/mu.hat[(i+5*nrows)]),rep(0,(nrows-1)),
                         (-1/mu.hat[(i+6*nrows)]),rep(0,(nrows-1)),
                         (1/mu.hat[(i+7*nrows)]),rep(0,(nrows-i))))
      q[(nrows+i),]<-delta
    }
    logOR.obs<-matrix(data = NA, nrow = 2*nrows, ncol = 3)
    var.logOR.obs<-matrix(data = NA, nrow = 2*nrows, ncol = 1)
    logOR.est<-matrix(data = NA, nrow = 2*nrows, ncol = 3)
  }
  #Asymptotic covariance matrix for model estimated ORs
  var.logOR.asymp<-q%*%tcrossprod(cov.mu,q)
  var.logOR.asymp<-as.matrix(diag(var.logOR.asymp))
  for(i in 1:nrows) {
    #Observed ORs and corresponding CIs
    logOR.obs[i,1]<-(log(model.data[i,(2*nvars+1)]) 
                     -log(model.data[(nrows+i),(2*nvars+1)]) 
                     -log(model.data[(2*nrows+i),(2*nvars+1)])
                     +log(model.data[(3*nrows+i),(2*nvars+1)]))
    var.logOR.obs[i,1]<-(1/model.data[i,(2*nvars+1)] 
                         + 1/model.data[(nrows+i),(2*nvars+1)] 
                         + 1/model.data[(2*nrows+i),(2*nvars+1)] 
                         + 1/model.data[(3*nrows+i),(2*nvars+1)])
    logOR.obs[i,2]<-(logOR.obs[i,1] - qnorm(1-alpha/2)*sqrt(var.logOR.obs[i,1]))
    logOR.obs[i,3]<-(logOR.obs[i,1] + qnorm(1-alpha/2)*sqrt(var.logOR.obs[i,1]))
    #Model estimated ORs and corresponding CIs based on asymptotic variance
    logOR.est[i,1]<-(log(mu.hat[i]) - log(mu.hat[(nrows+i)])
                     - log(mu.hat[(2*nrows+i)]) + log(mu.hat[(3*nrows+i)]))
    #For SPMI, prevent warnings/NAs by setting lower and upper CI bounds to 1
    if (model=="spmi") {
      logOR.est[i,2]<-log(1)
      logOR.est[i,3]<-log(1)
    }
    if (model!="spmi") {
      logOR.est[i,2]<-(logOR.est[i,1] - qnorm(1-alpha/2)*sqrt(var.logOR.asymp[i,1]))
      logOR.est[i,3]<-(logOR.est[i,1] + qnorm(1-alpha/2)*sqrt(var.logOR.asymp[i,1]))
    }
  }
  if (nvars==3) {
    for(i in (4*nrows+1):(5*nrows)) {
      #Observed ORs and corresponding CIs
      logOR.obs[(i-3*nrows),1]<-(log(model.data[i,(2*nvars+1)]) 
                                 -log(model.data[(nrows+i),(2*nvars+1)]) 
                                 -log(model.data[(2*nrows+i),(2*nvars+1)])
                                 +log(model.data[(3*nrows+i),(2*nvars+1)]))
      var.logOR.obs[(i-3*nrows),1]<-(1/model.data[i,(2*nvars+1)] 
                                     + 1/model.data[(nrows+i),(2*nvars+1)] 
                                     + 1/model.data[(2*nrows+i),(2*nvars+1)] 
                                     + 1/model.data[(3*nrows+i),(2*nvars+1)])
      logOR.obs[(i-3*nrows),2]<-(logOR.obs[(i-3*nrows),1] 
                                 - qnorm(1-alpha/2)*sqrt(var.logOR.obs[(i-3*nrows),1]))
      logOR.obs[(i-3*nrows),3]<-(logOR.obs[(i-3*nrows),1] 
                                 + qnorm(1-alpha/2)*sqrt(var.logOR.obs[(i-3*nrows),1]))
      #Model estimated ORs and corresponding CIs based on asymptotic variance
      logOR.est[(i-3*nrows),1]<-((log(mu.hat[i]) - log(mu.hat[(nrows+i)]) 
                                  - log(mu.hat[(2*nrows+i)]) + log(mu.hat[(3*nrows+i)])))
      logOR.est[(i-3*nrows),2]<-(logOR.est[(i-3*nrows),1] 
                                 - qnorm(1-alpha/2)*sqrt(var.logOR.asymp[(i-3*nrows),1]))
      logOR.est[(i-3*nrows),3]<-(logOR.est[(i-3*nrows),1] 
                                 + qnorm(1-alpha/2)*sqrt(var.logOR.asymp[(i-3*nrows),1]))
    }
    for (i in 1:(2*nrows)) {
      if (abs(exp(logOR.est[i,1])-1)<.00000001) {
        logOR.est[i,2]<-log(1)
        logOR.est[i,3]<-log(1)
      }
    }
  }
  OR.obs<-exp(logOR.obs) 
  OR.model.asymp<-exp(logOR.est)
  counter<-0
  if (nvars==2) {
    rowlabels<-matrix(data=NA, nrow=1, ncol=nrows)
    for (i in 1:I){
      for (j in (I+1):(I+J)) {
        counter<-counter+1
        cell<-paste(colnames(data)[i], colnames(data)[j], sep = "")
        rowlabels[,counter]<-cell
      }
    }
  }
  if (nvars==3) {
    rowlabels<-matrix(data=NA, nrow=1, ncol=2*nrows)
    if (any(pair == "YZ" | pair == "ZY"| pair == "yz"| pair == "zy")) {
      for (h in 1:0) {
        for (i in 1:I){
          for (j in (I+1):(I+J)) {
            for (k in (I+J+1):(I+J+K)) {
              counter<-counter+1
              cell<-paste(colnames(data)[i], "=", h, ",", colnames(data)[j], colnames(data)[k], sep = "")
              rowlabels[,counter]<-cell
            }
          }
        }
      }
    }
    if (any(pair == "WZ" | pair == "ZW"| pair == "wz"| pair == "zw")) {
      for (h in 1:0) {
        for (j in (I+1):(I+J)){
          for (i in 1:I) {
            for (k in (I+J+1):(I+J+K)) {
              counter<-counter+1
              cell<-paste(colnames(data)[j], "=", h, ",", colnames(data)[i], colnames(data)[k], sep = "")
              rowlabels[,counter]<-cell
            }
          }
        }
      }
    }
    if (any(pair == "WY" | pair == "YW"| pair == "wy"| pair == "yw")) {
      for (h in 1:0) {
        for (k in (I+J+1):(I+J+K)){
          for (i in 1:I) {
            for (j in (I+1):(I+J)) {
              counter<-counter+1
              cell<-paste(colnames(data)[k], "=", h, ",", colnames(data)[i], colnames(data)[j], sep = "")
              rowlabels[,counter]<-cell
            }
          }
        }
      }
    }
  }
  colnames(OR.obs)<-c("OR", "lower.bound", "upper.bound")
  rownames(OR.obs)<-rowlabels
  colnames(OR.model.asymp)<-c("OR", "lower.bound", "upper.bound")
  rownames(OR.model.asymp)<-rowlabels
  original.arg<-list(data = data, I = I, J = J, K = K, nvars = nvars, alpha = alpha)
  output<-list(original.arg = original.arg, OR.obs = OR.obs, 
               OR.model.asymp = OR.model.asymp)
  
  if (boot) {
    B.use<-object$boot.results$B.use
    B.discard<-object$boot.results$B.discard
    if (model=="spmi") {
      bca.ci.lower<-matrix(data = 1, nrow = 1, ncol = nrows)
      bca.ci.upper<-matrix(data = 1, nrow = 1, ncol = nrows)
    }
    if (model!="spmi") {
      mu.hat.star<-object$boot.results$mod.fit.star[,(1:((2^nvars)*nrows))]
      #Perform jackknife calculations
      data.n_1<-merge(data, 1:n, by=NULL, all=TRUE)
      data.n_1<-data.n_1[-seq(from = 1, to = n*n, by = (n+1)),]
      #Create progress bar for bootstrapping
      if (print.status) {
        cat("Bootstrap Progress:", "\n")
        MRCV_globals$pb <- txtProgressBar(min = 0, max = (1.7*n+n), style = 3)
        MRCV_globals$pb.counter<-0
      }
      if (nvars==2) {
        model.data.unsorted.n_1<-by(data = data.n_1[,1:(I+J)], 
                                    INDICES = data.n_1[,(I+J+1)], 
                                    FUN = data.format, I = I, J = J, K = K, 
                                    nvars = nvars, add.constant = 
                                    add.constant, predict.func = TRUE)
      }
      if (nvars==3) {
        model.data.unsorted.n_1<-by(data = data.n_1[,1:(I+J+K)], 
                                    INDICES = data.n_1[,(I+J+K+1)], 
                                    FUN = data.format, I = I, J = J, K = K, 
                                    nvars = nvars, add.constant = 
                                    add.constant, predict.func = TRUE)
      }
      model.data.unsorted.n_1<-as.data.frame(do.call(rbind, model.data.unsorted.n_1))
      model.data.unsorted.n_1<-data.frame(model.data.unsorted.n_1, index = rep(c(1:n),
                               each=((2^nvars)*nrows)))      
      if (class(model)=="formula") {
        for (i in 1:I) {
          parm<-paste("W", i, sep = "")
          if (length(agrep(parm, model, max.distance=0)) > 0) {
            model.data.unsorted.n_1<-data.frame(model.data.unsorted.n_1, 
                                                as.numeric((model.data.unsorted.n_1[,1]==names(data)[i])))
            colnames(model.data.unsorted.n_1)[ncol(model.data.unsorted.n_1)]<-parm
          } 
        }
        for (j in 1:J) {
          parm<-paste("Y", j, sep = "")
          if (length(agrep(parm, model, max.distance=0)) > 0) {
            model.data.unsorted.n_1<-data.frame(model.data.unsorted.n_1, 
                                                as.numeric((model.data.unsorted.n_1[,2]==names(data)[(I+j)])))
            colnames(model.data.unsorted.n_1)[ncol(model.data.unsorted.n_1)]<-parm
          } 
        }
        if (nvars==3) {
          for (k in 1:K) {
            parm<-paste("Z", k, sep = "")
            if (length(agrep(parm, model, max.distance=0)) > 0) {
              model.data.unsorted.n_1<-data.frame(model.data.unsorted.n_1, 
                                                  as.numeric((model.data.unsorted.n_1[,3]==names(data)[(I+J+k)])))
              colnames(model.data.unsorted.n_1)[ncol(model.data.unsorted.n_1)]<-parm
            } 
          }
        }
      }
      var.cols<-c(1:(2*nvars+1))
      if (ncol(model.data.unsorted.n_1) > (2*nvars+2)) {
        var.cols<-c(1:(2*nvars+1),(2*nvars+3):ncol(model.data.unsorted.n_1))
      }
      mu.hat.n_1<-by(data = model.data.unsorted.n_1[,var.cols], 
                     INDICES = model.data.unsorted.n_1[,(2*nvars+2)], 
                     FUN = genloglin.fit, model = model, nvars = nvars, 
                     limit.output = TRUE)
      if (print.status) {
        setTxtProgressBar(MRCV_globals$pb, (1.7*n+n))
        close(MRCV_globals$pb)
      }
      model.data.unsorted.n_1<-model.data.unsorted.n_1[,1:(2*nvars+2)]
      mu.hat.n_1<-t(as.data.frame(matrix(do.call(rbind, mu.hat.n_1), nrow = n, 
                                         ncol = ((2^nvars)*nrows))))
      mu.hat.n_1<-data.frame(model.data.unsorted.n_1[(1:((2^nvars)*nrows)),1:(2*nvars)], 
                             mu.hat.n_1)
      if (nvars==2) {
        mu.hat.n_1<-mu.hat.n_1[order(-mu.hat.n_1$wi, -mu.hat.n_1$yj), ]
      }
      if (nvars==3) {
        mu.hat.n_1<-mu.hat.n_1[order(-mu.hat.n_1$wi, -mu.hat.n_1$yj, -mu.hat.n_1$zk,
                                     mu.hat.n_1$W, mu.hat.n_1$Y, mu.hat.n_1$Z), ]
        if (any(pair == "WZ" | pair == "ZW"| pair == "wz"| pair == "zw")) {
          mu.hat.n_1<-mu.hat.n_1[order(-mu.hat.n_1$yj, -mu.hat.n_1$wi, -mu.hat.n_1$zk,
                                       mu.hat.n_1$Y, mu.hat.n_1$W, mu.hat.n_1$Z), ]
          model.data<-model.data[order(-model.data$wi,-model.data$yj,-model.data$zk,
                                       model.data$W, model.data$Y, model.data$Z),]
          mu.hat.star<-cbind(t(mu.hat.star),model.data[,1:6])
          mu.hat.star<-mu.hat.star[order(-mu.hat.star$yj,-mu.hat.star$wi,-mu.hat.star$zk,
                                         mu.hat.star$Y, mu.hat.star$W, mu.hat.star$Z),]
          mu.hat.star<-t(mu.hat.star[,1:(ncol(mu.hat.star)-6)])
        }
        if (any(pair == "WY" | pair == "YW"| pair == "wy"| pair == "yw")) {
          mu.hat.n_1<-mu.hat.n_1[order(-mu.hat.n_1$zk, -mu.hat.n_1$wi, -mu.hat.n_1$yj,
                                       mu.hat.n_1$Z, mu.hat.n_1$W, mu.hat.n_1$Y), ]
          model.data<-model.data[order(-model.data$wi,-model.data$yj,-model.data$zk,
                                       model.data$W, model.data$Y, model.data$Z),]
          mu.hat.star<-cbind(t(mu.hat.star),model.data[,1:6])
          mu.hat.star<-mu.hat.star[order(-mu.hat.star$zk,-mu.hat.star$wi,-mu.hat.star$yj,
                                         mu.hat.star$Z, mu.hat.star$W, mu.hat.star$Y),]
          mu.hat.star<-t(mu.hat.star[,1:(ncol(mu.hat.star)-6)])
        }
      }
      mu.hat.n_1<-mu.hat.n_1[,-c(1:(2*nvars))]
      #Calculate ORs for bootstrap resamples and jackknife samples
      if (nvars==2) {
        or.star<-matrix(data = NA, nrow = B.use, ncol = nrows)
        or.est<-matrix(data = NA, nrow = 1, ncol = nrows)
        or.est.jack<-matrix(data = NA, nrow = n, ncol = nrows)
      }
      if (nvars==3) {
        or.star<-matrix(data = NA, nrow = B.use, ncol = 2*nrows)
        or.est<-matrix(data = NA, nrow = 1, ncol = 2*nrows)
        or.est.jack<-matrix(data = NA, nrow = n, ncol = 2*nrows)
      }
      for (i in 1:nrows) {
        or.star[,i]<-exp(log(mu.hat.star[,i]) - log(mu.hat.star[,(nrows+i)]) 
                         - log(mu.hat.star[, (2*nrows+i)]) 
                         + log(mu.hat.star[,(3*nrows+i)]))
        or.star[,i]<-or.star[do.call("order", as.data.frame(or.star[,i])),i]
        or.est.jack[,i]<-apply(X = mu.hat.n_1, MARGIN = 2, FUN = est.jack, i = i, 
                               I = I, J = J, K = K)
      }
      if (nvars==3) {
        for(i in (4*nrows+1):(5*nrows)) {
          or.star[,(i-3*nrows)]<-exp(log(mu.hat.star[,i]) - log(mu.hat.star[,(nrows+i)]) 
                                     - log(mu.hat.star[, (2*nrows+i)]) 
                                     + log(mu.hat.star[,(3*nrows+i)]))
          or.star[,(i-3*nrows)]<-or.star[do.call("order", as.data.frame(or.star[,(i-3*nrows)])),(i-3*nrows)]
          or.est.jack[,(i-3*nrows)]<-apply(X = mu.hat.n_1, MARGIN = 2, 
                                           FUN = est.jack, i = i, 
                                           I = I, J = J, K = K)
        }
      }
      #Make matrices conformable for subsequent calculations
      if (nvars==2) {
        or.est.rep.n<-merge(t(exp(logOR.est)[,1]), 1:n, by=NULL, 
                            all=TRUE)[,1:nrows]
        or.est.rep.B<-merge(t(exp(logOR.est)[,1]), 1:B.use, by=NULL, 
                            all=TRUE)[,1:nrows]
      }
      if (nvars==3) {
        or.est.rep.n<-merge(t(exp(logOR.est)[,1]), 1:n, by=NULL, 
                            all=TRUE)[,1:(2*nrows)]
        or.est.rep.B<-merge(t(exp(logOR.est)[,1]), 1:B.use, by=NULL, 
                            all=TRUE)[,1:(2*nrows)]
      }
      #Compute BCa confidence intervals
      l.jack<-(n-1)*(or.est.rep.n-or.est.jack)
      a.jack<-1/6*colSums(l.jack^3)/colSums(l.jack^2)^(3/2)
      w.p<-colSums(or.star<=or.est.rep.B)/(B.use+1)  
      w<-qnorm(p = w.p)
      alpha.bca<-c(alpha/2, 1-alpha/2)
      z.tilde.1<-w + qnorm(p = alpha.bca[1])
      z.tilde.2<-w + qnorm(p = alpha.bca[2])
      alpha.tilde.1<-pnorm(q = w + z.tilde.1/(1-a.jack*z.tilde.1))
      alpha.tilde.2<-pnorm(q = w + z.tilde.2/(1-a.jack*z.tilde.2))
      if (nvars==2) {
        bca.ci.lower<-matrix(data = NA, nrow = 1, ncol = nrows)
        bca.ci.upper<-matrix(data = NA, nrow = 1, ncol = nrows)
      }
      if (nvars==3) {
        bca.ci.lower<-matrix(data = NA, nrow = 1, ncol = 2*nrows)
        bca.ci.upper<-matrix(data = NA, nrow = 1, ncol = 2*nrows)
      }
      for (i in 1:nrows) {
        bca.ci.lower[1,i]<-quantile(or.star[,i], alpha.tilde.1[i])
        bca.ci.upper[1,i]<-quantile(or.star[,i], alpha.tilde.2[i])
      }
      if (nvars==3) {
        for (i in (nrows+1):(2*nrows)) {
          bca.ci.lower[1,i]<-quantile(or.star[,i], alpha.tilde.1[i])
          bca.ci.upper[1,i]<-quantile(or.star[,i], alpha.tilde.2[i])
        }
        for (i in 1:(2*nrows)) {
          if (abs(exp(logOR.est[i,1])-1)<.00000001) {
            bca.ci.lower[1,i]<-1
            bca.ci.upper[1,i]<-1
          }
        }
      }
    }
    OR.model.BCa<-as.matrix(data.frame(exp(logOR.est)[,1], 
                                       t(bca.ci.lower), t(bca.ci.upper)))
    colnames(OR.model.BCa)<-c("OR", "lower.bound", "upper.bound")
    rownames(OR.model.BCa)<-rowlabels
    boot.results<-list(B.use = B.use, B.discard = B.discard, OR.model.BCa = 
                         OR.model.BCa)
    output<-list(original.arg = original.arg, OR.obs = OR.obs, OR.model.asymp = 
                   OR.model.asymp, boot.results = boot.results)
  }
  class(output)<-"predict.genloglin"
  output
}

#####################################################################################