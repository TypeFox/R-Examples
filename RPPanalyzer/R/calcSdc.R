`calcSdc` <-
function(x,
         ## arguments to select and define samples
         sample.id=c("sample","sample.n"),sel=c("measurement","control"),dilution="dilution",
         ## argumets for plot.dilution.series and protein.con
         D0=2,sensible.min=5, sensible.max=1.e9,minimal.err=5, plot=T,
         ## argument only for protein.conc
         r=1.2){

     ## select samples that should be analyzed
     #xi <- select.sample.group(x,param="sample_type",sel=sel)
     xi <- select.sample.group(x,params=list("sample_type"=sel))

     ## add column for sample identification
     xi <- create.ID.col(xi,sample.id)

     ## generate ID vector for measurements
     identifier <- dil.ser.id (x,sample.id)

     ## vector for target identification
     targets <- x[[3]]["target",]

     ## open an empty list to store formated data
     dat.list <- vector("list",length(targets))

     dilutions <- sort(unique(xi[[4]][,dilution]),decreasing=TRUE)

     ## format RPPA list (output sample.median(...)) to make it useable for sdc...
     for (i in seq(along=targets)){

         dat <- x[[1]][,i]       ## calculate Data matrix column wise

         sample.mat <- matrix(0,nrow=length(identifier),ncol=length(dilutions),dimnames=list(identifier,dilutions))

         for (j in seq(along=identifier)){

             sample.lines <- which(xi[[4]][,"identifier"]==identifier[j])

             sample.n <- dat[sample.lines]
             sample.c <- x[[4]][sample.lines,dilution]

             names(sample.n) <- sample.c

             o <- order(names(sample.n),decreasing = T)
             sample.n <- sample.n[o]                     

             sample.mat[identifier[j],] <- sample.n      
         }

         dat.list[[i]] <- sample.mat
     }


     vals.list <- vector("list",length(targets)) 

     for (i in seq(along=dat.list)) {

         vals <-  plotDilutionSeries(D0=D0, data.dilutes=dat.list[[i]]
                 ,sensible.min=sensible.min, sensible.max=sensible.max
                 ,minimal.err=minimal.err
                 , plot=plot, title=targets[i])

         vals.list[[i]] <- vals
     }  


     conc.list <- vector("list",length(targets))

     for (i in seq(along=conc.list)){

         D <- vals.list[[i]]$D
         c <- vals.list[[i]]$c
         a <- vals.list[[i]]$a
         d.D <- vals.list[[i]]$d.D
         d.c <- vals.list[[i]]$d.c
         d.a <- vals.list[[i]]$d.a

         conc <- protein.con(D0=D0,D,c,a,d.D,d.c, d.a,data.dilutes=dat.list[[i]]
                 ,r=1.2,minimal.err=5)

         conc.list[[i]] <- conc
     }

     ## conc.list object in RPPA List format
     concentration <- matrix(0,nrow=nrow(conc.list[[1]]),ncol=length(conc.list),
             dimnames=list(rownames(conc.list[[1]]),names(targets)))
     error <- matrix(0,nrow=nrow(conc.list[[1]]),ncol=length(conc.list),
             dimnames=list(rownames(conc.list[[1]]),names(targets)))

     for (i in seq(along=conc.list)){
         concentration[,i] <- conc.list[[i]][,1]
         error[,i] <- conc.list[[i]][,2]
     }

     id.lines <- match (rownames( conc.list[[1]] ),xi[[4]][,"identifier"])
     sample.description <- xi[[4]][id.lines,]

     data.list <- list(concentration_sdc=concentration
             ,error_sdc=error
             ,arraydescription_sdc=xi[[3]]
             ,sampledescription_sdc=sample.description)


     return(data.list)
}

