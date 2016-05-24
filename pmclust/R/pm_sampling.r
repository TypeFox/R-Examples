
### Assign sample size to all N.spmd and N.allspmds.
# assign.N.sample <- function(total.sample = 5000, N.org.spmd,
#     var.table.spmd = NULL){
assign.N.sample <- function(total.sample = 5000, N.org.spmd){
  N.new.spmd <- N.org.spmd

  ### Check for SS.
  # if(!is.null(var.table.spmd)){
  #   N.new.spmd <- N.new.spmd - nrow(var.table.spmd)
  #   var.table.spmd <- var.table.spmd[order(var.table.spmd$id),]
  # }
  my.size <- spmd.comm.size()
  N.org.allspmds <- spmd.allgather.integer(as.integer(N.org.spmd),
                                           integer(my.size))
  N.new.allspmds <- spmd.allgather.integer(as.integer(N.new.spmd),
                                           integer(my.size))
  if(any(N.new.allspmds <= 0)){
    stop("N.org.spmd is too small.")
  }

  ### Sampling start.
  if(any(N.new.allspmds < total.sample / my.size)){
    N.sample.spmd <- N.org.spmd
    N.sample.allspmds <- N.org.allspmds
    ID.sample.spmd <- 1:N.org.spmd

    ### Check for SS.
    # if(!is.null(var.table.spmd)){
    #   var.table.spmd$id.sample <-
    #     which(ID.sample.spmd %in% var.table.spmd$id)
    # }
  } else{
    N.sample.allspmds <- floor(N.org.allspmds / sum(N.org.allspmds) *
                                 total.sample)
    remainder <- total.sample - sum(N.sample.allspmds)
    if(remainder > 0){
      N.sample.allspmds[(my.size - remainder + 1):my.size] <-
        N.sample.allspmds[(my.size - remainder + 1):my.size] + 1 
    }

    N.sample.spmd <- N.sample.allspmds[spmd.comm.rank() + 1]

    ### Check for SS.
    # if(!is.null(var.table.spmd)){
    #   ID.sample.spmd <- sample((1:N.org.spmd)[-var.table.spmd$id],
    #                              N.sample.spmd)
    #   ID.sample.spmd <- sort(c(var.table.spmd$id, ID.sample.spmd))

    #   var.table.spmd$id.sample <-
    #     which(ID.sample.spmd %in% var.table.spmd$id)
    # } else{
      ID.sample.spmd <- sort(sample(1:N.org.spmd, N.sample.spmd))
    # }
  }

  N.sample <- sum(N.sample.allspmds)

  list(N = N.sample, N.spmd = N.sample.spmd,
       N.allspmds = N.sample.allspmds,
       ID.spmd = ID.sample.spmd)
#       SS.table.spmd = var.table.spmd)
} # End of assign.N.sample().

