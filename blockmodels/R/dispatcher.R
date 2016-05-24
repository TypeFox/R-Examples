
dispatcher <- function(membership_name,membership_init,model_name,network,real_EM){
	.Call( "dispatcher",
           membership_name,
           membership_init,
           model_name,
           network,
           real_EM,
          PACKAGE = "blockmodels" )
}

