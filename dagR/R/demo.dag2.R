demo.dag2 <-
function()
{ # re-creates DAG in figure 2a of Shrier & Platt,
  # BMC Med Res Methodol 2008;8:20;
  dag<-dag.init(y.name="injury", x.name="warm-up exercises",
       covs=c(1,2,1,1,1,1,1,1,1,1,1),
       cov.names=c("coach", "genetics",
                   "fitness", "connective tissue disorder",
                   "pre-game proprioception", "motivation, aggression",
                   "fatigue", "contact sport",
                   "previous injury", "tissue weakness",
                   "intra-game proprioception"));
  dag$x<-c(0.000, 0.261, 0.722, 0.494, 0.995, 0.257,
           0.002, 0.723, 0.505, 0.305, 0.998, 0.502, 1.000);
  dag$y<-c(0.000, 0.852, 0.862, 0.761, 0.735, 0.595,
           0.527, 0.611, 0.449, 0.304, 0.401, 0.149, 0.000);

  dag<-add.arc(dag, c(1,12));
  dag<-add.arc(dag, c(2,7));
  dag<-add.arc(dag, c(2,4));
  dag<-add.arc(dag, c(3,4));
  dag<-add.arc(dag, c(3,8));
  dag<-add.arc(dag, c(3,5));
  dag<-add.arc(dag, c(4,6));
  dag<-add.arc(dag, c(4,8));
  dag<-add.arc(dag, c(5,8));
  dag<-add.arc(dag, c(5,11));
  dag<-add.arc(dag, c(6,1));
  dag<-add.arc(dag, c(7,1));
  dag<-add.arc(dag, c(7,10));
  dag<-add.arc(dag, c(8,12));
  dag<-add.arc(dag, c(8,13));
  dag<-add.arc(dag, c(9,12));
  dag<-add.arc(dag, c(9,10));
  dag<-add.arc(dag, c(11,13));
  dag<-add.arc(dag, c(12,13));

  return(dag);
}

