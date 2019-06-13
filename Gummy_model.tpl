GLOBALS_SECTION
  #include <iostream>
  #include <fstream>
  #undef REPORT
  #define REPORT(object) report<<#object "\n" << object << endl;
  #include <time.h>
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

DATA_SECTION
  //conventional
  init_number nzones;
  init_number nrec_c;
  init_matrix Rec_dat_c(1,nrec_c,1,3);


  //acoustic
  init_number nrec_a;
  init_number nATag;
  init_ivector nATag_obs(1,nATag);
  init_ivector TagID_start(1,nATag);
  init_ivector TagID_end(1,nATag);
  init_3darray Rec_dat_a(1,nATag,TagID_start,TagID_end,1,4); //ragged array of recaptures by tagID, (1st dim), recapture events (2nd dim) and release, recapture info (3rd dim)
  
  int t;
  int y;
  int n;
  //!!cout<<Rec_dat_a<<endl;exit(1);


PARAMETER_SECTION
  init_bounded_number log_p11(-15,-1e-06,1);          //Zone 1
  //init_bounded_number log_p11(-0.2,-1e-06,1);
  init_bounded_number log_p22(-15,-1e-06,1);          //Zone 2
  init_bounded_number log_p21(-15,-1e-06,1);
  //init_bounded_number log_p21(-15,-2,1);
  init_bounded_number log_p33(-15,-1e-06,1);          //Zone 3             
  
  objective_function_value f

  //Declare objects used in Move_rate()
  matrix Mov_mat(1,nzones,1,nzones);
  matrix Move_c(1,nzones,1,nzones);
  matrix Move_a(1,nzones,1,nzones);    
  
  vector Days_libery_c(1,nrec_c);
  vector Rels_c(1,nrec_c);
  vector Recs_c(1,nrec_c);
    
  number p11;
  number p22;
  number p21;
  number p33;
  
  number Tag_NLL_c;
  number Tag_NLL_a;
  number Tag_like;
  number fpen_tag;



PRELIMINARY_CALCS_SECTION
   Days_libery_c=column(Rec_dat_c,1);
   Rels_c=column(Rec_dat_c,2);
   Recs_c=column(Rec_dat_c,3);


PROCEDURE_SECTION
  Move_rate();
  Calc_obj_function();


FUNCTION Move_rate
  Tag_NLL_c.initialize();
  Tag_NLL_a.initialize();  
  fpen_tag.initialize();
  
   // extract movement parameters
  p11=mfexp(log_p11);
  p22=mfexp(log_p22);
  p21=mfexp(log_p21);
  p33=mfexp(log_p33);


  //1. fill in daily Movement matrix
  Mov_mat(1,1)=p11;
  Mov_mat(1,2)=1-p11;
  Mov_mat(1,3)=0;

  Mov_mat(2,1)=p21;
  Mov_mat(2,2)=p22;
  Mov_mat(2,3)=1-(p22+p21);

  Mov_mat(3,1)=0;
  Mov_mat(3,2)=1-p33;
  Mov_mat(3,3)=p33;

  //penalty to keep all row elements positive
  Mov_mat(1,2)=posfun(Mov_mat(1,2),0.000001,fpen_tag);
  fpen_tag+=fpen_tag;
 
  Mov_mat(2,3)=posfun(Mov_mat(2,3),0.000001,fpen_tag);
  fpen_tag+=fpen_tag;
 
  Mov_mat(3,2)=posfun(Mov_mat(3,2),0.000001,fpen_tag);
  fpen_tag+=fpen_tag;

  //penalty to keep all diagonal elements <=1
  Mov_mat(1,1)=1-posfun((1-Mov_mat(1,1)),0.001,fpen_tag);
  fpen_tag+=fpen_tag;
 
  Mov_mat(2,2)=1-posfun((1-Mov_mat(2,2)),0.0001,fpen_tag);
  fpen_tag+=fpen_tag;
 
  Mov_mat(3,3)=1-posfun((1-Mov_mat(3,3)),0.0001,fpen_tag);
  fpen_tag+=fpen_tag;

 
  //2. Predict location of recaptured shark

       //2.1 Conventional tags 
  for(t=1; t<=nrec_c;t++)
  {
     Move_c=Mov_mat;
     if(Days_libery_c(t)>1)
     {    
      for(y=2; y<=Days_libery_c(t);y++) Move_c=Move_c*Mov_mat;      
     }
      int from=value(Rels_c(t));          //C++ syntax: need value() to convet to int
      int to=value(Recs_c(t));
      dvariable Pred_Prob_c=Move_c(from,to);     
      Tag_NLL_c += - log(Pred_Prob_c);  // negative log likelihood of observations     
   }
      
      //2.2 Acoustic tags    
  for(t=1; t<=nATag;t++)                          // loop over each tag ID
  {
     //extract release and recapture info for each TagID
   dvector Days_libery_a=column(Rec_dat_a(t),1);
   dvector Rels_a=column(Rec_dat_a(t),2);
   dvector Recs_a=column(Rec_dat_a(t),3);
    
   Tag_like.initialize();
       
    //loop over observations for each tag ID
   for(n=TagID_start(t); n<=TagID_end(t);n++)                 
   {
     //get position and time per observation
     int N_days_t_a=Days_libery_a(n);
     int REL_a=Rels_a(n);
     int REC_a=Recs_a(n);
      
     //re set matrix
     Move_a=Mov_mat;

     //calculate position after N_days
     if(N_days_t_a>1)
     {
       for(y=2; y<=N_days_t_a;y++) Move_a=Move_a*Mov_mat;
     }
     dvariable Pred_Prob_a=Move_a(REL_a,REC_a);
     Tag_like += - log(Pred_Prob_a);      //negative log likelihood of observations     
    }

    //Calculate likelihood of each tag ID
   // Tag_NLL_a+=Tag_like; 
    Tag_NLL_a+=(Tag_like/nATag_obs(t));   //weight by number of observations per tagID 
   }


FUNCTION Calc_obj_function
 f = Tag_NLL_a + Tag_NLL_c + 1e6*fpen_tag;


 // cout<<"Tag_NLL_a "<< Tag_NLL_a<<endl;
 // cout<<"Tag_NLL_c "<<Tag_NLL_c<<endl;
 // cout<<"fpen_tag "<<1e6*fpen_tag<<endl;



REPORT_SECTION
  //REPORT(Pred_Prob);               //predicted probabilities
