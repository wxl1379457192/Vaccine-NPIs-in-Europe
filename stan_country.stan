//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
data{
  int <lower=1> m; // number of all records
  int <lower=1> k1; // number of interest variables
  int <lower=1> k2; // number of environment variables
  matrix[m, k1] X1; // interest variables
  matrix[m, k2] X2; // control variables
  vector[m] Rt;
  vector[m] R0;
}
parameters{
  real<lower=0> alpha_hier[k1]; // real priors
  real<lower=0> beta_hier[k2]; // real priors
 // real<lower=0> sigma2;
}
transformed parameters {
  vector[k1] alpha;
  vector[k2] beta;
  vector[m] mu;
//  matrix[m,k1+k2] fluctuation;
  {
  // calculation of variants-based real R0
  //  for(c in 1:n){
   //   for(i in ((c-1)*Cnum+1):(c*Cnum)){  
   //     for(j in 1:k1+k2){
   //       fluctuation[i,j] = country_variation[c,j];
   //     }
   //   }
   // } 
 
    for(i in 1:k1){
      alpha[i] = alpha_hier[i] - log(1.000382)/7;
    // alpha[i] = alpha_hier[i];
    }
    for(i in 1:k2){
      beta[i] = beta_hier[i] - log(1.09)/5;
    }
    //mu = -X1*alpha;
    mu = -X1*alpha-X2*beta;
    

  }
}


model{
//////////////////// main ///////////////////////////////////
  //alpha_hier ~ gamma(.1667,1);
  alpha_hier ~ normal(0,0.3);
  //sigma2 ~normal(0.5,1);
  beta_hier ~ gamma(.1667,1);
  for(i in 1:m){
    Rt[i] ~ gamma((exp(mu[i])*R0[i]),0.5);
  }
}
