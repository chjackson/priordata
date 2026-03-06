functions {
  // Model 1 for human growth from Preece and Baines (1978)
  vector PB1(vector age, real h0, real h1, real s0, real s1, real gamma){
    return h1 - 2.*(h1 - h0) ./ (exp(s0.*(age - gamma)) + exp(s1.*(age - gamma)));
  }
}
data {
  int<lower=0> n;
  vector[n] height;
  vector[n] age;

  int<lower=0> n_elic;
  vector[n_elic] age_elic;
  vector[n_elic] ybar_elic;
  vector[n_elic] se_elic;
  int<lower=0,upper=1> pred; // 1: elicited information is predictive dist for an observable. 0: uncertainty dist for a mean

  int<lower=0> n_pred;
  vector[n_pred] age_pred;

  real ms0; real<lower=0> ss0;
  real mds; real<lower=0> sds;
  real mh1; real<lower=0> sh1;
  real<lower=0> ph0a; real<lower=0> ph0b;
  real<lower=0> mingamma; real<lower=0> maxgamma;
  real msd; real<lower=0> ssd;
}

parameters {
  real h1_std;
  real s0_std;
  real ds_std;
  real<lower=0,upper=1> ph0;
  real<lower=0,upper=1> gamma_std;
  real sd_height_std;
  vector[n_elic] alpha;
}

transformed parameters {
  real ds = exp(ds_std*sds + mds);
  real h1 = exp(h1_std*sh1 + mh1);
  real h0 = h1*ph0;
  real s0 = exp(s0_std*ss0 + ms0);
  real s1 = s0 + ds;
  real gamma = mingamma + gamma_std*(maxgamma - mingamma);
  real sd_height = exp(sd_height_std*ssd + msd);
}

model {
  vector[n] mean_height;
  vector[n_elic] mean_elic;
  h1_std ~ std_normal();
  ph0 ~ beta(ph0a, ph0b);
  s0_std ~ std_normal();
  ds_std ~ std_normal();
  gamma_std ~ uniform(0, 1);
  sd_height_std ~ std_normal();

  if (n > 0){
    mean_height = PB1(age, h0, h1, s0, s1, gamma);
    height ~ normal(mean_height, sd_height);
  }

  if (n_elic > 0){
    mean_elic = PB1(age_elic, h0, h1, s0, s1, gamma);
    if (pred==0){
      ybar_elic ~ normal(mean_elic, se_elic); 
    } else if (pred==1){
      alpha ~ normal(mean_elic, sd_height);
      ybar_elic ~ normal(alpha, se_elic);
    }
  }

}

generated quantities {
  vector[n_pred] height_pred; 
  vector[n_pred] mean_pred; 
  if (n_pred > 0){
    mean_pred = PB1(age_pred, h0, h1, s0, s1, gamma);
    for (i in 1:n_pred){
      height_pred[i] = normal_rng(mean_pred[i], sd_height);
    }
  }
}
