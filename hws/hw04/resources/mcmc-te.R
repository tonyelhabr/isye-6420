
# prep ----
# Data.
n_obs <- 100
set.seed(2)
y <- 2 * rnorm(n_obs) + 1 
n_mcmc <- 10000 
y_sum <- sum(y)

# Hypterparamters.
mu_0 <- 0.5
tau_0 <- 1/100 
a_0 <- 1/2 
b_0 <-  2 

mu_tau_0 <- mu_0 * tau_0

# Initial values.
mu_i <- 0.5
tau_i <- 0.5 

# main ----
# Output.
vec_num_n <- vector(mode = 'numeric', length = nn)
mu <- vec_num_n
tau <- vec_num_n

compute_mu_new <- function(mu_i, tau_i, mu_0, tau_0, y_sum, n_obs) {
   mu_tau_0 <- mu_0 * tau_0
   mu_rnorm_num <- tau_i * y_sum + mu_tau_0
   mu_rnorm_den <- tau_0 + n_obs * tau_i
   mu_rnorm <- mu_rnorm_num / mu_rnorm_den
   sigma2_rnorm <- 1 / (tau_0 + n_obs * tau_i)
   sigma_rnorm <- sqrt(sigma2_rnorm)
   rnorm(1, mu_rnorm, sigma_rnorm)
}

compute_tau_new <- function(mu_new, y, a_0, b_0, n_obs) {
   shape_rgamma <- a_0 + 0.5 * n_obs
   rate_rgamma  <- b_0 + 0.5 * sum((y - mu_new) ^ 2)
   rgamma(1, shape = shape_rgamma, rate = rate_rgamma) 
}

for (i in 1:n_mcmc) {
   # mu_rnorm_num <- tau_i * y_sum + mu_tau_0
   # mu_rnorm_den <- tau_0 + n_obs * tau_i
   # mu_rnorm <- mu_rnorm_num / mu_rnorm_den
   # sigma2_rnorm <- 1 / (tau_0 + n_obs * tau_i)
   # sigma_rnorm <- sqrt(sigma2_rnorm)
   # mu_new <- rnorm(1, mu_rnorm, sigma_rnorm)
   mu_new <-
      compute_mu_new(
         mu_i = mu_i,
         tau_i = tau_i,
         mu_0 = mu_0,
         tau_0 = tau_0,
         y_sum = y_sum,
         n_obs = n_obs
      )
   
   # shape_rgamma <- a_0 + 0.5 * n_obs
   # rate_rgamma  <- b_0 + 0.5 * sum((y - mu_new) ^ 2)
   # tau_new <- rgamma(1, shape = shape_rgamma, rate = rate_rgamma)
   tau_new <-
      compute_tau_new(
         mu_new = mu_new,
         y = y,
         a_0 = a_0,
         b_0 = b_0,
         n_obs = n_obs
      )
   
   mu[i] <- mu_new
   tau[i] <- tau_new
   mu_i <- mu_new
   tau_i <- tau_new
}

# postprocess ----
n_burnin <- 200 
idx_final <- n_burnin:n_mcmc
mu_final <- mu[idx_final] 
tau_final <- tau[idx_final] 

mean(mu_final) # 0.9404001
mean(tau_final) # 0.186471
mean(1 / tau_final)
