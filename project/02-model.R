
library(tidyverse)

data <- 'data/epl.csv' %>% read_csv()
data

tms <- data %>% tetidy::pull_distinctly(tm_h)
tms

seasons <- data %>% tetidy::pull_distinctly(season)
seasons

pull2 <- function(data, ...) {
  data %>%
    pull(...) %>% 
    as.factor() %>% 
    as.integer()
}

data_list <-
  list(
    g_h = data %>% pull(g_h),
    g_a = data %>% pull(g_a),
    tm_h = data %>% pull2(tm_h),
    tm_a = data %>% pull2(tm_a),
    season = data %>% pull2(season),
    n_tm = tms %>% length(),
    n_gm = data %>% nrow(),
    n_season = seasons %>% length()
  )
data_list

model <- glue::glue_collapse('model {
  for(g in 1:n_gm) {
    g_h[g] ~ dpois(lambda_h[season[g], tm_h[g], tm_a[g]])
    g_a[g] ~ dpois(lambda_a[season[g], tm_h[g], tm_a[g]])
  }
  
  for(s in 1:n_season) {
    for(h in 1:n_tm) {
      for(a in 1:n_tm) {
        lambda_h[s, h, a] <- exp(lvl_h[s] + z[s, h] - z[s, a])
        lambda_a[s, h, a] <- exp(lvl_a[s] + z[s, a] - z[s, h])
      }
    }
  }
    
  z[1, 1] <- 0 
  for(t in 2:n_tm) {
    z[1, t] ~ dnorm(z_all, tau_all)
  }
    
  z_all ~ dnorm(0, 0.0625)
  tau_all <- 1 / pow(sigma_all, 2)
  sigma_all ~ dunif(0, 3)
  
  lvl_h[1] ~ dnorm(0, 0.0625)
  lvl_a[1] ~ dnorm(0, 0.0625)
    
  for(s in 2:n_season) {
    z[s, 1] <- 0 
    for(t in 2:n_tm) {
      z[s, t] ~ dnorm(z[s - 1, t], tau_s)
    }
    lvl_h[s] ~ dnorm(lvl_h[s - 1], tau_s)
    lvl_a[s] ~ dnorm(lvl_a[s - 1], tau_s)
  }
  
  tau_s <- 1 / pow(sigma_s, 2) 
  sigma_s ~ dunif(0, 3) 
}')


path_model <- 'model.txt'
write_lines(model, path_model)

path_res_sim <- 'output/res_sim.rds'
path_res_sim_jags <- 'output/res_sim_jags.rds'
if(fs::file_exists(path_res_sim)) {
  
} else {
  
  inits <- NULL
  params <-
    c(
      paste0('lvl_', c('h', 'a')),
      paste0(c('', 'grp_'), 'z'),
      paste0(c('grp_', 's_'), 'sigma')
    )
  # params <- 
  #   c(
  #     paste0('lvl_', c('h', 'a'))
  #   )
  params
  
  res_sim <-
    R2OpenBUGS::bugs(
      # debug = TRUE,
      data = data_list,
      inits = inits,
      model.file = path_model,
      parameters.to.save = params,
      DIC = FALSE,
      n.chains = 1,
      n.iter = 10000,
      n.burnin = 1000
    )
  res_sim$summary
  
  teproj::export_path(res_sim, path_res_sim)
}

if(FALSE) {
  model_jags <-
    rjags::jags.model(textConnection(model), data = data_list, n.chains = 1, n.adapt = 5000)
  update(model_jags, 10000)
  res_sim_summ_jags <-
    coda::coda.samples(model_jags, variable.names = params, n.iter = 10000, thin = 2)
  res_sim_summ_jags <-
    res_sim_jags %>% 
    as.matrix() %>% 
    as_tibble()

  res_sim_summ_jags
  teproj::export_path(res_sim, path_res_sim_jags)
}


res_sim_summ <-
  res_sim$summary %>% 
  as.matrix() %>% 
  # as_tibble()
  as_tibble(rownames = 'row') %>% 
  remove_rownames() %>% 
  rename_at(vars(matches('[%]')), ~str_replace_all(., '(^.*)([%]$)', '\\1') %>% paste0('q', .))
res_sim_summ

# plot(res_sim$sims.list$lvl_h)
hist(res_sim$sims.list$lvl_h)
hist(res_sim$sims.list$sigma_s)
