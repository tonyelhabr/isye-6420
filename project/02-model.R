
data <- 'data/epl.csv' %>% read_csv()
data
# tms <- data %>% tetidy::pull_distinctly(tm_h)
# tms

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
    n_season = .seasons %>% length()
  )
data_list

model <- 'model {
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
    z[1, t] ~ dnorm(grp_z, grp_tau)
  }
    
  grp_z ~ dnorm(0, 0.0625)
  grp_tau <- 1 / pow(grp_sigma, 2)
  grp_sigma ~ dunif(0, 3)
  
  lvl_h[1] ~ dnorm(0, 0.0625)
  lvl_a[1] ~ dnorm(0, 0.0625)
    
  for(s in 2:n_season) {
    z[s, 1] <- 0 
    for(t in 2:n_tm) {
      z[s, t] ~ dnorm(z[s - 1, t], s_tau)
    }
    lvl_h[s] ~ dnorm(lvl_h[s - 1], s_tau)
    lvl_a[s] ~ dnorm(lvl_a[s - 1], s_tau)
  }
  
  s_tau <- 1 / pow(s_sigma, 2) 
  s_sigma ~ dunif(0, 3) 
}'

path_model <- 'model.txt'
write_lines(model, path_model)


path_res_sim <- 'output/res_sim.rds'
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
}

teproj::export_path(res_sim, path_res_sim)

res_sim$summary %>% as.matrix()

plot(res_sim$sims.list$lvl_h)
