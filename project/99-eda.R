
library(tidyverse)

data <- 'data/epl.csv' %>% read_csv()
data

# tms <- data %>% tetidy::pull_distinctly(tm_h)
# tms
# 
# seasons <- data %>% tetidy::pull_distinctly(season)
# seasons

data_stack <-
  bind_rows(
    data %>%
      mutate(tm = tm_h, g = g_h),
    data %>%
      mutate(tm = tm_a, g = g_a)
  )
data_stack
# data_stack %>% 
#   # group_by(tm) %>% 
#   count(tm, season, sort = T) %>% 
#   filter(n < 38)

data_fct <-
  data %>% 
  mutate_at(vars(matches('^tm_')), as.factor)
data_fct

fit_g_partially <- purrr::partial(glm, data = data_fct, family = 'poisson', ... = )
fit_g_h <- formula(g_h ~ tm_h + tm_a) %>% fit_g_partially()
fit_g_a <- formula(g_a ~ tm_h + tm_a) %>% fit_g_partially()
fit_g_h %>% broom::tidy()
fit_g_a %>% broom::tidy()

fit_g_diff <-
  data_fct %>% 
  glm(formula(g_diff ~ tm_h + tm_a), data = ., family = 'poisson')
fit_g_diff

fit_h <-
  data_fct %>% 
  glm(formula(g_h ~ tm_h), data = ., family = 'poisson')
fit_h
fit_h %>% 
  broom::tidy() %>% 
  mutate_at(vars(estimate), list(estimate_trans = ~exp(.))) %>% 
  arrange(-estimate_trans)
fit_h %>% 
  broom::augment(type.predict = 'response') %>% 
  filter(tm_h == 'Man City')

fit_a <-
  data_fct %>% 
  glm(formula(g_a ~ tm_a), data = ., family = 'poisson')
fit_a

fits_by_tm_nested <-
  data_stack %>% 
  select(tm, g) %>% 
  group_by(tm) %>% 
  nest() %>% 
  mutate(fit = purrr::map(data, ~glm(formula(g ~ 1), data = .x, family = 'poisson')))
fits_nested

terms_by_tm_nested <-
  fits_nested %>% 
  select(-data) %>% 
  mutate(terms = purrr::map(fit, ~broom::tidy(.x))) %>% 
  select(-fit)

terms_by_tm_nested %>%
  unnest(terms) %>% 
  arrange(desc(estimate)) %>% 
  mutate(est_trans = exp(estimate)) %>% 
  mutate_at(vars(tm), ~forcats::fct_reorder(., est_trans)) %>% 
  ggplot() +
  aes(x = tm, y = est_trans) +
  geom_col() +
  coord_flip()
  
summ_by_tm_nested <-
  fits_nested %>% 
  select(-data) %>% 
  mutate(summ = purrr::map(fit, ~broom::glance(.x))) %>% 
  select(-fit)
summ_by_tm_nested
summ_by_tm_nested %>% unnest(summ)
  