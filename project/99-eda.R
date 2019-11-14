
library(tidyverse)

data <- 'data/epl.csv' %>% read_csv()
data

tms <- data %>% tetidy::pull_distinctly(tm_h)
tms

seasons <- data %>% tetidy::pull_distinctly(season)
seasons

data_stack <-
  bind_rows(
    data %>% 
      mutate(tm = tm_h, g = g_h),
    data %>% 
      mutate(tm = tm_a, g = g_a)
  )
data_stack
data_stack %>% 
  # group_by(tm) %>% 
  count(tm, season, sort = T) %>% 
  filter(n < 38)

fits_nested <-
  data_stack %>% 
  select(tm, g) %>% 
  group_by(tm) %>% 
  nest() %>% 
  mutate(fit = purrr::map(data, ~glm(formula(g ~ 1), data = .x, family = 'poisson')))
fits_nested

terms_nested <-
  fits_nested %>% 
  select(-data) %>% 
  mutate(terms = purrr::map(fit, ~broom::tidy(.x))) %>% 
  select(-fit)

terms_nested %>%
  unnest(terms) %>% 
  arrange(desc(estimate)) %>% 
  mutate(est_trans = exp(estimate)) %>% 
  mutate_at(vars(tm), ~forcats::fct_reorder(., est_trans)) %>% 
  ggplot() +
  aes(x = tm, y = est_trans) +
  geom_col() +
  coord_flip()
  

summ_nested <-
  fits_nested %>% 
  select(-data) %>% 
  mutate(summ = purrr::map(fit, ~broom::glance(.x))) %>% 
  select(-fit)
summ_nested
summ_nested %>% 
  unnest(summ)
  