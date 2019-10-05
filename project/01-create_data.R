
library(tidyverse)

# data_raw_old <- 'https://raw.githubusercontent.com/jalapic/engsoccerdata/2a8f61005f0b89e797e8fbecd063f36fe394c2c5/data-raw/england.csv' %>% read_csv()
# data_raw_old %>% count(Season) %>% arrange(-Season)
# 
# 
# data_raw <- 'https://raw.githubusercontent.com/jalapic/engsoccerdata/master/data-raw/england.csv' %>% read_csv()
# data_raw
# data_filt <-
#   data_raw_old %>%
#   janitor::clean_names() %>%
#   filter(season >= 2015, season <= 2016) %>%
#   filter(division == 1) %>%
#   rename(
#     tm_h = home,
#     tm_v = visitor,
#     g_h = hgoal,
#     g_v = vgoal,
#     g_total = totgoal,
#     g_diff = goaldif
#   )
# data_filt
# data_filt %>% count(season)
# ? engsoccerdata::england

# Reference: https://github.com/jalapic/engsoccerdata/blob/master/R/england_current.R
scrape_epl_data <- function(season = lubridate::year(Sys.Date()) - 1L) {
  # season <- 2018L
  s1 <- season %>% str_sub(3, 4) %>% as.integer()
  s2 <- s1 + 1L
  path <- sprintf('http://www.football-data.co.uk/mmz4281/%2d%2d/E0.csv', s1, s2)
  data_raw <- path %>% read_csv()
  data <-
    data_raw %>% 
    janitor::clean_names() %>% 
    mutate_at(vars(date), ~as.Date(., '%d/%m/%y')) %>% 
    select(
      date, 
      tm_h = home_team, 
      tm_a = away_team,
      g_h = fthg,
      g_a = ftag
    ) %>% 
    mutate(
      g_total = g_h + g_a,
      g_diff = g_h - g_a,
      result = 
        case_when(
          g_h > g_a ~ 'h', 
          g_h < g_a ~ 'a', 
          TRUE ~ 't'
        ),
      tm_winner = 
        case_when(
          g_h > g_a ~ tm_h, 
          g_h < g_a ~ tm_a, 
          TRUE ~ NA_character_
        )
    ) %>% 
    mutate_at(vars(matches('season|^g_')), as.integer)
}

.seasons <- 2015L:2017L
data <-
  tibble(season = .seasons) %>% 
  mutate(data = purrr::map(season, scrape_epl_data)) %>% 
  unnest(data)
data
teproj::export_path(data, 'data/epl.csv')

