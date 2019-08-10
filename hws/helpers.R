
library(tidyverse)
create_kable_md <- function (data) {
  knitr::kable(data, format = 'markdown', escape = FALSE)
}

convert_page_to_lines <- function(page) {
  page %>%
    str_split("\\n") %>%
    purrr::map(str_squish) %>%
    unlist() %>%
    enframe(name = "idx_line", value = "line")
}

path <- 'hws/Homework1f.pdf'
text_raw <- path %>% pdftools::pdf_text()
text_raw
lines_raw <- text_raw %>% convert_page_to_lines()
lines_raw