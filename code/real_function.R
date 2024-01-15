
library(tidyverse)

# real data for cjs model ---------------------------------------

df_0 <- read_csv(here::here("data_raw/formatted_cmr_data.csv")) %>% 
  subset(species == "creek_chub") %>%   # subset out 1 species for practice 
  dplyr::select(-c(w,
                   fin_recap,
                   f_occasion,
                   section,
                   species,
                   tag,
                   length,
                   weight,
                   date,
                   '...1'))   # remove extraneous columns


df_0$recap<- as.numeric(c("y" = "1", "n" = "0")[df_0$recap]) #convert character to binomial
df_0$mortality<- as.numeric(c("y" = "1", "n" = "0")[df_0$mortality]) #convert character to binomial

recap_mat <- df_0 %>% 
  select(-mortality) %>% 
  pivot_wider(names_from = occasion,
            values_from = recap)
