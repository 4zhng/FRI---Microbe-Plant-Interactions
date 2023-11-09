library(tidyverse)
library(growthcurver)

## The map of which isolates are in which wells
isolate_locations <- data.frame(column = as.character(1:12), isolates = c("P500", "KBS25", "P5", "AAF232" ,"ColR100", "P560", "P3A6", "P3B10", "275T1", "AAF210", "GP1C1", "P2B5"))

# Stats taken from the software of the plate reader machine. We really care about the time at max velocity
growth_stats <- growthstats %>% 
  mutate(tmax = period_to_seconds(hms(tatMaxV)) / 60 / 60) %>% 
  mutate(column = gsub("[A-Z]", "", Well),
         row = gsub("[0-9]+", "", Well)) %>% 
  left_join(isolate_locations) %>% 
  mutate(isolates = ifelse(row == "C", "Control", isolates)) 

# These are controls
  gc_control <- growth_curves %>% 
  select_if(grepl("C", names(.))) 

# These are the inoculated wells. Need to manipulate the time column to actually be readable by R
  gc <- growth_curves %>% 
  select_if(!grepl("C", names(.))) %>% 
  as_tibble() %>% 
  mutate(Time = period_to_seconds(hms(Time)) / 60 / 60) %>% 
  rename(time = Time) %>% 
  mutate(blank = rowMeans(gc_control)) %>% 
  as.data.frame()

# Tidy up the data to make it long
gc_tidy <- gc %>% 
  gather(Well, OD600, -time) %>% 
  mutate(column = gsub("[ABCDEFGH]", "", Well)) %>% 
  inner_join(isolate_locations) %>% 
  mutate(type = ifelse(grepl("C", Well), "Control", "Inoculated"))

# Plot out growth curves for things less than 40 hours
gc_tidy %>%
  filter(time < 40) %>% 
  ggplot(aes(time, OD600, group = Well)) +
  geom_line() +
  facet_grid(.~isolates) +
  geom_vline(data = growth_stats %>% filter(isolates != "Control" & tmax < 40), aes(xintercept = tmax, color = isolates)) +
  theme_bw()




         