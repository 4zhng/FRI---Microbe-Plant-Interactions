library(tidyverse)
library(Bolstad2)

plate_map <- colr100_p560_map %>% 
  mutate(Column = as.character(Column))
lum <- colR100_p560_exp %>% 
  fill(Time) %>% 
  mutate(Row = rep(LETTERS[1:8], 31)) %>% 
  gather(Column, LUM, -c(Time, Row)) %>% 
  left_join(plate_map)

# Take mean luminescence at each time point for each treatment
lum_means <- lum %>%
  group_by(Trt, Isolate, Time) %>% 
  summarise(meanlum = mean(LUM), se = sd(LUM) / sqrt(n()))

# Make a graph to look at the curves
lum_means  %>% 
  na.omit() %>% 
  ungroup() %>% 
  mutate(sample = paste(Isolate, Trt),
         sample = fct_relevel(sample, "Control Mock", "Control Flg22", "ColR100 Alive", "ColR100 Dead", "P560 Alive", "P560 Dead"),
         Isolate = fct_relevel(Isolate, "Control", "ColR100", "P560")) %>% 
  ggplot(aes(Time, meanlum, color = sample)) +
  geom_line() +
  geom_linerange(aes(ymin = meanlum - se, ymax = meanlum + se)) +
  facet_grid(.~Isolate) +
  scale_color_brewer(palette = "Paired") +
  theme_bw()


## Area under the curve
lum %>% 
  na.omit() %>% 
  mutate(sample = paste(Isolate, Trt),
         sample = fct_relevel(sample, "Control Mock", "Control Flg22", "ColR100 Alive", "ColR100 Dead", "P560 Alive", "P560 Dead")) %>% 
  group_by(sample, Row, Column) %>% 
  summarise(auc = sintegral(Time, LUM)$int) %>% 
  ggplot(aes(sample, auc/100000, fill = sample)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes(color = sample), width = 0.1) +
  scale_fill_brewer(palette = "Paired", direction = 1) +
  scale_color_brewer(palette = "Paired", direction = 1) +
  theme_bw() +
  labs(x = "Treatment", y = "Area Under the Curve") +
  theme(axis.text.x = element_text(angle = 30, hjust =1))

## Do Stats for area under the curve
lum %>% 
  mutate(sample = paste(Isolate, Trt),
         sample = fct_relevel(sample, "Control Mock", "Control Flg22", "ColR100 Alive", "ColR100 Dead", "P560 Alive", "P560 Dead")) %>% 
  group_by(sample, Row, Column) %>% 
  summarise(auc = sintegral(Time, LUM)$int) %>% 
  ungroup() %>% 
  filter(!grepl("NA", sample)) %>% 
  aov(auc ~ sample, .) %>% 
  TukeyHSD() %>% 
  broom::tidy()