library(tidyverse)
library(Bolstad2)

exp1_plate_map <- colr100_p560_map %>% 
  mutate(Column = as.character(Column))
exp1_lum <- colR100_p560_exp %>% 
  fill(Time) %>% 
  mutate(Row = rep(LETTERS[1:8], 31)) %>% 
  gather(Column, LUM, -c(Time, Row)) %>% 
  left_join(exp1_plate_map) %>% 
  mutate(experiment = "Exp1")



exp2_plate_map <- GP1C1_AAF232_AAF210_map %>% 
  mutate(Column = as.character(Column))
exp2_lum <- GP1C1_AAF232_AAF210_exp %>% 
  fill(Time) %>% 
  dplyr::select(-Temperature) %>% 
  mutate(Row = rep(LETTERS[1:8], 31)) %>% 
  gather(Column, LUM, -c(Time, Row)) %>% 
  left_join(exp2_plate_map)%>% 
  mutate(experiment = "Exp2")


exp3_plate_map <- P3B10_P500_P3A6_map %>% 
  mutate(Column = as.character(Column))
exp3_lum <- P3B10_P500_P3A6_exp %>% 
  fill(Time) %>% 
  dplyr::select(-Temperature) %>% 
  mutate(Row = rep(LETTERS[1:8], 31)) %>% 
  gather(Column, LUM, -c(Time, Row)) %>% 
  left_join(exp3_plate_map)%>% 
  mutate(experiment = "Exp3")



# Take mean luminescence at each time point for each treatment
lum_means <- bind_rows(exp1_lum, exp2_lum, exp3_lum) %>%
  group_by(Trt, Isolate, Time, experiment) %>% 
  summarise(meanlum = mean(LUM), se = sd(LUM) / sqrt(n()))

# Make a graph to look at the curves
lum_means  %>% 
  na.omit() %>% 
  ungroup() %>% 
  mutate(sample = paste(Isolate, Trt),
         Isolate = fct_relevel(Isolate, "Control", "GP1C1", "AAF232", "P560", "AAF210", "ColR100", "P3B10", "P500", "P3A6")) %>% 
  ggplot(aes(Time, meanlum, color = Trt)) +
  geom_line(size = 1) +
  geom_linerange(aes(ymin = meanlum - se, ymax = meanlum + se)) +
  facet_grid(.~Isolate + experiment) +
  scale_color_manual(values = c("black", "grey70", "red", "dodgerblue")) +
  theme_bw()


## Area under the curve
bind_rows(exp1_lum, exp2_lum, exp3_lum) %>%
  na.omit() %>% 
  mutate(sample = paste(Isolate, Trt),
         sample = fct_relevel(sample, "Control Mock", "Control Flg22", "GP1C1 Alive", "GP1C1 Dead", "AAF232 Alive", "AAF232 Dead", "AAF210 Alive", "AAF210 Dead", "P3B10 Alive", "P3B10 Dead", "P500 Alive", "P500 Dead", "P3A6 Alive", "P3A6 Dead"),
         Isolate = fct_relevel(Isolate, "Control", "GP1C1", "AAF232", "P560", "AAF210", "P500", "P3B10", "ColR100", "P3A6")) %>% 
  group_by(sample, Row, Column, Isolate, Trt) %>% 
  summarise(auc = sintegral(LUM, Time)$int) %>% 
  ggplot(aes(sample, auc/100000, fill = Trt)) +
  geom_boxplot(alpha = 0.8, outlier.color = "white") +
  geom_jitter(width = 0.1, color = "black", size = 0.5) +
  #scale_color_manual(values = c( "magenta2", "grey50", "#EAC451", "#C3822C", "#A7C4E5", "#295278",  "#DE9D9B", "#E05241", "#bcd6ac", "#9dcc3f",  "#cea8bc", "#B82674")) +
  theme_bw() +
  scale_fill_manual(values = c("black", "grey50", "red", "dodgerblue")) +
  labs(x = "Treatment", y = "Area Under the Curve") +
  facet_grid(.~Isolate, scales = "free") +
  theme(axis.text.x = element_text(angle = 30, hjust =1))

## Do Stats for area under the curve
# First way is to do linear models
bind_rows(exp1_lum, exp2_lum, exp3_lum)  %>% 
  mutate(sample = paste(Isolate, Trt)) %>% 
  group_by(sample, Row, Column) %>% 
  summarise(auc = sintegral(LUM, Time)$int) %>% 
  ungroup() %>% 
  filter(!grepl("NA", sample)) %>% 
  lm(log10(auc) ~ sample, .) %>% 
  emmeans::emmeans(., specs = pairwise ~ sample, adjust = "none") %>% 
  .$contrasts %>% 
  broom::tidy() %>% 
  filter(grepl("Control.*Control", contrast) | grepl("AAF232.*AAF232", contrast) |grepl("AAF210.*AAF210", contrast) | grepl("GP1C1.*GP1C1", contrast) | grepl("ColR100.*ColR100", contrast) | grepl("P560.*P560", contrast) | grepl("P3B10.*P3B10", contrast) | grepl("P500.*P500", contrast) | grepl("P3A6.*P3A6", contrast)) %>% 
  mutate(padj = p.adjust(p.value))

# Second way is to do t-test
bind_rows(exp1_lum, exp2_lum, exp3_lum)  %>% 
  mutate(sample = paste(Isolate, Trt)) %>% 
  group_by(sample, Row, Column) %>% 
  summarise(auc = sintegral(LUM, Time)$int) %>% 
  ungroup() %>% 
  filter(!grepl("NA", sample)) %>% 
  separate(sample, into = c("Isolate", "Status"), sep = " ", remove = F) %>% 
  group_by(Isolate) %>% nest() %>% 
  mutate(mod = map(data, ~broom::tidy(t.test(log10(auc) ~ Status, .)))) %>% 
  unnest(mod) %>% ungroup() %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  View()





#EXPERIMENT 1 AND 2 ONLY (3rd experiment had issue with plant tissue and control)
 
library(tidyverse)
library(Bolstad2)

exp1_plate_map <- colr100_p560_map %>% 
  mutate(Column = as.character(Column))
exp1_lum <- colR100_p560_exp %>% 
  fill(Time) %>% 
  mutate(Row = rep(LETTERS[1:8], 31)) %>% 
  gather(Column, LUM, -c(Time, Row)) %>% 
  left_join(exp1_plate_map) %>% 
  mutate(experiment = "Exp1")



exp2_plate_map <- GP1C1_AAF232_AAF210_map %>% 
  mutate(Column = as.character(Column))
exp2_lum <- GP1C1_AAF232_AAF210_exp %>% 
  fill(Time) %>% 
  dplyr::select(-Temperature) %>% 
  mutate(Row = rep(LETTERS[1:8], 31)) %>% 
  gather(Column, LUM, -c(Time, Row)) %>% 
  left_join(exp2_plate_map)%>% 
  mutate(experiment = "Exp2")


# Take mean luminescence at each time point for each treatment
lum_means <- bind_rows(exp1_lum, exp2_lum) %>%
  group_by(Trt, Isolate, Time, experiment) %>% 
  summarise(meanlum = mean(LUM), se = sd(LUM) / sqrt(n()))

# Make a graph to look at the curves
lum_means  %>% 
  na.omit() %>% 
  ungroup() %>% 
  mutate(sample = paste(Isolate, Trt),
         Isolate = fct_relevel(Isolate, "Control", "GP1C1", "AAF232", "P560", "AAF210", "ColR100")) %>% 
  ggplot(aes(Time, meanlum, color = Trt)) +
  geom_line(size = 1) +
  geom_linerange(aes(ymin = meanlum - se, ymax = meanlum + se)) +
  facet_grid(.~Isolate + experiment) +
  scale_color_manual(values = c("black", "grey70", "red", "dodgerblue")) +
  theme_bw()


## Area under the curve
bind_rows(exp1_lum, exp2_lum) %>%
  na.omit() %>% 
  mutate(sample = paste(Isolate, Trt),
         sample = fct_relevel(sample, "Control Mock", "Control Flg22", "GP1C1 Alive", "GP1C1 Dead", "AAF232 Alive", "AAF232 Dead", "AAF210 Alive", "AAF210 Dead"),
         Isolate = fct_relevel(Isolate, "Control", "GP1C1", "AAF232", "P560", "AAF210","ColR100")) %>% 
  group_by(sample, Row, Column, Isolate, Trt) %>% 
  summarise(auc = sintegral(LUM, Time)$int) %>% 
  ggplot(aes(sample, auc/100000, fill = Trt)) +
  geom_boxplot(alpha = 0.8, outlier.color = "white") +
  geom_jitter(width = 0.1, color = "black", size = 0.5) +
  #scale_color_manual(values = c( "magenta2", "grey50", "#EAC451", "#C3822C", "#A7C4E5", "#295278",  "#DE9D9B", "#E05241", "#bcd6ac", "#9dcc3f",  "#cea8bc", "#B82674")) +
  theme_bw() +
  scale_fill_manual(values = c("black", "grey50", "red", "dodgerblue")) +
  labs(x = "Treatment", y = "Area Under the Curve") +
  facet_grid(.~Isolate, scales = "free") +
  theme(axis.text.x = element_text(angle = 30, hjust =1))

## Do Stats for area under the curve
# First way is to do linear models
bind_rows(exp1_lum, exp2_lum)  %>% 
  mutate(sample = paste(Isolate, Trt)) %>% 
  group_by(sample, Row, Column) %>% 
  summarise(auc = sintegral(LUM, Time)$int) %>% 
  ungroup() %>% 
  filter(!grepl("NA", sample)) %>% 
  lm(log10(auc) ~ sample, .) %>% 
  emmeans::emmeans(., specs = pairwise ~ sample, adjust = "none") %>% 
  .$contrasts %>% 
  broom::tidy() %>% 
  filter(grepl("Control.*Control", contrast) | grepl("AAF232.*AAF232", contrast) |grepl("AAF210.*AAF210", contrast) | grepl("GP1C1.*GP1C1", contrast) | grepl("ColR100.*ColR100", contrast) | grepl("P560.*P560", contrast)) %>% 
  mutate(padj = p.adjust(p.value))

# Second way is to do t-test
bind_rows(exp1_lum, exp2_lum)  %>% 
  mutate(sample = paste(Isolate, Trt)) %>% 
  group_by(sample, Row, Column) %>% 
  summarise(auc = sintegral(LUM, Time)$int) %>% 
  ungroup() %>% 
  filter(!grepl("NA", sample)) %>% 
  separate(sample, into = c("Isolate", "Status"), sep = " ", remove = F) %>% 
  group_by(Isolate) %>% nest() %>% 
  mutate(mod = map(data, ~broom::tidy(t.test(log10(auc) ~ Status, .)))) %>% 
  unnest(mod) %>% ungroup() %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  View()