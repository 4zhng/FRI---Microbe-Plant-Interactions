install.packages("tidyverse")
install.packages("Bolstad2")

library(tidyverse)
library(Bolstad2)

ros_data <- 
  
  
install.packages("tidyverse")
install.packages("Bolstad2")

library(tidyverse)
library(Bolstad2)

# Down load ROS data from
ros_data <- read_tsv("~/Library/CloudStorage/Box-Box/FRI/ImmunityGroup2023/leaf_zone.txt") %>% 
  fill(Time, Temperature)
plate_map <- read_tsv("~/Library/CloudStorage/Box-Box/FRI/ImmunityGroup2023/PlateMap.txt")

ros_data_tidy <- ros_data %>% 
  select(-Temperature) %>% # remove temperature column
  mutate(Row = rep(LETTERS[1:8], 21)) %>% # Adds row labels
  gather(Column, LUM, -c(Row, Time)) %>% 
  mutate(Column = as.numeric(Column)) %>% 
  mutate(well = paste0(Row, Column)) %>% 
  inner_join(plate_map) %>% 
  mutate(Zone = fct_relevel(Zone, "Tip", "Middle", "Stem"))

leaf_colors <- RColorBrewer::brewer.pal(4, "BuPu")[2:4]


ggplot(ros_data_tidy, aes(x = Time, y = LUM, group = well, color = Zone)) +
  geom_line() +
  facet_grid(Leaf ~ Trt) +
  scale_color_manual(values = leaf_colors) +
  theme_bw()

# Make graphs with median (or mean) of each treatmentment level at each timepoint
ros_data_tidy %>% 
  group_by(Time, Zone, Trt, Leaf) %>% 
  summarise(meanlum = median(LUM), se = sd(LUM)/sqrt(n())) %>% 
  ggplot(aes(Time, meanlum, color = Zone)) +
  geom_line() +
  facet_grid(Leaf ~ Trt) +
  scale_color_manual(values = leaf_colors) +
  theme_bw()

## Calculate are under the curve
ros_data_tidy %>% 
  group_by(Zone, Trt, Leaf, well) %>% 
  summarise(auc = sintegral(Time, LUM)$int) %>% 
  ggplot(aes(Zone, auc, fill = Trt)) +
  geom_boxplot() +
  facet_grid(.~Leaf)

## Run anova on the AUC
ros_data_tidy %>% 
  group_by(Zone, Trt, Leaf, well) %>% 
  filter(Trt == "Flg22") %>% 
  summarise(auc = sintegral(Time, LUM)$int) %>% 
  ungroup() %>% 
  aov(auc ~ Zone*Leaf, .) %>% 
  broom::tidy() %>% 
  mutate(percent_var = sumsq / sum(sumsq) * 100)