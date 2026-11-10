library(dplyr)

tp_levels <- c("T0","T1","T2","T3")

state_all <- dat_long %>%
  mutate(Timepoint = factor(as.character(Timepoint),
                            levels = tp_levels,
                            ordered = TRUE)) %>%
  filter(!is.na(Timepoint), Timepoint %in% tp_levels)

tab_all <- table(state_all$Timepoint, state_all$Assigned_CST_uncertain)

tab_all

set.seed(123)
chisq_all <- chisq.test(tab_all, simulate.p.value = TRUE, B = 20000)
chisq_all

prop_tab_all <- prop.table(tab_all, margin = 1) 
round(prop_tab_all, 3)

library(ggplot2)

prop_df <- as.data.frame(prop_tab_all)
colnames(prop_df) <- c("Timepoint","CST","Proportion")
dev.new()
ggplot(prop_df, aes(x = Timepoint, y = Proportion, fill = CST)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 12) +
  labs(title = "CST distribution across timepoints",
       y = "Proportion")

