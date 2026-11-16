if (!exists("save_panel")) source("code/01_plot_abundance.R")
diversity <- read_plot_data("figure1C_diversity.csv")
annotations <- read_plot_data("figure1C_annotations.csv")
annotations$Kingdom <- factor(annotations$Kingdom, c("Bacteria", "Fungi", "Archaea"))
annotations$x_start <- 1
annotations$x_end <- 2
annotations$bracket_y <- 6.4
diversity$Kingdom <- factor(diversity$Kingdom, c("Bacteria", "Fungi", "Archaea"))
diversity$Period <- factor(diversity$Period, c("Pre", "Post"))
diversity$Fill <- paste(diversity$Kingdom, diversity$Period)
colours <- c("Bacteria Pre" = "#F2CD89", "Bacteria Post" = "#B56A3A",
             "Fungi Pre" = "#CFDD97", "Fungi Post" = "#929A3E",
             "Archaea Pre" = "#C6E3F3", "Archaea Post" = "#1488B7")
# Labels are saved model results. No hypothesis test is run here.
pC <- ggplot(diversity, aes(Period, Shannon)) +
  geom_violin(aes(fill = Fill), trim = FALSE, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.13, fill = "white", linewidth = 0.25, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.18, seed = 16), size = 0.15, alpha = 0.7) +
  geom_segment(data = annotations, aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = annotations, aes(x = 1.5, y = 7, label = Label),
            inherit.aes = FALSE, size = 2.5) +
  facet_wrap(~Kingdom, nrow = 1) + scale_fill_manual(values = colours, guide = "none") +
  scale_y_continuous(breaks = c(0, 4, 8)) + coord_cartesian(ylim = c(-3.5, 8)) +
  labs(x = NULL, y = "Shannon index", title = "(C)") + panel_theme()
save_panel(pC, "figure1C", 90, 40)
