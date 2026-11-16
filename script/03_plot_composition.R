if (!exists("save_panel")) source("code/01_plot_abundance.R")
library(ggalluvial)
relative <- read_plot_data("figure1DE_composition.csv")
absolute <- read_plot_data("figure1FG_abundance.csv")
palette <- read_plot_data("figure1_palette.csv")
plot_labels <- read_plot_data("figure1_annotations.csv")

composition_panel <- function(values, kingdom, group, is_relative) {
  d <- values[values$Kingdom == kingdom & values$Group == group, ]
  pal <- palette[palette$Kingdom == kingdom, ]
  pal <- pal[order(pal$Bottom_to_top), ]
  d$Genus <- factor(d$Genus, levels = rev(pal$Genus))
  times <- if (group == "w/o") c("T0", "T1", "T3") else c("T0", "T1", "T2", "T3")
  d$Timepoint <- factor(d$Timepoint, levels = times)
  p <- ggplot(d, aes(Timepoint, Value, stratum = Genus, alluvium = Genus, fill = Genus)) +
    geom_flow(width = 0.5, alpha = 0.45, colour = "transparent") +
    geom_stratum(width = 0.5, colour = NA) +
    scale_fill_manual(values = setNames(pal$Color, pal$Genus),
                      breaks = pal$Genus[order(pal$Legend_order)], drop = FALSE) +
    scale_x_discrete(expand = expansion(add = 0.25)) +
    labs(x = NULL, y = if (is_relative) "Relative abundance" else "Absolute abundance",
         title = paste(kingdom, group), fill = NULL) + panel_theme() +
    theme(plot.title = element_text(hjust = 0.5, size = 7),
          legend.text = element_text(size = 4.3, face = "italic"),
          legend.key.height = grid::unit(1.6, "mm"),
          legend.key.width = grid::unit(1.6, "mm"), legend.spacing.y = grid::unit(0, "mm"))
  if (is_relative) p <- p + scale_y_continuous(breaks = seq(0, 1, 0.25),
                                                limits = c(0, 1.000001), expand = c(0, 0))
  else p <- p + scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                                   breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0, 0.04)))
  label <- plot_labels[plot_labels$Kingdom == kingdom & plot_labels$Group == group &
                         plot_labels$Measure == if (is_relative) "Relative" else "Absolute", ]
  p + geom_text(data = label, aes(Text_x, Text_y, label = Label), inherit.aes = FALSE,
                fontface = "italic", size = 1.8, hjust = 0) +
    geom_segment(data = label, aes(x = Arrow_x, y = Arrow_y,
                  xend = Arrow_end_x, yend = Arrow_end_y), inherit.aes = FALSE,
                  colour = "grey45", linewidth = 0.2,
                  arrow = grid::arrow(length = grid::unit(1, "mm")))
}
paired_panels <- function(values, kingdom, is_relative, tag) {
  left <- composition_panel(values, kingdom, "w/o", is_relative)
  right <- composition_panel(values, kingdom, "w/", is_relative) + labs(y = NULL)
  (left | right) + plot_layout(guides = "collect") +
    plot_annotation(title = tag)
}
pD <- paired_panels(relative, "Bacteria", TRUE, "(D)")
pE <- paired_panels(relative, "Fungi", TRUE, "(E)")
pF <- paired_panels(absolute, "Bacteria", FALSE, "(F)")
pG <- paired_panels(absolute, "Fungi", FALSE, "(G)")
for (letter in c("D", "E", "F", "G")) save_panel(get(paste0("p", letter)), paste0("figure1", letter), 100, 47)
