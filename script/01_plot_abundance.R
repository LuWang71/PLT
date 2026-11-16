library(ggplot2)
library(patchwork)

# Run from the Figure1 folder, or pass that folder to Rscript.
if (!exists("figure1_root")) {
  args <- commandArgs(trailingOnly = TRUE)
  figure1_root <- if (length(args)) normalizePath(args[1]) else normalizePath(".")
}
data_dir <- file.path(figure1_root, "data")
output_dir <- file.path(figure1_root, "plots")
dir.create(output_dir, showWarnings = FALSE)
read_plot_data <- function(name) read.csv(file.path(data_dir, name), check.names = FALSE)
save_panel <- function(plot, name, width, height) {
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "mm", dpi = 600,
         device = ragg::agg_png, bg = "white")
}
panel_theme <- function() theme_classic(base_size = 7, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.35),
        axis.line = element_blank(), axis.text = element_text(colour = "black"),
        axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = grid::unit(1, "mm"),
        strip.background = element_blank(), strip.text = element_text(size = 7),
        plot.margin = margin(2, 2, 2, 2))

abundance <- read_plot_data("figure1B_abundance.csv")
abundance$Kingdom <- factor(abundance$Kingdom, c("Bacteria", "Fungi", "Archaea"))
pB <- ggplot(abundance, aes(Day, Log10Abundance)) +
  geom_point(size = 0.18, alpha = 0.45, colour = "grey30") +
  facet_wrap(~Kingdom, nrow = 1) +
  scale_x_continuous(breaks = seq(0, 40, 10)) +
  scale_y_continuous(breaks = seq(0, 12.5, 2.5)) +
  coord_cartesian(xlim = c(-1, 40), ylim = c(0, 12.5)) +
  labs(x = "Day", y = expression(Log[10]~"(Total abundance)"), title = "(B)") +
  panel_theme()
save_panel(pB, "figure1B", 90, 40)
