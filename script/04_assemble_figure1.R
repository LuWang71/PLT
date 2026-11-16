# Rscript code/04_assemble_figure1.R [path/to/Figure1]
args <- commandArgs(trailingOnly = TRUE)
figure1_root <- if (length(args)) normalizePath(args[1]) else normalizePath(".")
for (script in c("01_plot_abundance.R", "02_plot_diversity.R", "03_plot_composition.R")) {
  source(file.path(figure1_root, "code", script))
}
# Panel A is the final, manually drawn study diagram, supplied as a 600-dpi asset.
diagram <- png::readPNG(file.path(data_dir, "figure1A_diagram.png"))
pA <- wrap_elements(full = grid::rasterGrob(diagram, interpolate = FALSE))
figure1 <- pA / (wrap_elements(full = pB) | wrap_elements(full = pC)) /
  (wrap_elements(full = pD) | wrap_elements(full = pE)) /
  (wrap_elements(full = pF) | wrap_elements(full = pG)) +
  plot_layout(heights = c(0.95, 0.65, 0.85, 0.85))
save_panel(figure1, "figure1", 190, 165)
