# Run from the repository root: Rscript script/render_all.R
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1],mustWork=TRUE) else normalizePath(".",mustWork=TRUE)
if(!requireNamespace("rmarkdown",quietly=TRUE)) stop("Install rmarkdown before rendering.")
if(!rmarkdown::pandoc_available()) stop("Pandoc is required. Use RStudio or configure RSTUDIO_PANDOC.")
for(i in 1:6) {
  input <- file.path(root,paste0("fig",i),paste0("fig",i,".Rmd"))
  ggplot2::theme_set(ggplot2::theme_gray())
  message("Rendering Figure ",i)
  rmarkdown::render(input,envir=new.env(parent=globalenv()),encoding="UTF-8",quiet=FALSE)
}
