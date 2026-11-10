# Figure 6. Recompute ROC, external validation intervals and Brier scores
# from the frozen patient-level predictions of the final model workflow.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "common.R"))
features <- c("Microbiome", "Metabolites", "Clinical")
model_levels <- c("RF", "XGB", "GBM", "Lasso", "Logistic Regression", "QDA", "SVM")
balance_levels <- c("Unbalanced", "Balanced")
d <- read_input("model_predictions.csv",
  c("FeatureSet", "Cohort", "PatientID", "Timepoint", "Outcome", "Model", "Balance", "Probability"))
spec <- read_input("validation_model_manifest.csv", c("FeatureSet", "Model", "Balance", "Threshold"))
stopifnot(all(d$FeatureSet %in% features), all(d$Model %in% model_levels),
          all(d$Balance %in% balance_levels), all(d$Cohort %in% c("DiscoveryCV", "Validation")),
          all(d$Timepoint == "T0"), all(d$Outcome %in% 0:1),
          all(is.finite(d$Probability) & d$Probability >= 0 & d$Probability <= 1),
          !anyNA(d$PatientID),
          !anyDuplicated(d[c("FeatureSet", "Cohort", "Model", "Balance", "PatientID")]),
          nrow(spec) == 3, setequal(spec$FeatureSet, features),
          all(spec$Model[match(features, spec$FeatureSet)] ==
                c("QDA", "Logistic Regression", "Logistic Regression")),
          all(spec$Balance %in% balance_levels),
          all(is.finite(spec$Threshold) & spec$Threshold >= 0 & spec$Threshold <= 1))
dc <- filter(d, Cohort == "DiscoveryCV")
stopifnot(nrow(distinct(dc, FeatureSet, Model, Balance)) ==
            length(features) * length(model_levels) * length(balance_levels),
          all(vapply(split(dc, dc$FeatureSet), function(x) {
            all(table(x$PatientID) == length(model_levels) * length(balance_levels))
          }, logical(1))),
          all((d |> distinct(Cohort, PatientID, Outcome) |> count(Cohort, PatientID))$n == 1))
validation <- filter(d, Cohort == "Validation") |>
  inner_join(spec, by = c("FeatureSet", "Model", "Balance"))
stopifnot(setequal(validation$FeatureSet, features))
fig_data <- bind_rows(dc, validation)
roc_data <- list()
summary_data <- list()
set.seed(42)
blocks <- split(fig_data, interaction(fig_data$FeatureSet, fig_data$Cohort, fig_data$Model,
                                      fig_data$Balance, drop = TRUE))
roc_object <- function(x) {
  pROC::roc(x$Outcome, x$Probability, levels = c(0, 1), direction = "<",
            quiet = TRUE, auc = FALSE)
}
fast_auc <- function(y, p) {
  n1 <- sum(y == 1)
  n0 <- length(y) - n1
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
for (key in names(blocks)) {
  x <- blocks[[key]]
  stopifnot(length(unique(x$Outcome)) == 2)
  a <- roc_object(x)
  labels <- x[1, c("FeatureSet", "Cohort", "Model", "Balance")]
  roc_data[[key]] <- cbind(labels[rep(1, length(a$specificities)), ],
                         FPR = 1 - a$specificities, TPR = a$sensitivities)
  outcome <- x$Outcome
  score <- x$Probability
  z <- c(AUROC = fast_auc(outcome, score), Specificity = NA_real_, Accuracy = NA_real_,
         Brier = mean((score - outcome)^2))
  lower <- upper <- setNames(rep(NA_real_, length(z)), names(z))
  if (labels$Cohort == "Validation") {
    positive <- which(outcome == 1)
    negative <- which(outcome == 0)
    predicted <- as.integer(score >= unique(x$Threshold))
    z["Specificity"] <- mean(predicted[negative] == 0)
    z["Accuracy"] <- mean(predicted == outcome)
    sim <- replicate(5000, {
      i0 <- sample(negative, length(negative), replace = TRUE)
      i1 <- sample(positive, length(positive), replace = TRUE)
      i <- c(i0, i1)
      c(fast_auc(outcome[i], score[i]), mean(predicted[i0] == 0),
        mean(predicted[i] == outcome[i]), mean((score[i] - outcome[i])^2))
    })
    lower[] <- apply(sim, 1, quantile, .025)
    upper[] <- apply(sim, 1, quantile, .975)
  }
  summary_data[[key]] <- cbind(labels, N = nrow(x), Cases = sum(outcome), AUROC = z["AUROC"],
    AUROC_Lower = lower["AUROC"], AUROC_Upper = upper["AUROC"], Specificity = z["Specificity"],
    Specificity_Lower = lower["Specificity"], Specificity_Upper = upper["Specificity"],
    Accuracy = z["Accuracy"], Accuracy_Lower = lower["Accuracy"], Accuracy_Upper = upper["Accuracy"],
    Brier = z["Brier"], Brier_Lower = lower["Brier"], Brier_Upper = upper["Brier"])
}
curves <- bind_rows(roc_data)
metrics <- bind_rows(summary_data)
write_output(curves, "figure06_roc_coordinates.csv")
write_output(metrics, "figure06_model_performance.csv")
ranked <- filter(metrics, Cohort == "DiscoveryCV") |> arrange(FeatureSet, desc(AUROC))
write_output(ranked, "figure06_discovery_model_ranking.csv")
best <- ranked |> group_by(FeatureSet) |> filter(AUROC == max(AUROC)) |> ungroup()
checked <- spec |> mutate(IsDiscoveryMaximum = paste(FeatureSet, Model, Balance) %in%
                            paste(best$FeatureSet, best$Model, best$Balance))
write_output(checked, "figure06_validation_selection_check.csv")
if (any(!checked$IsDiscoveryMaximum)) {
  warning("At least one frozen validation model differs from the top plotted discovery AUROC.",
          " See figure06_validation_selection_check.csv.")
}
model_colors <- setNames(c("#E64943", "#F2AC32", "#944CA1", "#62A74A",
                          "#267E78", "#4A6CB4", "#666666"), model_levels)
plot_roc <- function(k, cohort) {
  a <- filter(curves, FeatureSet == k, Cohort == cohort)
  b <- filter(metrics, FeatureSet == k, Cohort == cohort)
  if (cohort == "Validation") {
    note <- sprintf("%s, AUROC %.3f (95%% CI %.3f–%.3f)", b$Model, b$AUROC,
                    b$AUROC_Lower, b$AUROC_Upper)
  } else note <- sprintf("Best discovery AUROC %.3f", max(b$AUROC))
  ggplot(a, aes(FPR, TPR, colour = Model, linetype = Balance, group = interaction(Model, Balance))) +
    geom_abline(intercept = 0, slope = 1, linetype = 3, colour = "grey65") +
    geom_path(linewidth = .8) +
    scale_colour_manual(values = model_colors) +
    scale_linetype_manual(values = c(Unbalanced = "solid", Balanced = "dashed")) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + figure_theme() +
    labs(title = paste(k, if (cohort == "DiscoveryCV") "discovery" else "validation"),
         subtitle = note, x = "1 − specificity", y = "Sensitivity")
}
panels <- unlist(lapply(features, function(k) list(plot_roc(k, "DiscoveryCV"), plot_roc(k, "Validation"))),
                 recursive = FALSE)
save_figure(wrap_plots(panels, ncol = 2) + plot_annotation(tag_levels = "A"), "figure06", 13, 16)
