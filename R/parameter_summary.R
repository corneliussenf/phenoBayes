# Get annual variation

parameter_summary <- function(x, type = "temporal") {
  if(type == "temporal") {
    pars <- as.data.frame(x)
    pars_extract <- pars[, names(pars) %in% grep("phi\\[", names(pars), value = TRUE)]
    summary <- as.data.frame(t(apply(pars_extract, 2, function(x) c(mean(x), sd(x), quantile(x, prob = c(0.025, 0.25, 0.5, 0.75, 0.975))))))
    names(summary) <- c("Mean", "Sd", "Q025", "Q25", "Q50", "Q75", "Q975")
    summary$year <- 1:nrow(summary)
  } else if(type == "spatial") {
    pars <- as.data.frame(x)
    summary <- lapply(1:5, function(i) {
      pars_extract <- pars[, names(pars) %in% grep(paste0("^beta\\[", i), names(pars), value = TRUE)]
      summary <- as.data.frame(t(apply(pars_extract, 2, function(x) c(mean(x), sd(x), quantile(x, prob = c(0.025, 0.25, 0.5, 0.75, 0.975))))))
      names(summary) <- c("Mean", "Sd", "Q025", "Q25", "Q50", "Q75", "Q975")
      summary$pixel <- 1:nrow(summary)
      summary$parameter <- c("minimum", "maximum", "changerate", "timing", "greendown")[i]
      return(summary)
    })
    summary <- do.call("rbind", summary)
  } else {
    stop("Type must be either 'temporal' or 'spatial'!")
  }
  return(summary)
}