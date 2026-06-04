set.seed(123)

obs <- read.csv("nsw_observational.csv")
ctrl <- obs[obs$treat == 0, ]
trt <- obs[obs$treat == 1, ]

cat("Original treated rows:", nrow(trt), "\n")
cat("Original control rows:", nrow(ctrl), "\n\n")

## ---- Step 1: Bootstrap-resample treated rows ----
n_trt <- nrow(trt)
boot_idx <- sample(seq_len(n_trt), size = n_trt,
                   replace = TRUE)
trt_boot <- trt[boot_idx, ]

continuous <- c("age", "educ", "re74", "re75", "re78")
for (v in continuous) {
  noise <- rnorm(n_trt, mean = 0, sd = sd(trt[[v]]) * 0.05)
  trt_boot[[v]] <- trt_boot[[v]] + noise
}

trt_boot$age <- pmax(round(trt_boot$age), 16)
trt_boot$educ <- pmax(round(trt_boot$educ), 0)
trt_boot$re74 <- pmax(trt_boot$re74, 0)
trt_boot$re75 <- pmax(trt_boot$re75, 0)
trt_boot$re78 <- pmax(trt_boot$re78, 0)

rownames(trt_boot) <- NULL

resampled <- rbind(trt_boot, ctrl)
resampled <- resampled[order(
  -resampled$treat, runif(nrow(resampled))
), ]
rownames(resampled) <- NULL

write.csv(resampled, "nsw_observational_resample.csv",
          row.names = FALSE)
cat("Wrote nsw_observational_resample.csv:",
    nrow(resampled), "rows\n\n")

## ---- Step 2: Detect bimodal valley ----
cat("=== Distribution diagnostics ===\n")
for (v in c("re74", "re75", "re78")) {
  vals <- resampled[[v]]
  cat(
    v, "| zero:", sum(vals == 0),
    "| positive:", sum(vals > 0),
    "| mean:", round(mean(vals)),
    "| median:", round(median(vals)),
    "| Q75:", round(quantile(vals, 0.75)),
    "| Q90:", round(quantile(vals, 0.90)), "\n"
  )
}

find_valley <- function(x, var_name) {
  d <- density(x, bw = "SJ", n = 1024)
  peaks <- which(diff(sign(diff(d$y))) == -2) + 1
  if (length(peaks) >= 2) {
    p1 <- peaks[1]
    p2 <- peaks[length(peaks)]
    region <- p1:p2
    valley_idx <- region[which.min(d$y[region])]
    cutoff <- d$x[valley_idx]
    cat(
      var_name, "| peaks at",
      round(d$x[p1]), "and", round(d$x[p2]),
      "| valley cutoff:", round(cutoff), "\n"
    )
    return(cutoff)
  }
  med <- median(x)
  cat(var_name, ": no clear bimodality, using median",
      round(med), "\n")
  return(med)
}

cat("\n=== Valley detection on all rows ===\n")
cut_74 <- find_valley(resampled$re74, "re74")
cut_75 <- find_valley(resampled$re75, "re75")
cut_78 <- find_valley(resampled$re78, "re78")

## ---- Step 3: Trim upper peak ----
in_upper <- (resampled$re74 > cut_74) &
  (resampled$re75 > cut_75) &
  (resampled$re78 > cut_78)

cat(
  "\nRows in upper peak (all three above cutoff):",
  sum(in_upper), "of", nrow(resampled), "\n"
)

in_upper_any <- (resampled$re74 > cut_74) |
  (resampled$re75 > cut_75) |
  (resampled$re78 > cut_78)

cat(
  "Rows in upper peak (any one above cutoff):",
  sum(in_upper_any), "of", nrow(resampled), "\n"
)

trimmed <- resampled[!in_upper_any, ]
rownames(trimmed) <- NULL

write.csv(
  trimmed,
  "nsw_observational_resampling_trimmed.csv",
  row.names = FALSE
)
cat(
  "\nWrote nsw_observational_resampling_trimmed.csv:",
  nrow(trimmed), "rows\n"
)

cat("\n=== Treat counts after trimming ===\n")
print(table(trimmed$treat))

cat("\n=== Trimmed distribution summary ===\n")
for (v in c("re74", "re75", "re78")) {
  cat(v, "| mean:", round(mean(trimmed[[v]])),
      "| max:", round(max(trimmed[[v]])), "\n")
}
