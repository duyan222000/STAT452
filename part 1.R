# Packages
# install.packages(c("tidyverse","glmnet","corrplot","pROC","gridExtra"))
library(tidyverse)
library(glmnet)
library(corrplot)
library(pROC)
library(gridExtra)
library(grid)

set.seed(123)
dir.create("plots", showWarnings = FALSE)

# Unified PNG savers
save_png <- function(file, width=1600, height=1100, res=220, plot_expr){
  png(file, width=width, height=height, res=res)
  on.exit(dev.off(), add = TRUE)
  eval.parent(substitute(plot_expr))
}
png_table <- function(df, file, width=1200, height=600, res=220, rows=NULL){
  tbl <- gridExtra::tableGrob(df, rows = rows)
  png(file, width=width, height=height, res=res); grid::grid.draw(tbl); dev.off()
}
save_table_png <- function(df, file, caption=NULL, width=1200, height=600, res=220){
  if (!is.null(caption)) {
    df <- rbind(setNames(as.list(rep("", ncol(df))), names(df)), df)
    grob_cap <- textGrob(caption, gp=gpar(fontface="bold", cex=1.1))
    tbl <- tableGrob(df)
    lay <- rbind(c(1), c(2))
    g <- gtable:::gtable_add_rows(tbl, heights=grobHeight(grob_cap) + unit(4,"pt"), pos=0)
    g <- gtable:::gtable_add_grob(g, grob_cap, 1, 1, 1, ncol(tbl))
    png(file, width=width, height=height, res=res); grid::grid.draw(g); dev.off()
  } else {
    png_table(df, file, width, height, res)
  }
}

# Helpers
metrics_reg <- function(y, yhat){
  rmse <- sqrt(mean((y - yhat)^2))
  mae  <- mean(abs(y - yhat))
  r2   <- 1 - sum((y - yhat)^2)/sum((y - mean(y))^2)
  c(RMSE=rmse, MAE=mae, R2=r2)
}
trainTestSplit <- function(data, seed=0, trainRatio=.8){
  set.seed(seed); n <- nrow(data)
  tr <- sample(seq_len(n), size = round(trainRatio*n), replace = FALSE)
  list(train=data[tr,,drop=FALSE], test=data[-tr,,drop=FALSE])
}
stratified_split <- function(df, y="target", train_ratio=.8, seed=123){
  set.seed(seed)
  i1 <- which(df[[y]]==1); i0 <- which(df[[y]]==0)
  n1 <- round(train_ratio*length(i1)); n0 <- round(train_ratio*length(i0))
  tr_idx <- c(sample(i1, n1), sample(i0, n0))
  te_idx <- setdiff(seq_len(nrow(df)), tr_idx)
  list(train=df[tr_idx,,drop=FALSE], test=df[te_idx,,drop=FALSE])
}
iqr_trim_mask <- function(df_num){
  trims <- lapply(df_num, function(x){
    q1 <- quantile(x, .25, na.rm=TRUE); q3 <- quantile(x, .75, na.rm=TRUE)
    i  <- IQR(x, na.rm=TRUE); lo <- q1 - 1.5*i; hi <- q3 + 1.5*i
    x >= lo & x <= hi
  })
  Reduce(`&`, trims)
}

# Load & basic transforms
forest <- read.csv("data/ForestFire.csv", stringsAsFactors = FALSE)

forest <- forest %>%
  mutate(
    log_area = log1p(area),
    month = factor(month, levels = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")),
    day   = factor(day,   levels = c("mon","tue","wed","thu","fri","sat","sun")),
    target = as.integer(area > 0)
  )

# Missingness
na_counts <- colSums(is.na(forest))
save_table_png(as.data.frame(na_counts) |> rownames_to_column("Column"),
               "plots/table_missing_values.png", "Missing values per column")

# EDA
p_hist_raw <- ggplot(forest, aes(area)) +
  geom_histogram(bins=50, fill="#4C78A8") +
  labs(title="Distribution of Burned Area (raw)", x="Area (ha)", y="Count") +
  theme_minimal(base_size=13)
save_png("plots/fig_hist_area_raw.png", plot_expr = { print(p_hist_raw) })

p_hist_log <- ggplot(forest, aes(log_area)) +
  geom_histogram(bins=50, fill="#72B7B2") +
  labs(title="Distribution of Burned Area (log scale)", x="log(1 + Area)", y="Count") +
  theme_minimal(base_size=13)
save_png("plots/fig_hist_area_log.png", plot_expr = { print(p_hist_log) })

p_box_month <- ggplot(forest, aes(month, log_area)) +
  geom_boxplot(outlier.alpha=.45, fill="#E45756") +
  labs(title="Burned Area (log) by Month", x="Month", y="log(1+Area)") +
  theme_minimal(base_size=13)
save_png("plots/fig_box_month.png", plot_expr = { print(p_box_month) })

p_box_day <- ggplot(forest, aes(day, log_area)) +
  geom_boxplot(outlier.alpha=.45, fill="#F1A208") +
  labs(title="Burned Area (log) by Day of Week", x="Day of Week", y="log(1+Area)") +
  theme_minimal(base_size=13)
save_png("plots/fig_box_day.png", plot_expr = { print(p_box_day) })

# Correlation heatmap
num_df <- forest |> select(where(is.numeric))
corr_mat <- cor(num_df, method="spearman", use="pairwise.complete.obs")
save_png("plots/fig_corr_heatmap.png", width=2000, height=1600, plot_expr = {
  corrplot(corr_mat, type="upper", order="hclust",
           tl.cex=.7, tl.col="black", tl.srt=45,
           mar=c(0,0,2,0), main="Spearman Correlations")
})

# IQR outlier cleaning + before/after boxplots
num_cols <- names(select(forest, where(is.numeric)))
mask <- iqr_trim_mask(forest[num_cols])
forest_iqr <- forest[mask, , drop=FALSE]
removed_iqr <- nrow(forest) - nrow(forest_iqr)
save_table_png(data.frame(Removed_IQR = removed_iqr,
                          Kept = nrow(forest_iqr),
                          Total = nrow(forest)),
               "plots/table_iqr_removed.png", "IQR trim — row counts")

# Drop 'rain' if constant post-clean
if ("rain" %in% names(forest_iqr) && length(unique(forest_iqr$rain)) <= 1L) {
  forest_iqr$rain <- NULL
}

p_before <- forest %>% select(where(is.numeric)) %>% pivot_longer(everything()) %>%
  ggplot(aes(name, value)) + geom_boxplot(outlier.alpha=.25) + coord_flip() +
  labs(title="Numeric variables BEFORE IQR cleaning", x=NULL, y=NULL) + theme_minimal(12)
save_png("plots/fig_box_before_iqr.png", width=1400, height=1200, plot_expr = { print(p_before) })

# Cook's distance based on a sensible baseline model
rhs <- intersect(c("X","Y","month","day","FFMC","DMC","DC","ISI","temp","RH","wind","rain"),
                 names(forest_iqr))
fml <- reformulate(rhs, response = "log_area")
lm0 <- lm(fml, data = forest_iqr)
cooks <- cooks.distance(lm0); cooks[!is.finite(cooks)] <- NA_real_
thr <- 4/length(cooks)
save_png("plots/fig_cooks.png", width=1600, height=900, plot_expr = {
  ymax <- max(cooks, na.rm=TRUE); if(!is.finite(ymax)) ymax <- 1
  plot(cooks, type="h", ylim=c(0, 1.08*ymax),
       main="Cook's Distance for Influential Observations",
       ylab="Cook's Distance", xlab="Index")
  abline(h=thr, col="red", lty=2)
})
infl <- which(is.finite(cooks) & cooks > thr)
forest_final <- if (length(infl)) forest_iqr[-infl,,drop=FALSE] else forest_iqr
save_table_png(data.frame(Cooks_Removed = length(infl)), "plots/table_cooks_removed.png",
               "Cook's distance — removed rows")

p_after <- forest_final %>% select(where(is.numeric)) %>% pivot_longer(everything()) %>%
  ggplot(aes(name, value)) + geom_boxplot(outlier.alpha=.25) + coord_flip() +
  labs(title="Numeric variables AFTER cleaning", x=NULL, y=NULL) + theme_minimal(12)
save_png("plots/fig_box_after_iqr_and_cook's_distance.png", width=1400, height=1200, plot_expr = { print(p_after) })

# =========================
# BOX–COX STEP (custom, no MASS)
# =========================

# Ensure response is strictly positive
eps <- 1e-6
y_pos <- forest_final$log_area + eps

# Model matrix for predictors (rhs was defined earlier)
f_rhs <- reformulate(rhs)
Xmat  <- model.matrix(f_rhs, data = forest_final)

# BOX–COX


# Ensure strictly-positive response for Box–Cox
eps <- 1e-6
y_raw <- forest_final$log_area
y_pos <- y_raw + eps
if (any(y_pos <= 0 | !is.finite(y_pos))) {
  shift <- abs(min(y_pos[is.finite(y_pos)], na.rm = TRUE)) + 1e-3
  y_pos <- y_pos + shift
  message(sprintf("Box–Cox: shifted response by +%.6g to ensure positivity.", shift))
}

# Build the predictor set
rhs <- intersect(
  c("X","Y","month","day","FFMC","DMC","DC","ISI","temp","RH","wind","rain"),
  names(forest_final)
)

# Drop predictors with zero variance
rhs <- rhs[vapply(forest_final[rhs], function(col) length(unique(col)) > 1L, logical(1))]

# Model matrices for fast fitting
fml_x <- reformulate(rhs)           
X <- model.matrix(fml_x, data = forest_final)
n <- NROW(X)

# Box–Cox profile log-likelihood:
bc_loglik <- function(lambda, y, X) {
  # Transform response
  if (abs(lambda) < 1e-12) {
    yt <- log(y)
  } else {
    yt <- (y^lambda - 1) / lambda
  }
  # OLS fit via lm.fit for speed/stability
  fit <- lm.fit(X, yt)
  rss <- sum(fit$residuals^2)
  (lambda - 1) * sum(log(y)) - 0.5 * n * log(rss / n)
}

# Coarse grid search
lambda_grid <- seq(-2, 2, by = 0.05)
ll_vals <- vapply(lambda_grid, bc_loglik, numeric(1), y = y_pos, X = X)
lambda_coarse <- lambda_grid[which.max(ll_vals)]

# Refine with 1D optimize() around the best coarse point
ref_lo <- max(-2, lambda_coarse - 0.25)
ref_hi <- min( 2, lambda_coarse + 0.25)
opt <- optimize(function(l) -bc_loglik(l, y_pos, X), interval = c(ref_lo, ref_hi))
lambda_opt <- opt$minimum

# Save the profile plot with chosen lambda marked
save_png("plots/fig_boxcox_profile.png", width = 1200, height = 900, plot_expr = {
  plot(lambda_grid, ll_vals, type = "l", lwd = 2,
       xlab = expression(lambda), ylab = "Profile log-likelihood (up to const)",
       main = sprintf("Box–Cox Profile (λ* ≈ %.3f)", lambda_opt))
  abline(v = lambda_opt, col = "red", lty = 2, lwd = 2)
})

# Apply the Box–Cox transform using lambda_opt
if (abs(lambda_opt) < 1e-12) {
  y_bc <- log(y_pos)
} else {
  y_bc <- (y_pos^lambda_opt - 1) / lambda_opt
}

# Fit baseline (original log-scale + eps) vs Box–Cox model
forest_final$log_area_pos <- y_pos  # keep for reproducibility
fml_log <- reformulate(rhs, response = "log_area_pos")
lm_log  <- lm(fml_log, data = forest_final)

forest_final$y_boxcox <- y_bc
fml_bc <- reformulate(rhs, response = "y_boxcox")
lm_bc  <- lm(fml_bc, data = forest_final)

# Compare in-sample fit (R^2 and overall F-test p-value)
sum_log <- summary(lm_log); sum_bc <- summary(lm_bc)
f_p_log <- pf(sum_log$fstatistic[1], sum_log$fstatistic[2], sum_log$fstatistic[3], lower.tail = FALSE)
f_p_bc  <- pf(sum_bc$fstatistic[1],  sum_bc$fstatistic[2],  sum_bc$fstatistic[3],  lower.tail = FALSE)

cmp_tbl <- data.frame(
  Metric            = c("Multiple R-squared", "F-statistic p-value"),
  `Original (log)`  = c(round(sum_log$r.squared, 4), signif(f_p_log, 4)),
  `Box–Cox`         = c(round(sum_bc$r.squared, 4),  signif(f_p_bc,  4))
)

save_table_png(cmp_tbl, "plots/table_boxcox_compare.png",
               caption = "Model performance before and after Box–Cox transform")

# Residual diagnostics for the Box–Cox model
save_png("plots/fig_residuals_boxcox.png", width = 1800, height = 1200, plot_expr = {
  op <- par(mfrow = c(2,2))
  plot(lm_bc)
  par(op)
})

# REGRESSION (glmnet)
split_r <- trainTestSplit(forest_final, seed = 0, trainRatio = .8)
tr_r <- split_r$train; te_r <- split_r$test

mm_reg <- model.matrix(~ . - area - log_area - target, data = rbind(tr_r, te_r))
ntr <- nrow(tr_r)
x_tr <- mm_reg[1:ntr, -1, drop=FALSE]
x_te <- mm_reg[(ntr+1):nrow(mm_reg), -1, drop=FALSE]
y_tr <- tr_r$log_area
y_te <- te_r$log_area

set.seed(123)
cv.ridge <- cv.glmnet(x_tr, y_tr, alpha=0,   standardize=TRUE)
cv.lasso <- cv.glmnet(x_tr, y_tr, alpha=1,   standardize=TRUE)
alpha_grid <- seq(0.1, 0.9, by=.2)
cv_list <- lapply(alpha_grid, function(a) cv.glmnet(x_tr, y_tr, alpha=a, standardize=TRUE))
best_en_idx <- which.min(sapply(cv_list, function(cv) min(cv$cvm)))
best_alpha <- alpha_grid[best_en_idx]; cv.elastic <- cv_list[[best_en_idx]]

ridge.model   <- glmnet(x_tr, y_tr, alpha=0,          lambda=cv.ridge$lambda.min,   standardize=TRUE)
lasso.model   <- glmnet(x_tr, y_tr, alpha=1,          lambda=cv.lasso$lambda.min,   standardize=TRUE)
elastic.model <- glmnet(x_tr, y_tr, alpha=best_alpha, lambda=cv.elastic$lambda.min, standardize=TRUE)

save_png("plots/fig_cv_ridge.png", 1400, 1000, plot_expr = { plot(cv.ridge);  title("Ridge CV",   line=2.5) })
save_png("plots/fig_cv_lasso.png", 1400, 1000, plot_expr = { plot(cv.lasso);  title("LASSO CV",   line=2.5) })
save_png("plots/fig_cv_elnet.png", 1400, 1000, plot_expr = { plot(cv.elastic); title(paste0("Elastic Net CV (alpha=",best_alpha,")"), line=2.5) })

pred_ridge   <- as.numeric(predict(ridge.model,   newx=x_te))
pred_lasso   <- as.numeric(predict(lasso.model,   newx=x_te))
pred_elastic <- as.numeric(predict(elastic.model, newx=x_te))

res_tbl <- rbind(
  LASSO   = metrics_reg(y_te, pred_lasso),
  RIDGE   = metrics_reg(y_te, pred_ridge),
  ELASTIC = metrics_reg(y_te, pred_elastic)
) %>% as.data.frame() %>% rownames_to_column("Model") %>% mutate(across(-Model, ~round(.,4)))
save_table_png(res_tbl, "plots/table_regression_metrics.png", "Predictive performance (test set)")

# Residual diagnostics for LASSO
lasso_resid <- y_te - pred_lasso
ok <- is.finite(lasso_resid) & is.finite(pred_lasso)
save_png("plots/fig_residuals_lasso.png", width=1800, height=600, plot_expr = {
  par(mfrow=c(1,3))
  if (any(ok)) {
    hist(lasso_resid[ok], breaks=30, main="LASSO Residuals (Test)", xlab="Residual")
    plot(pred_lasso[ok], lasso_resid[ok], xlab="Predicted log(1+Area)", ylab="Residual",
         main="Residuals vs Predicted"); abline(h=0, col="red", lty=2)
    qqnorm(lasso_resid[ok], main="Q–Q Plot of Residuals"); qqline(lasso_resid[ok], col="red")
  } else {
    plot.new(); title("No finite residuals"); plot.new(); title("No finite residuals"); plot.new(); title("No finite residuals")
  }
  par(mfrow=c(1,1))
})

# CLASSIFICATION (glmnet)
forest_final$target <- as.integer(forest_final$area > 0)
split_c <- stratified_split(forest_final, y="target", train_ratio=.8, seed=1111)
tr_c <- split_c$train; 
te_c <- split_c$test

mm_clf <- model.matrix(~ . - area - log_area - target, data = rbind(tr_c, te_c))
ntrc <- nrow(tr_c)
x_tr_c <- mm_clf[1:ntrc, -1, drop=FALSE]
x_te_c <- mm_clf[(ntrc+1):nrow(mm_clf), -1, drop=FALSE]
y_tr_c <- tr_c$target; y_te_c <- te_c$target

# class weights (prioritize recall)
weight_ratio <- sum(y_tr_c == 0) / sum(y_tr_c == 1)
w <- ifelse(y_tr_c == 1, weight_ratio, 1)

set.seed(1)
cv.logit <- cv.glmnet(x_tr_c, y_tr_c, family="binomial", weights=w, standardize=TRUE)
logit    <- glmnet(x_tr_c, y_tr_c, family="binomial", weights=w, lambda=cv.logit$lambda.min, standardize=TRUE)

prob <- as.numeric(predict(logit, newx=x_te_c, type="response"))

# Recall vs threshold
th <- seq(0,1,by=.01)
rec <- sapply(th, function(t){
  lab <- as.integer(prob >= t)
  TP <- sum(lab==1 & y_te_c==1); FN <- sum(lab==0 & y_te_c==1)
  ifelse(TP+FN==0, NA, TP/(TP+FN))
})
best_rec_th <- th[ which.max(replace(rec, is.na(rec), -Inf)) ]

df_rec <- tibble(threshold = th, recall = rec)
p_rec <- ggplot(df_rec, aes(threshold, recall)) +
  geom_line(linewidth=1) +
  geom_vline(xintercept = best_rec_th, linetype=2, colour="red") +
  labs(title="Recall vs Threshold (Fire Occurrence)", x="Threshold", y="Recall") +
  theme_minimal(base_size=13)
save_png("plots/fig_recall_threshold.png", 1200, 800, plot_expr = { print(p_rec) })

# ROC + AUC
roc_obj <- pROC::roc(response=y_te_c, predictor=prob, na.rm=TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))
save_png("plots/fig_roc.png", width=900, height=800, plot_expr = {
  plot(roc_obj, col="blue", lwd=2, main=paste0("ROC (AUC = ", round(auc_val,3), ")"))
  abline(a=0, b=1, lty=2, col="red")
})

# Confusion matrices + metrics at t=0.50, Max Recall, Best F1
f1_at <- function(r,p){ ifelse(is.na(p) || (p+r)==0, NA, 2*p*r/(p+r)) }
cm_metrics <- function(prob, y_true, thr){
  pred <- as.integer(prob >= thr)
  TP <- sum(pred==1 & y_true==1); FP <- sum(pred==1 & y_true==0)
  FN <- sum(pred==0 & y_true==1); TN <- sum(pred==0 & y_true==0)
  Recall <- TP/(TP+FN)
  Precision <- ifelse((TP+FP)==0, NA, TP/(TP+FP))
  F1 <- f1_at(Recall, Precision)
  list(TP=TP, FP=FP, FN=FN, TN=TN, Recall=Recall, Precision=Precision, F1=F1)
}

# Best-F1 threshold search
pre <- sapply(th, function(t){
  lab <- as.integer(prob >= t); TP <- sum(lab==1 & y_te_c==1); FP <- sum(lab==1 & y_te_c==0)
  ifelse((TP+FP)==0, NA, TP/(TP+FP))
})
f1 <- mapply(f1_at, rec, pre)
best_f1_th <- th[ which.max(replace(f1, is.na(f1), -Inf)) ]

m05 <- cm_metrics(prob, y_te_c, 0.50)
mBR <- cm_metrics(prob, y_te_c, best_rec_th)
mBF <- cm_metrics(prob, y_te_c, best_f1_th)
m04 <- cm_metrics(prob, y_te_c, 0.45)

metrics_df <- tibble(
  Setting   = c("t=0.50 (default)", "t=0.45", sprintf("t=%.2f (Max Recall)",best_rec_th), sprintf("t=%.2f (Best F1)",best_f1_th)),
  Threshold = c(0.50, 0.45, best_rec_th, best_f1_th),
  Recall    = round(c(m05$Recall, m04$Recall, mBR$Recall, mBF$Recall), 3),
  Precision = round(c(m05$Precision, m04$Precision, mBR$Precision, mBF$Precision), 3),
  F1        = round(c(m05$F1, m04$F1, mBR$F1, mBF$F1), 3),
  AUC       = round(auc_val, 3)
)
save_table_png(metrics_df, "plots/table_classification_metrics.png", "Classification metrics (test set)")

cm_to_df <- function(stat){
  data.frame(`Pred \\ True`=c("0","1"),
             `0`=c(stat$TN, stat$FP),
             `1`=c(stat$FN, stat$TP))
}
save_table_png(cm_to_df(m05), "plots/table_confusion_matrix_t05.png", caption="Confusion Matrix (t=0.50)")
save_table_png(cm_to_df(m04), "plots/table_confusion_matrix_t045.png", caption="Confusion Matrix (t=0.45)")
save_table_png(cm_to_df(mBR), "plots/table_confusion_matrix_maxRecall.png",
               caption=sprintf("Confusion Matrix (t=%.2f, Max Recall)", best_rec_th))
save_table_png(cm_to_df(mBF), "plots/table_confusion_matrix_bestF1.png",
               caption=sprintf("Confusion Matrix (t=%.2f, Best F1)", best_f1_th))

message(
  "Saved to 'plots/':\n",
  "- fig_hist_area_raw.png, fig_hist_area_log.png\n",
  "- fig_box_month.png, fig_box_day.png\n",
  "- fig_corr_heatmap.png\n",
  "- table_missing_values.png, table_iqr_removed.png, table_cooks_removed.png\n",
  "- fig_box_before_iqr.png, fig_box_after_iqr.png\n",
  "- fig_cooks.png\n",
  "- fig_cv_ridge.png, fig_cv_lasso.png, fig_cv_elnet.png\n",
  "- table_regression_metrics.png, fig_residuals_lasso.png\n",
  "- fig_recall_threshold.png, fig_roc.png\n",
  "- table_classification_metrics.png\n",
  "- table_confusion_matrix_t05.png, table_confusion_matrix_maxRecall.png, table_confusion_matrix_bestF1.png\n"
)