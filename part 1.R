# install.packages("tidyverse")
# install.packages("glmnet")
# install.packages("corrplot")

library(tidyverse)
library(glmnet)
library(corrplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pROC)


forest <- read.csv("data/ForestFire.csv")

# Relevant data visualizations of features of this dataset
str(forest)
summary(forest)

# Check distribution of burned area
hist(forest$area, breaks = 50, main = "Distribution of Burned Area", xlab = "Area")

# We can see that most are close to zero but some are huge surpassing 1000 area
# use log_area instead cuz most forest fire are small but some are huge which make it skew
# toward the right
forest$log_area <- log1p(forest$area)
model <- lm(log_area ~ ., data = forest %>% select(-area))

hist(forest$log_area,
     breaks = 50,
     main = "Distribution of Burned Area (Log Transformed)",
     xlab = "log(1 + Area)")

# Cleanup dataset. Are there any outliers ? Are there any missing values in any of the features ?
# Explain how you handle categorial features in the dataset.
# Check for missing values
colSums(is.na(forest))

# remove outliers with cook distance
# Calculate cook distance
cooksd <- cooks.distance(model)

cooksd[!is.finite(cooksd)] <- 0

plot(
  cooksd, type = "h",
  ylim = c(0, max(cooksd) * 1.1),
  main = "Cook's Distance for Influential Observations",
  ylab = "Cook's Distance"
)
abline(h = 4 / length(cooksd), col = "red", lty = 2)

# Identify influential points
threshold <- 4 / length(cooksd)
influential <- as.numeric(names(cooksd)[(cooksd > threshold)])
influential <- influential[!is.na(influential)]
influential
# Show how many are removed
cat("Removed:", length(influential), " (",round(100*length(influential)/nrow(forest), 1), "%)\n", sep="")

# Remove influential observations
forest_clean <- forest[-influential, ]

# Handle categorial features
# https://www.r-bloggers.com/2022/01/handling-categorical-data-in-r-part-1/
# We convert it to factor so when needed, R can create dummy tables
forest_clean$month <- factor(forest_clean$month,
                       levels = c("jan", "feb", "mar", "apr", "may", "jun",
                                  "jul", "aug", "sep", "oct", "nov", "dec"))

forest_clean$day <- factor(forest_clean$day,
                     levels = c("mon", "tue", "wed", "thu", "fri", "sat", "sun"))

str(forest_clean)

num_df <- forest_clean |> select(where(is.numeric))
corr_mat <- cor(num_df, method = "spearman", use = "pairwise.complete.obs")
# Correlation heatmap
corrplot(
  corr_mat,
  type = "upper",
  order = "hclust",
  tl.cex = 0.6,            # smaller text size
  tl.col = "black",        # text color
  tl.srt = 45,             # rotate labels 45 degrees
  main = "Correlation of Fire-Related Variables (Montesinho Park)",
  mar = c(0,0,2,0)         # adjust plot margins
)

# Boxplot of log_area by month
ggplot(forest_clean, aes(x = month, y = log_area)) + geom_boxplot(outlier.alpha = 0.4) +
  labs(title = "Burned Area (log scale) by Month in Montesinho Park",
       x = "Month", y = "log(1+Area)") +
  theme_minimal(base_size = 12)

# Boxplot of log_area by day of week
ggplot(forest_clean, aes(x = day, y = log_area)) +
  geom_boxplot(outlier.alpha = 0.4) +
  labs(title = "Burned Area (log scale) by Day of Week in Montesinho Park",
  x = "Day of Week", y = "log(1+Area)") +
  theme_minimal(base_size = 12)

# Boxplots for numeric predictors
forest_clean |>
  select(where(is.numeric), -area) |>
  pivot_longer(cols = everything()) |>
  ggplot(aes(x = name, y = value)) +
  geom_boxplot(outlier.alpha = 0.3) +
  coord_flip() +
  labs(title = "Distribution of Fire-Related Predictors (Montesinho Park)",
       x = "Variables", y = "Value") +
  theme_minimal(base_size = 12)

# Apply appropriate transformations of the features and/or output.
trainTestSplit = function(data, seed, trainRatio = 0.8){
  set.seed(seed)
  dataSize = nrow(data)
  trainSize = round(trainRatio * dataSize)
  trainIndex = sample(1:dataSize, replace = FALSE, size = trainSize)
  
  trainData = data[trainIndex,]
  testData = data[-trainIndex,]
  return(list(train = trainData, test = testData))
  }
fire.split = trainTestSplit(forest_clean, 0, trainRatio = 0.8)
fire.train = fire.split$train
fire.test = fire.split$test

mm_formula <- as.formula("~ . - area")  # exclude raw area (use log_area instead)
mm_all <- model.matrix(mm_formula, data = rbind(fire.train, fire.test))

n_tr <- nrow(fire.train)
x_all   <- mm_all[, -1, drop = FALSE]      # drop intercept column
x_train <- x_all[1:n_tr, , drop = FALSE]
x_test  <- x_all[(n_tr + 1):nrow(x_all), , drop = FALSE]

y_train <- fire.train$log_area
y_test  <- fire.test$log_area
# Use appropriate feature selection methods. Explain how you setup Ridge/Lasso/Elastic net 
# regularization method and interpret the result.

set.seed(123)

# Ridge Regression=
cv.ridge <- cv.glmnet(x_train, y_train, alpha = 0)
ridge.model <- glmnet(x_train, y_train, lambda = cv.ridge$lambda.min, alpha = 0)
plot(cv.ridge);   title("Ridge CV",   line = 2.5)

# Lasso Regression
cv.lasso <- cv.glmnet(x_train, y_train, alpha = 1, standardize = TRUE)
lasso.model <- glmnet(x_train, y_train, lambda = cv.lasso$lambda.min, alpha = 1)
plot(cv.lasso);   title("LASSO CV", line = 2.5)

nz <- as.matrix(coef(lasso.model))
nz <- nz[nz[,1] != 0, , drop = FALSE] # keep non-zero
nz <- nz[order(abs(nz[,1]), decreasing = TRUE), , drop = FALSE] # sort by magnitude
print(nz)

# Elastic Net
cv.elastic <- cv.glmnet(x_train, y_train, alpha = 0.5)
elastic.model <- glmnet(x_train, y_train, lambda = cv.elastic$lambda.min, alpha = 0.5)
plot(cv.elastic); title("Elastic Net CV", line = 2.5)

# Coefficients from best models
coef(ridge.model)
coef(lasso.model)
# Certain months like September, december have high influence on the burned area
# The DMC value and temp are also picked through the lasso model
coef(elastic.model)
# Elastic model select the same predictor as lasso meaning that they are important across models
# These results suggest that fire seasonality (December and September, etc) 
# and dryness indicators (DMC)
# play key roles in explaining fire area variation

# Cross validation
# Prepare test data

# Predictions
y_pred <- predict(lasso.model, s = cv.lasso$lambda.min, newx = x_test)

pred_ridge  <- predict(ridge.model,  newx = x_test, s = cv.ridge$lambda.min)
pred_elastic<- predict(elastic.model,newx = x_test, s = cv.elastic$lambda.min)

metrics <- function(y, yhat){
  rmse <- sqrt(mean((y - yhat)^2))
  mae  <- mean(abs(y - yhat))
  r2   <- 1 - sum((y - yhat)^2)/sum((y - mean(y))^2)
  c(RMSE=rmse, MAE=mae, R2=r2)
}
res <- rbind(
  LASSO   = metrics(y_test, as.numeric(y_pred)),
  RIDGE   = metrics(y_test, as.numeric(pred_ridge)),
  ELASTIC = metrics(y_test, as.numeric(pred_elastic))
)
print(round(res, 4))

# --- Residual diagnostics for LASSO on the test set ---
lasso_resid <- as.numeric(y_test - y_pred)
par(mfrow = c(1,3))

hist(lasso_resid, breaks = 30, main = "LASSO Residuals (Test)", xlab = "Residual")

plot(as.numeric(y_pred), lasso_resid,
     xlab = "Predicted log(1+Area)", ylab = "Residual",
     main = "Residuals vs Predicted")
abline(h = 0, col = "red", lty = 2)

qqnorm(lasso_resid, main = "Q–Q Plot of Residuals")
qqline(lasso_resid, col = "red")

par(mfrow = c(1,1))

# Evaluation
rmse <- sqrt(mean((y_test - y_pred)^2))
mae <- mean(abs(y_test - y_pred))
r2 <- 1 - sum((y_test - y_pred)^2) / sum((y_test - mean(y_test))^2)

cat("RMSE:", rmse, "\n")
cat("MAE:", mae, "\n")
cat("R²:", r2, "\n")

# Classification task
# We should use Recall score as classification metric because:
# we want to focus on avoiding false negative cases when we predict that 
# there is no fire but a fire occur
# === CLASSIFICATION (fixed: exclude 'target' from predictors) ===

# Ensure target exists (0/1)
forest_clean$target <- ifelse(forest_clean$area > 0, 1, 0)

# Split
fire.split <- trainTestSplit(forest_clean, seed = 0, trainRatio = 0.8)
fire.train <- fire.split$train
fire.test  <- fire.split$test

# Build aligned design matrices; EXCLUDE area, log_area, AND target
clf_formula <- as.formula("~ . - area - log_area - target")
mm_all_clf  <- model.matrix(clf_formula, data = rbind(fire.train, fire.test))

n_tr <- nrow(fire.train)
x_all_clf   <- mm_all_clf[, -1, drop = FALSE] # drop intercept
x_train_clf <- x_all_clf[1:n_tr, , drop = FALSE]
x_test_clf  <- x_all_clf[(n_tr + 1):nrow(x_all_clf), , drop = FALSE]

y_train_clf <- fire.train$target
y_test_clf  <- fire.test$target

# Fit logistic regression
set.seed(123)
train_df_for_glm <- data.frame(target = y_train_clf, x_train_clf)
test_df_for_glm  <- data.frame(x_test_clf)

fire.logit <- glm(target ~ ., data = train_df_for_glm, family = "binomial")
summary(fire.logit)

# Predict + metrics
test_prob  <- predict(fire.logit, newdata = test_df_for_glm, type = "response")
test_label <- ifelse(test_prob > 0.5, 1, 0)

TP <- sum((test_label == 1) & (y_test_clf == 1))
FP <- sum((test_label == 1) & (y_test_clf == 0))
FN <- sum((test_label == 0) & (y_test_clf == 1))
TN <- sum((test_label == 0) & (y_test_clf == 0))

recall    <- TP / (TP + FN)
precision <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
f1        <- ifelse(is.na(precision) || (precision + recall) == 0, NA, 2*precision*recall/(precision+recall))
cat(sprintf("Recall: %.3f   Precision: %.3f   F1: %.3f\n", recall, precision, f1))

# --- Confusion matrix at the default 0.5 threshold ---
cm <- table(Pred = test_label, True = y_test_clf)
print(cm)

# --- Threshold vs Recall curve ---
th <- seq(0, 1, by = 0.01)
rec_curve <- sapply(th, function(t) {
  mean((test_prob >= t) & (y_test_clf == 1)) / mean(y_test_clf == 1)
})
plot(th, rec_curve, type = "l",
     main = "Recall vs Threshold (Fire Occurrence Prediction)",
     xlab = "Classification Threshold", ylab = "Recall")



# Optional AUC
# install.packages("pROC")
library(pROC)
auc <- roc(response = y_test_clf, predictor = as.numeric(test_prob))$auc
cat(sprintf("AUC: %.3f\n", auc))
