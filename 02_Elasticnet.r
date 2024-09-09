## Elasticnet machine learning

library(glmnet)
y <- as.factor(tumor$distant_metastasis.ch1)
X<-as.matrix(df)

library(caret)
trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[trainIndex,]
y_train <- y[trainIndex]
X_test <- X[-trainIndex,]
y_test <- y[-trainIndex]


# sequence for tuning of lambda and alpha
lambda_grid <- 10^seq(3, -2, by = -0.1)
alpha_grid <- seq(0.1, 0.9, by = 0.1)

# initialisation of variables
results <- expand.grid(alpha = alpha_grid, lambda = lambda_grid)
results$auc <- NA
library(pROC)
# loop for lambda and alpha tuning
for (i in 1:nrow(results)) {
  alpha <- results$alpha[i]
  lambda <- results$lambda[i]
  
  model <- glmnet(X_train, y_train, alpha = alpha, lambda = lambda, family = "binomial")
  
  # Prédictions on test data
  probs <- predict(model, s = lambda, newx = X_test, type = "response")
  
  # compute AUC
  roc_obj <- roc(y_test, as.vector(probs))
  auc <- auc(roc_obj)
  
  # data results
  results$auc[i] <- auc
}


# Selection of best parameters
best_params <- results[which.max(results$auc),]
best_params
best_alpha <- best_params$alpha

best_lambda <- best_params$lambda
best_auc <- best_params$auc

library(ggplot2)
# tuning visualisation
results <- results %>%
  mutate(log_lambda = log10(lambda))
library(pals)
ggplot(results, aes(x = log_lambda, y = auc, color = as.factor(alpha))) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = log10(best_lambda), linetype = "dashed", color = "red",size=1) +
	scale_color_manual(values=glasbey())+
  labs(title = "AUC depending of lambda and alpha",
       x = "log10(lambda)",
       y = "AUC",
       color = "Alpha") +
  theme_minimal()+theme(text = element_text(size = 16))



### elastic net 0.2 alpha 

fit <- glmnet(X, y, family = "binomial",alpha=0.2)
plot(fit)

cvfit <- cv.glmnet(X, y, family = "binomial",alpha=0.2)
plot(cvfit)
coefficients<-coef(fit, s = cvfit$lambda.min)
coefficients<-as.matrix(coefficients)
coefficients<-as.data.frame(coefficients)

library(dplyr)
coefficients %>% dplyr::rename(coef = "s1") %>% filter (coef != 0)->selcoef
selcoef$gene<-row.names(selcoef)
selcoef%>%arrange(desc(coef))->selcoef
selcoef%>%filter(gene != "(Intercept)")->selcoef
selcoef
selcoef%>%filter(coef>0)->pos
pos
dim(pos)

write.table(selcoef,file="selcoef.tsv",row.names=T,sep="\t")

library(ggplot2)
ggplot(data=pos,aes(x=reorder(gene,coef),y=coef))+geom_bar(stat="identity",fill="plum2")+
coord_flip()+theme_minimal()+xlab("Genesets") + ylab("Normalized enrichment score")+
geom_text(aes(label=round(coef,3)),hjust=0, vjust=0.5,color="darkblue",position= position_dodge(0),size=4,angle=0)+
xlab("Genes") + ylab("Elasticnet coefficients")+
ggtitle("") +theme(text = element_text(size = 14))+
theme(legend.position = "none")