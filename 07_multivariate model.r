## multivariate model
model<-glm(metastasis~meta.score+histological_type+age_months+gender+pretext_stage,data=tumor,family="binomial")
summary(model)

library(broom)
library(broom.helpers)
library(GGally)

ggcoef_model(model,exponentiate=T,colour_guide=TRUE)+
scale_color_brewer(palette="Set2")+theme(text=element_text(size=18))+theme(legend.position="none")


# Coefficient for the predictor
coef_predictor <- coef(model)["meta.score"]

# Odds Ratio
odds_ratio <- exp(coef_predictor)
odds_ratio

library(regplot)
regplot(model, clickable=F, points=T, droplines=T,rank="sd") 