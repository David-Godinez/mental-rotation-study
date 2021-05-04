## ---
## MRT Long RQ 2
## David Godinez
## ---


# Prelim ------------------------------------------------------------------

library(tidyverse)

library(lme4)

library(performance)

library(factoextra)

# Data Import -------------------------------------------------------------

mrt_data <- read_csv("mrt_data.csv")

mrt_data <- mrt_data %>% select(-c(X1, cr_true_bin))




# Type FIlter -------------------------------------------------------------

## Filters out FALSE trials:

mrt_data_true <- mrt_data %>% filter(CorrectResponse == TRUE)



## Includes total_change variable to model:


rt_change_model <-lmer(Stimuli.RT ~ Angle + bin_dir_x + eq_armlength + total_change + 
                         (1 + Angle + bin_dir_x + eq_armlength + total_change| Subject),
                       data = mrt_data_true,
                       REML = TRUE)


summary(rt_change_model)

check_model(rt_change_model)



acc_change_model <- glmer(Accuracy ~ Angle + bin_dir_x + eq_armlength + total_change +
                            (1 + Angle + bin_dir_x + eq_armlength + total_change | Subject),
                          data = mrt_data_true,
                          family = "binomial",
                          control=glmerControl(optimizer = "bobyqa"))


summary(acc_change_model)


write_csv(mrt_data_true, "mrt_data_match.csv")




# Scaling and K-means -----------------------------------------------------

set.seed(123)

## RT

randos <- ranef(rt_change_model)

rdf <- randos$Subject

scaled_ranef <- scale(rdf)

scaled_ranef

fviz_nbclust(scaled_ranef, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)


first_km <- kmeans(scaled_ranef, centers = 3)

first_km

fviz_cluster(first_km, data = scaled_ranef, 
             show.clust.cent = TRUE, repel = TRUE, stand = TRUE)



## ACC


randos_acc <- ranef(rt_change_model)

rdf_acc <- randos_acc$Subject

scaled_rdf_acc <- scale(rdf_acc)

fviz_nbclust(scaled_rdf_acc, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)

first_km_acc <- kmeans(scaled_rdf_acc, centers = 3)

first_km_acc

fviz_cluster(first_km_acc, data = scaled_rdf_acc, show.clust.cent = TRUE, 
             repel = TRUE, stand = TRUE)













