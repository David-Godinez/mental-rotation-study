## ---
## MRT Long ME Analysis
## David Godinez
## ---


# Prelim ------------------------------------------------------------------

library(tidyverse)

library(lme4)

library(lmerTest)

library(readxl)

library(jtools)

library(performance)


# Data Import -------------------------------------------------------------

mrt_long <- read_excel("MRT_T1_longform_cleaned.xlsx")


# EDA ---------------------------------------------------------------------

length(unique(mrt_long$Subject)) # Unique subjects = 209

length(unique(mrt_long$Trial)) # Total trials = 84

hist(mrt_long$Stimuli.RT) # Shows  RT == 0

mrt_long$Stimuli.RT[mrt_long$Stimuli.RT == 0] <- NA # Changes RT == 0 to NA


## Gives average accuracy and response time by trial type. Output is in
## descending order of `avg_rt`. NOTE: This is grouping by `Trial`, so we wont
## get individual-specific effects, but rather the average output for each trial
## across all subjects:

by_trial <- mrt_long %>%
  select(c(Trial, CorrectResponse, Angle, Direction, Accuracy, Stimuli.RT)) %>%
  group_by(Trial, CorrectResponse, Angle, Direction) %>%
  summarize(avg_accuracy = mean(Accuracy, na.rm = TRUE),
            avg_rt = mean(Stimuli.RT, na.rm = TRUE)) %>%
  arrange(desc(avg_rt))

view(by_trial)

table(by_trial$Direction) # x = 42; z = 42

table(by_trial$CorrectResponse) # 30 False and 54 True unique trials

## Shows positive linear relationship between angle difference and response time
## (i.e. response time increases with higher angle differences). This is
## consistent with Shepard & Metzler's (1971) results:

ggplot(data = by_trial) +
  stat_summary(mapping = aes(x = Angle, y = avg_rt, color = Direction)) +
  geom_smooth(aes(x = Angle, y = avg_rt, color = Direction), method = 'lm', se = FALSE) +
  xlab("Angle Difference") + ylab("Average Response Time in MS") + 
  ggtitle("Response Time by Angle Difference")


## Linear relationship is also positive when faceted by `CorrectResponse`. Angle
## diff is a better indicator of average response time for true trials; Zero
## angle difference weighs down average response times for true trials more-so
## than for false trials:

ggplot(data = by_trial) +
  stat_summary(mapping = aes(x = Angle, y = avg_rt, color = Direction)) +
  geom_smooth(aes(x = Angle, y = avg_rt, color = Direction), method = lm, se = FALSE) +
  facet_grid(. ~ CorrectResponse ) +
  xlab("Angle Difference") + ylab("Avg. Response Time in MS") + 
  ggtitle("Response Time by Angle Difference")


## Shows negative linear relationship between angle difference and average
## accuracy. For true trials, subjects tend to be less accurate when the object
## is rotated along the Z axis than when on the X axis. This is not the case for
## False trials, as they tend to be more accurate when the object is rotated
## along the Z axis:


ggplot(by_trial) +
  geom_point(mapping = aes(x = Angle, y = avg_accuracy, color = Direction)) +
  geom_smooth(aes(x = Angle, y = avg_accuracy, color = Direction), se = FALSE, method = "lm") +
  facet_grid(. ~ CorrectResponse) +
  xlab("Angle Difference") + ylab("Avg. Accuracy") + 
  ggtitle("Accuracy v. Angle, Grouped by Trial")


## Creates separate subsets of the data by correct response type:

false_trials <- by_trial %>%
  filter(CorrectResponse == "FALSE")


true_trials <- by_trial %>%
  filter(CorrectResponse == "TRUE")


# Mixed Effect Modeling ---------------------------------------------------

## Creates binary variable for x-direction trials:

mrt_long <- mrt_long %>%
  mutate(bin_dir_x = ifelse(Direction == 'x', 1, 0), .after = Direction)

mrt_long$Subject <- as.factor(mrt_long$Subject)

model_rt <- lmer(Stimuli.RT ~ Angle + bin_dir_x + CorrectResponse + (1 + Angle + bin_dir_x + CorrectResponse | Subject), 
              data = mrt_long, 
              REML = TRUE)

ranef(model_rt)



summary(model_rt)

##  Scaled residuals: 
##    Min      1Q  Median      3Q     Max 
##  -4.2018 -0.7134 -0.0854  0.6649  3.8816 
##  
##  Random effects:
##    Groups             Name        Variance    Std.Dev.  Corr             
##  Subject             (Intercept)  3.257e+05   570.710                  
##                      Angle        8.191e+00   2.862   -0.24            
##                      bin_dir_x    1.249e+03   35.346 -0.54  0.14      
##            CorrectResponseTRUE   1.843e+05  429.321 -0.22 -0.46  0.13
##  Residual                         1.286e+06 1134.168                  
##  Number of obs: 16687, groups:  Subject, 209
##  
##  Fixed effects:
##                      Estimate       Std. Error     df  t value Pr(>|t|)    
##  (Intercept)           3759.4550    45.1495  191.3015  83.267   <2e-16 ***
##    Angle                  7.1584    0.2614   211.0198  27.380   <2e-16 ***
##    bin_dir_x            -18.2298    17.7436 1282.4196  -1.027    0.304    
##  CorrectResponseTRUE   -553.9603    34.9981  192.8203 -15.828   <2e-16 ***
##    ---
##    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##  
##  Correlation of Fixed Effects:
##    (Intr) Angle  bn_dr_
##  Angle       -0.352              
##  bin_dir_x   -0.257  0.016       
##  CrrctRsTRUE -0.289 -0.326  0.013



summ(model_rt)

ranef(model_rt)




# Image characteristics ---------------------------------------------------

img_eq_arm <- c(8:11, 13:16) # Image numbers with equal arm lengths

table(mrt_long$BaseImage) # All image number counts


# Creates `eq_armlength` variable where equal image arm length == 1,
# otherwise == 0; rearranges columns for convenience:

mrt_long <- mrt_long %>% mutate(eq_armlength = ifelse(BaseImage %in% img_eq_arm, 1, 0), .after = BaseImage)

head(mrt_long)


summary(lmer(Stimuli.RT ~ Angle + bin_dir_x + CorrectResponse + eq_armlength + 
       (1 + Angle + bin_dir_x + CorrectResponse + eq_armlength | Subject), 
     data = mrt_long, 
     REML = TRUE))


## Collapsed analysis: 

eq_arm_trials <- mrt_long %>% 
  filter(eq_armlength == 1)

trial_nums <- unique(eq_arm_trials$Trial)

by_trial <- by_trial %>%
  mutate(eq_arm = ifelse(Trial %in% trial_nums, "Equal", "Not Equal"))

table(by_trial$eq_arm) # 63 Equal and 21 Not Equal


ggplot(by_trial, aes(x = Angle, y = avg_rt, color = Direction)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  facet_grid(eq_arm ~ CorrectResponse) +
  xlab("Angle Difference") + ylab("Avg. Response Time in MS") + 
  ggtitle("RT v. Angle, Grouped by Trial and Arm Length")


ggplot(by_trial, aes(x = Angle, y = avg_accuracy, color = Direction)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  facet_grid(eq_arm ~ CorrectResponse) +
  xlab("Angle Difference") + ylab("Avg. Accuracy") + 
  ggtitle("Accuracy v. Angle, Grouped by Trial and Arm Length")



# Model Troubleshooting ---------------------------------------------------


## Angle removed as RE:

model_noint <- lmer(Stimuli.RT ~ Angle + CorrectResponse + bin_dir_x + eq_armlength + 
               (1 | Subject),
     data = mrt_long,
     REML = TRUE)



## Angle as factor variable:

df <- mrt_long


df$Angle <- as.factor(df$Angle)

lmer(Stimuli.RT ~ Angle + CorrectResponse + bin_dir_x + (1 + Angle + CorrectResponse + bin_dir_x | Subject),
     data = df,
     REML = TRUE)



## Step-wise model selection:

model_inter <- lmer(Stimuli.RT ~ Angle + bin_dir_x + CorrectResponse + 
                      Angle*bin_dir_x + Angle*CorrectResponse + bin_dir_x*CorrectResponse + 
                      (1 + bin_dir_x + CorrectResponse+ 
                         Angle*bin_dir_x + Angle*CorrectResponse + bin_dir_x*CorrectResponse | Subject), 
                    data = mrt_long, 
                    REML = TRUE)

summary(model_inter)
model_step <- step(model_inter)
final <- get_model(model_step)
anova(final)




# Changes in Problem Attributes -------------------------------------------

# Creates binary variable for trial type; we need this to create `type_change`:

mrt_long <- mrt_long %>% 
  mutate(cr_true_bin = ifelse(mrt_long$CorrectResponse == "TRUE", 1, 0), .after = CorrectResponse)



# Creates binary variables for each attribute (1 = different from previous
# trial, 0 = same as previous trial):


mrt_long <- mrt_long %>%
  mutate(axis_change = ifelse(bin_dir_x != lag(bin_dir_x), 1, 0),
         armlength_change = ifelse(eq_armlength != lag(eq_armlength), 1, 0),
         type_change = ifelse(cr_true_bin != lag(cr_true_bin), 1, 0),
         angle_change = ifelse(Angle != lag(Angle), 1, 0))




# Recodes "change" variables to NA when `Trial` == 1

mrt_long$armlength_change[mrt_long$Trial == 1] <- NA
mrt_long$axis_change[mrt_long$Trial == 1] <- NA
mrt_long$type_change[mrt_long$Trial == 1] <- NA
mrt_long$angle_change[mrt_long$Trial == 1] <- NA


# Gets 'total_change' through row sums:

mrt_long <- mrt_long %>% mutate(total_change = select(., ends_with("change")) %>% rowSums)


write.csv(mrt_long, "mrt_data.csv")


