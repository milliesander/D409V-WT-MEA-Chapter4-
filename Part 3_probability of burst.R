# D409V/WT medial septum MEA analysis 
# Author: Millie Sander
# Part 3: Probability of burst 

# loading required packages
library(xlsx)
library(lme4)
library(emmeans)
library(car)
library(fitdistrplus)

p_burst_HQ_NAremoved <- read_excel("D:/Thesis/MEA/p_burst_HQ_NAremoved.xlsx")

data <- p_burst_HQ_NAremoved 

data$Genotype <- as.factor(data$Genotype)
data$Recording <- as.factor(data$Recording)
data$Cell_ID <- as.factor(data$Cell_ID)
data$Animal_ID <- as.factor(data$Animal_ID)
data$Condition <- as.factor(data$Condition)
data$Sex <- as.factor(data$Sex)
data$pb_bl <- as.numeric(data$pb_bl)

m1 <- lmer(pb_bl ~ Genotype*Condition + (1|Animal_ID/Recording/Cell_ID), data = data)
m2 <- lmer(pb_bl ~ Genotype*Condition + Sex + (1|Animal_ID/Recording/Cell_ID), data = data)

# seeing if there is a sig diff between models, and which is the best fit 
anova(m1, m2)

# testing assumptions of chosen model (m1)
qqnorm(residuals(m1))
qqline(residuals(m1))
hist(residuals(m1))
plot(m1)

vif(m1)

leveneTest(residuals(m1)~data$Genotype)
shapiro.test(residuals(m1))

# checking theoretical distribution 
library(fitdistrplus)
descdist(data$pb_bl)

# beta glmer using scaled Hz variable and glmmTMB package 

data$pb_bl[data$pb_bl == "0"] <- 0.00000001
data$pb_bl[data$pb_bl == "1"] <- 0.99

m1b <- glmmTMB(pb_bl ~ Genotype*Condition + (1|Animal_ID/Recording/Cell_ID), data = data, family = beta_family())

joint_tests(m1b)

emmeans(m1b, list(pairwise ~ Genotype*Condition), adjust = "bonferroni")
 
